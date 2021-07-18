%% Script to calculate the phase-synchronization measures corrected with the surrogates
%This script calculates phase synchronization measures to study brain
%dynamics. The goal here is to describe the synchornization dynamics based
%on phase statistics. The present script computes four basic measures;
%integration, segregation, phase-interaction fluctuations and FCD.
%These measures are applied mainly in fMRI BOLD signals.

%BASIC PHASE SYNCHRONIZATION MEASURES
%------------------------------
%Integration
%    
%Segregation
%     
%Phase-interaction fluctuations
%    
%FCD
%    
%

%Ane López-González 25/01/2021
clc
close all
clear all

%% Loading data and basic parameters:

load('data_UWS_all1') %timeseries of each group Nsubjects*nodes*time
ts=ts_all;
NSUB = size(ts,1);
N = size(ts,2);

%Basic filtering parameters
delta= 2; %TR
flp = .04;     % lowpass frequency of filter
fhi = .07;    % highpass
k = 2;                  % 2nd order butterworth filter
fnq = 1/(2*delta);       % Nyquist frequency
Wn = [flp/fnq fhi/fnq]; % butterworth bandpass non-dimensional frequency
[bfilt2,afilt2] = butter(k,Wn);   % construct the filter



%% Surrogate analysis 
%account for the phase-effects that can be removed by the surrogates.
%In this part of the code, we create surrogates of the BOLD timeseries and
%from them, the mean surrogate Phase-Interaction matrix, that will be used
%in the following section to clean the empirical matrices.


for nsub = 1:NSUB
xs = squeeze(ts(nsub,:,:));
Tmax = size(xs,2); %%%The time is calculated inside the subjects loop 
                   %%%in case any subject has different time duration.
T = 1:Tmax;
timeseriedata = zeros(N,Tmax);%length(xs)

for surr=1:100 %we used 100 surrogates as default
    
 for seed = 1:N
      x = demean(detrend(xs(seed,:)));
      x_surr=surrogates(x); % Create the phase-randomized copy of x
      timeseriedata(seed,:) = filtfilt(bfilt2,afilt2,x_surr);    % zero phase filter the data
      Xanalytic = hilbert(demean(timeseriedata(seed,:))); %Hilbert transformation
      Phases(seed,:) = angle(Xanalytic); %%% calculating phase of each ROI for each signal
                                         %%% which will use to compute the
                                         %%% phase synchronization measures
                                       
 end
 
% Phase-interaction matrix: At each time point, the phase difference between two regions was calculated: 

 for t = T
  for i = 1:N
    for j = 1:N
     dM_surr(i,j,t) = cos(adif(Phases(i,t),Phases(j,t))); 
    end
  end  
   %mPLM(nsub,t)=mean(mean(dM_surr(:,:,t))); 
 end
mdM_surr=mean(dM_surr,3); 
mdM(surr,:,:)=mdM_surr; %Mean phase-interaction matrix across time
surr
end

save(sprintf('surrogate_UWS_%02d',nsub),'mdM') %save the surrogate analysis
end

clear Phases x xs timeseriedata Xanalytic

%% Calculate the phase-synchronization measured from the
%% real BOLD timeseries and then corrected with the surrogates


for nsub = 1:NSUB
xs = squeeze(ts(nsub,:,:));
Tmax = size(xs,2); %%%The time is calculated inside the subjects loop 
                   %%%in case any subject has different time duration.
T = 1:Tmax;
timeseriedata = zeros(N,Tmax);%length(xs)
for seed = 1:N
  x = demean(detrend(xs(seed,:)));
  timeseriedata(seed,:) = filtfilt(bfilt2,afilt2,x);    % zero phase filter the data
  Xanalytic = hilbert(demean(timeseriedata(seed,:))); %Hilbert transformation
  Phases(seed,:) = angle(Xanalytic); %%% calculating phase of each ROI for each signal
                                     %%% which will use to compute the
                                     %%% phase synchronization measures
                                       
end

T = 1:Tmax;
sync = zeros(1, Tmax);
for t = T
    ku = sum(complex(cos(Phases(:,t)),sin(Phases(:,t))))/N;
    sync(t) = abs(ku);
end

fluct(nsub) = std(sync(:)); 
sync_all (nsub,:) =sync; 
%Phase-interaction matrix: 
%At each time point, the phase difference between two regions was calculated: 
 for t = T
  for i = 1:N
    for j = 1:N
     dM(i,j,t) = cos(adif(Phases(i,t),Phases(j,t)));  % computes dynamic matrix/ dFC
    end
  end 

 end
            
fluct_data(nsub,:)=std(sync);
   
load(sprintf('surrogate_UWS_%02d',nsub))
mmdM=squeeze(mean(mdM,1));

%Integration: 

cc = mean(dM,3)-mmdM; % Correct with the mean matrices calculated with the surrogates
%cc = mean(dM,3);
cc = cc-eye(N);
pp = 1;
PR = 0:0.01:0.99;
cs=zeros(1,length(PR));

for p = PR
A = (cc)>p;
[~, csize] = get_components(A);
cs(pp) = max(csize);
pp = pp+1;
end

integ(nsub) = sum(cs)*0.01/N;

% The modularity (as a measure of segregation) is calculated in the mean matrix and corrected with the 
% bined matrix given by the surrogate and imposing a threhsold of the
% 99% percentile

meandM=mean(dM,3);
for i=1:N
  for j=1:N
      Y=prctile(mdM(:,i,j),99);
      if meandM(i,j)>Y
      bin_M(i,j)=1;
      else 
      bin_M(i,j)=0;
      end
  end
end

[~, Q(nsub)] = community_louvain((bin_M));


%FCD (phases)
Isubdiag = find(tril(ones(N),-1));
for t=T
    for i=1:N
    for j=1:N
     dM(i,j,t) = cos(adif(Phases(i,t),Phases(j,t)));  % computes dynamic matrix/ dFC
    end
    end 
    patt=dM(:,:,t);
    pattern(t,:)=patt(Isubdiag);
end

npattmax=size(pattern,1);
kk3=1;
for t=1:npattmax-30
    p1=mean(pattern(t:t+30,:));
    for t2=t+1:npattmax-30
        p2=mean(pattern(t2:t2+30,:));
        phfcddata(kk3)=dot(p1,p2)/norm(p1)/norm(p2);
        matrixcdc(t,t2)=dot(p1,p2)/norm(p1)/norm(p2);
        matrixcdc(t2,t)=dot(p1,p2)/norm(p1)/norm(p2);
        kk3=kk3+1;
    end
end 
CDC(nsub,:)=phfcddata;
matrix_CDC(nsub,:,:)=matrixcdc;

display(nsub);
end

save('phase_synchronization_UWS','integ','Q','fluct','CDC','matrix_CDC')