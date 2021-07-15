%% Script for the calculation of the surropgates
%that we will use to correct the matrices of the itnegration and segregation

%Ane López-González

clear all
load('data_UWS_all') %timeseries of each group Nsubjects*nodes*time
ts=ts_all;
NSUB = size(ts,1);
    
N = size(ts,2);

%Basic filtering parameters
%%%%%%%%%%%%%
delta=2;
flp = .04;     % lowpass frequency of filter
fhi = .07;    % highpass
k = 2;                  % 2nd order butterworth filter
fnq = 1/(2*delta);       % Nyquist frequency
Wn = [flp/fnq fhi/fnq]; % butterworth bandpass non-dimensional frequency
[bfilt2,afilt2] = butter(k,Wn);   % construct the filter


for nsub = 1:NSUB
    clear timeseriedata events Phases bpl bplth;
    xs = squeeze(ts(nsub,:,:));
    Tmax = size(xs,2); 
    T = 1:Tmax;
    timeseriedata = zeros(N,Tmax);%length(xs)
    
for surr=1:100 %we used 100 surrogates (this value depends on the number of subjects)
    
 for seed = 1:N
      x = demean(detrend(xs(seed,:)));
      x_surr=surrogates(x);
      timeseriedata(seed,:) = filtfilt(bfilt2,afilt2,x_surr);    % zero phase filter the data
      Xanalytic = hilbert(demean(timeseriedata(seed,:)));
      Phases(seed,:) = angle(Xanalytic); %%% calculating phase of each ROI for each signal
                                         %%% which will use to compute
                                         %%% metastability and other
                                         %%% parameters
 end
 
% Phase-interaction matrix 

 for t = T
  for i = 1:N
    for j = 1:N
     dM_surr(i,j,t) = cos(adif(Phases(i,t),Phases(j,t)));  % computes dynamic matrix/ dFC
     %dM(i,j,t)=cos(abs(Phases(i,t)-Phases(j,t))); %dM/dFC/phasematrix all same 
     %phasematrix(i,j)=abs(x(i,t)-x(j,t));
    end
  end 
  
   mPLM(nsub,t)=mean(mean(dM(:,:,t))); %mean phase-interaction matrix across time

 end
 
mdM_surr=mean(dM_surr,3);
mdM(surr,:,:)=mdM_surr;

end
save(sprintf('m_surrogate_UWS_%02d',subj),'mdM')
end




%% Calculate the integration and correct with the surrogates


for nsub = 1:NSUB
    clear timeseriedata events Phases bpl bplth;
    xs = squeeze(ts(nsub,:,:));
    Tmax = size(xs,2); 
    T = 1:Tmax;
    timeseriedata = zeros(N,Tmax);%length(xs)
    for seed = 1:N
      x = demean(detrend(xs(seed,:)));
      timeseriedata(seed,:) = filtfilt(bfilt2,afilt2,x);    % zero phase filter the data
      Xanalytic = hilbert(demean(timeseriedata(seed,:)));
      Phases(seed,:) = angle(Xanalytic); %%% calculating phase of each ROI for each signal
                                         %%% which will use to compute
                                         %%% metastability and other
                                         %%% parameters
 end


%Phase-interaction matrix: 
 for t = T
  for i = 1:N
    for j = 1:N
     dM(i,j,t) = cos(adif(Phases(i,t),Phases(j,t)));  % computes dynamic matrix/ dFC
     %dM(i,j,t)=cos(abs(Phases(i,t)-Phases(j,t))); %dM/dFC/phasematrix all same 
     %phasematrix(i,j)=abs(x(i,t)-x(j,t));
    end
  end 
 end
 
load(sprintf('m_surrogate_W1_%02d',nsub))
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

% The segregation is calculated in the mean matrix and corrected with the 
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
display(nsub);
end

save('integration_segregation_surrogates','integ','Q')