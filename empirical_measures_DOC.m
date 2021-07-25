%% Script to calculate the phase-synchronization measures corrected with the surrogates
%This script calculates phase synchronization measures to study brain
%dynamics. The goal here is to describe the synchornization dynamics based
%on phase statistics. The present script computes four basic measures;
%integration, segregation, phase-interaction fluctuations and FCD.
%These measures are applied mainly in fMRI BOLD signals.

%Ane Lopez-Gonzalez (18-7-2021)

%BASIC PHASE SYNCHRONIZATION MEASURES
%See paper:'Loss of consciousness reduces the stability of brain hubs and
%the heterogeneity of brain dynamics' by Lopez-Gonzalez, A and Panda, R et al. 
%for detailed information
%------------------------------
%Integration  
%Segregation    
%Phase-interaction fluctuations
%FCD


clc
close all
clear all



%% Calculate the phase-synchronization measured from the
%real BOLD timeseries and then corrected with the surrogates
N=214;%number of regions
T=295;%time-points


for group=1:3
    if group==1
       SUBJ=1:35;
       n_group1='CNT'
       n_group2='controls'
    elseif group==2
       SUBJ=[1:8 10:13 15:18 20:31 33:37]
       n_group1='MCS'
       n_group2='MCS'
    elseif group==3
       SUBJ=[1:7 11:18];
       n_group1='UWS'
       n_group2='UWS'
       
    end
    nsub1=1;
    mPLM=[];
    meta=[];
    integ=[];
    Q=[];
    matrix_CDC=[];
    CDC=[];
     for nsub = SUBJ
        load(sprintf('/home/ane/Escritorio/DOC_anesthesia/Figura_2/new_Adrian/phase_interaction_matrix_time_%s_subj_%02d',n_group1, nsub))
        meandM=squeeze(mean(dM,3));
        load(sprintf('/home/ane/Escritorio/DOC_anesthesia/Figura_2/new_Adrian/surrogates_results/m_surrogate_%s_%02d',n_group2, nsub))
        mmdM=squeeze(mean(mdM,1));
        
        %% Phase-interaction fluctuations
         for t = 1:T
          mPLM(nsub1,t)=mean(mean(dM(:,:,t)));
         end
        meta(nsub1)=std(mPLM(nsub1,:));

        %% Integration: 

        cc = meandM-mmdM; % Correct with the mean matrices calculated with the surrogates
        %cc = mean(dM,3);
        cc = cc-eye(N);
        pp = 1;
        PR = 0:0.01:0.99;
        cs=zeros(1,length(PR));

        for p = PR
        A = (cc)>p;
        [~, csize] = get_components(double(A));
        cs(pp) = max(csize);
        pp = pp+1;
        end

        integ(nsub1) = sum(cs)*0.01/N;

        %% The modularity (as a measure of segregation) is calculated in the mean matrix and corrected with the 
        % bined matrix given by the surrogate and imposing a threhsold of the
        % 99% percentile

        meandM=squeeze(mean(dM,3));
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

        [~, Q(nsub1)] = community_louvain((bin_M));
        
        %% FCD (phases)
        %this part takes longer than the rest, so if you are doing a quick
        %test you can comment it
        Isubdiag = find(tril(ones(N,N),-1));
        for t=1:T
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
        CDC(nsub1,:)=phfcddata;
        matrix_CDC(nsub1,:,:)=matrixcdc;
    nsub1=1+nsub1;
     end
save(sprintf('phase_synchronization_results_%s',n_group1 ),'integ','Q','meta','CDC')    
end


