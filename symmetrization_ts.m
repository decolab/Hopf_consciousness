function totalH2=symmetrization_ts(ts_all,Nodes,NSUB)
%This function orders the timeseries so it is oredered by ROIs and no by
%hemispheres.
nn=1;
mid=Nodes/2;
for n=1:NSUB
    k=1;
    for n2=1:mid
        ts_all1(n,k,:)=ts_all(n,n2,:);
	k=k+1;
       ts_all1(n,k,:)=ts_all(n,n2+mid,:);
        k=k+1;
    end
   totalH2{nn}=squeeze(ts_all1(n,:,:)) ;
   nn=nn+1;
end
