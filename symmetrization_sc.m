function C=symmetrization_sc(C_old,Nodes)
%This function orders the SC so it is oredered by ROIs and no by
%hemispheres.
k=1;
mid=Nodes/2;
C=zeros(Nodes,Nodes);
for n2=1:mid
    bb=C_old(n2,:);
    k2=1;
    for n3=1:mid
        aa(k2)=bb(n3);
        k2=k2+1;
        aa(k2)=bb(n3+mid);
        k2=k2+1;

    end
    C(k,:)=aa;
    C(:,k)=aa';
    
    k=k+1;
    bb=C_old(n2+mid,:);
    k2=1;
    for n3=1:mid
        aa(k2)=bb(n3);
        k2=k2+1;
        aa(k2)=bb(n3+mid);
        k2=k2+1;
    end
    C(k,:)=aa;
    C(:,k)=aa';
    k=k+1;
end