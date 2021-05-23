function [sumx_grace]=grace2gravity(Lmax,k,Nlmx,Detagc,Detags,fir,n_c,n_f,nceta,nfir)

%--------
a=6.378136460E+06;% 单位m
GM=3.98603E+14;%单位m3/s2

mfir=zeros(Lmax+1,n_f);
for j=1:n_f
    for m=0:Lmax
        mfir(m+1,j)=m*fir(j);
    end
end

cosdmf=cosd(mfir);
sindmf=sind(mfir);
%************************************
clear  fir mfir;


sumg_grace=zeros(n_c,n_f,k);
sumx_grace=zeros(n_c*n_f,k);


n=0:Lmax;
kk=n-1;

for i=1:k
    for nn=1:n_c
        for j=1:n_f
            sumg_grace(nn,j,i)=GM/(a^2)*kk*(Nlmx(:,:,nn).*Detagc(:,:,i)*cosdmf(:,j)+Nlmx(:,:,nn).*Detags(:,:,i)*sindmf(:,j))*1E08;
        end
    end
end

for i=1:n_c*n_f
    sumx_grace(i,:)=sumg_grace(nceta(i),nfir(i),:);
end