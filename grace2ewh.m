function [sumx_grace]=grace2ewh(Lmax,k,Nlmx,D_Detagc,D_Detags,fir,n_c,n_f,nceta,nfir)

loveN_k0=[0,0.027,-0.303,-0.194,-0.132,-0.104,-0.089,-0.081,-0.076,-0.072,-0.069,-0.064,-0.058,-0.051,-0.040,-0.033,-0.027,-0.020,-0.014,-0.010,-0.017];
n_loveN=[0,1,2,3,4,5,6,7,8,9,10,12,15,20,30,40,50,70,100,150,200];

n=0:Lmax;
loveN_k=-interp1(n_loveN,loveN_k0,n);

loveN=(2*n+1)./(1+loveN_k);

%--------
a=6.378136460E+06;% µ•Œªm
pave=5507; %√‹∂»  

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
for i=1:k
    for nn=1:n_c
        for j=1:n_f
            sumg_grace(nn,j,i)=a/3*pave*loveN*(Nlmx(:,:,nn).*D_Detagc(:,:,i)*cosdmf(:,j)+Nlmx(:,:,nn).*D_Detags(:,:,i)*sindmf(:,j))/10;
        end
    end
end

for i=1:n_c*n_f
    sumx_grace(i,:)=sumg_grace(nceta(i),nfir(i),:);
end