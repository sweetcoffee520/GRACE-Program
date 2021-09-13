function vary_jy=water(x,y,c,d,Lmax,k,Nlmx,loveN_k1,dc,ds)
% Lmax=60;
% k=12;
% Filename='E:\GLDAS\water\datafile\Love_n.txt';
% loveN_k1=-textread(Filename,'%n');
% griddata=read_GLDAS_nc4('E:\GLDAS\2水文数据\');
% [dc,ds]=calculate_CS(Lmax,loveN_k1,1,griddata,k);
a=6371004;  %%地球平均半径m
pave=5507.85;         %%地球平均密度kg/m^3
% c
loveN_k=loveN_k1(1:Lmax+1,1)';
n=0:Lmax;
% n1=2*n;
% loveN_kk=1+loveN_k;
loveN=(2*n+1)./(1+loveN_k);

mfir=zeros(Lmax+1,d);
for j=1:d
    for m=0:Lmax
        mfir(m+1,j)=m*y(j);
    end
end
sumg_gracex=zeros(c,d,k);
cosdmf=cosd(mfir);
sindmf=sind(mfir);
for i=1:k
    for nn=1:c
        for j=1:d
            sumg_gracex(nn,j,i)=a/3*pave*loveN*(Nlmx(:,:,nn).*dc(:,:,i)*cosdmf(:,j)+Nlmx(:,:,nn).*ds(:,:,i)*sindmf(:,j))/10;
        end
    end
end
% temp=sumg_gracex(:,181:360,:);
% average_num=importdata('E:\GLDAS\water\datafile\GIA_n100_mass_0km.txt',' ',18);
% B=average_num.data/120;
% C=reshape(B(:,3),[c,d]);
% for i=1:k
% sumg_gracex(:,:,i)=sumg_gracex(:,:,i)-C;
% end
% % % 如是1-360则需替换
% sumg_grace(:,181:360,:)=sumg_gracex(:,1:180,:);
% sumg_grace(:,1:180,:)=sumg_gracex(:,181:360,:);
lat=ones(d,c);
for i=1:c
    lat(:,i)=lat(:,i)*x(i);
end
lon=zeros(d,c);
% j1=-179.5:1:179.5;
for j=1:c
    lon(:,j)=lon(:,j)+y';
end
latx=reshape(lat,c*d,1);
lonx=reshape(lon,c*d,1);
vary_jy3=zeros(c*d,k,'double');
vary_jy=zeros(c*d,3,k,'double');
for r=1:k
    vary_jy2=sumg_gracex(:,:,r);
    vary_jy1=vary_jy2';
    vary_jy3(:,r)=reshape(vary_jy1,c*d,1);
%     vary_jy3(:,r)=vary_jy3(:,r)-B(:,3);
    vary_jy(:,:,r)=[lonx,latx,vary_jy3(:,r)];
end