function [ave_vary_ax,vary_ax4,ave_vary_axx,vary_axx4]=calculate_ax(x,y,DeltaGCX,DeltaGSX,DeltaGCXX,DeltaGSXX,Lmax,Nlmx,k,loveN_k1)

a=6371004;  %%地球平均半径m
pave=5507.85;         %%地球平均密度kg/m^3
xx=length(x);
yy=length(y);

w1=zeros(1,Lmax+1);
radius=400;
% a=6.378136460E+06;% 单位m
r1=radius*1000;%%Unit:km
b1=log(2)/(1-cos(r1/a));
w1(1)=1;
w1(2)=(1+exp(-2*b1))/(1-exp(-2*b1))-1/b1;
for l=1:(Lmax-1)
    w1(l+2)=-(2*l+1)/b1*w1(l+1)+w1(l);
end

loveN_k=loveN_k1(1:Lmax+1,1)';
n=0:Lmax;
% n1=2*n;
% loveN_kk=1+loveN_k;
loveN=(2*n+1)./(1+loveN_k);

mfir=zeros(Lmax+1,yy);
for j=1:yy
    for m=0:Lmax
        mfir(m+1,j)=m*y(j);
    end
end

cosdmf=cosd(mfir);
sindmf=sind(mfir);
%************************************
clear  fir mfir;

sumg_grace=zeros(xx,yy,k);
sumg_gracex=zeros(xx,yy,k);

result=importdata('C:\GLDAS\water\datafile\GIA_n100_uplift_0km.txt',' ',18);
B=result.data(:,3);

for i=1:yy
    for j=1:xx
        C(i,j)=B((y(i)+0.5)*180+x(j)+0.5+90);
    end
end
% C=reshape(B,xx,yy);

%高斯滤波400km
for i=1:k
    for nn=1:xx
        for j=1:yy
            sumg_grace(nn,j,i)=2*a/3*pave*loveN.*w1*(Nlmx(:,:,nn).*DeltaGCX(:,:,i)*cosdmf(:,j)+Nlmx(:,:,nn).*DeltaGSX(:,:,i)*sindmf(:,j))/10;
        end
    end
    vary_ax(:,:,i)=sumg_grace(:,:,i)-C'./120;
end

for i=1:k
    for nn=1:xx
        for j=1:yy
            sumg_gracex(nn,j,i)=2*a/3*pave*loveN.*w1*(Nlmx(:,:,nn).*DeltaGCXX(:,:,i)*cosdmf(:,j)+Nlmx(:,:,nn).*DeltaGSXX(:,:,i)*sindmf(:,j))/10;
        end
    end
    vary_axx(:,:,i)=sumg_gracex(:,:,i)-C'./120;
end



lat=ones(yy,xx);
for i=1:xx
    lat(:,i)=lat(:,i)*x(i);
end
lon=zeros(yy,xx);
for j=1:xx
    lon(:,j)=lon(:,j)+y';
end
latx=reshape(lat,xx*yy,1);
lonx=reshape(lon,xx*yy,1);
vary_ax3=zeros(xx*yy,k,'double');
vary_ax4=zeros(xx*yy,3,k,'double');
vary_axx3=zeros(xx*yy,k,'double');
vary_axx4=zeros(xx*yy,3,k,'double');
for r=1:k
%     vary_ax2=sumg_grace(:,:,r);
    vary_ax2=vary_ax(:,:,r);
    vary_ax1=vary_ax2';
    vary_ax3(:,r)=reshape(vary_ax1,xx*yy,1);
    vary_ax4(:,:,r)=[lonx,latx,vary_ax3(:,r)];
%     vary_axx2=sumg_gracex(:,:,r);
    vary_axx2=vary_axx(:,:,r);
    vary_axx1=vary_axx2';
    vary_axx3(:,r)=reshape(vary_axx1,xx*yy,1);
    vary_axx4(:,:,r)=[lonx,latx,vary_axx3(:,r)];
end


lat=ones(yy,xx);
for i=1:xx
    lat(:,i)=lat(:,i)*x(i);
end

latxx=cosd(lat);
latx=reshape(latxx,1,xx*yy);

for i=1:k
ave_vary_ax(i)=latx*reshape(sumg_grace(:,:,i)',xx*yy,1)./sum(latx(:));
ave_vary_axx(i)=latx*reshape(sumg_gracex(:,:,i)',xx*yy,1)./sum(latx(:));
end
