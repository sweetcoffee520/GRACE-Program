% function [vary_f,vary_s,vary_fa,vary_sa,vary_fas,vary_chazhix]=gravityx(x,y,DeltaGC,DeltaGS,Lmax,Nlmx,k,loveN_k1)
Lmax=60;
[x,y,c,d,e,f,g,h]=region_grid(-90,90,-180,180,1);
[DeltaGC,DeltaGS,k]=read_CRCOF2_CSR('I:\86442\Level-2\GRACE-FO\Level-2\CSR\RL06\',60);
Nlmx = Nlmx_v3(60,x);
Filename='D:\code\GRACE_py\data\loadLove.txt';
loveN_k1=-textread(Filename,'%n');
a=6371004;  %%地球平均半径m
pave=5507.85;         %%地球平均密度kg/m^3
xx=length(x);
yy=length(y);

w1=zeros(1,Lmax+1);
radius=300;
% a=6.378136460E+06;% 单位m
r1=radius*1000;%%Unit:km
b1=log(2)/(1-cos(r1/a));
w1(1)=1;
w1(2)=(1+exp(-2*b1))/(1-exp(-2*b1))-1/b1;
for l=1:(Lmax-1)
    w1(l+2)=-(2*l+1)/b1*w1(l+1)+w1(l);
end

r2=radius*1000;
w2=zeros(1,Lmax+1);
b2=log(2)/(1-cos(r2/a));
w2(1)=1;
w2(2)=(1+exp(-2*b2))/(1-exp(-2*b2))-1/b2;
for l=1:(Lmax-1)
w2(l+2)=-(2*l+1)/b2*w2(l+1)+w2(l);
end

w3=zeros(1,Lmax+1);
radius=400;
% a=6.378136460E+06;% 单位m
r1=radius*1000;%%Unit:km
b1=log(2)/(1-cos(r1/a));
w3(1)=1;
w3(2)=(1+exp(-2*b1))/(1-exp(-2*b1))-1/b1;
for l=1:(Lmax-1)
    w3(l+2)=-(2*l+1)/b1*w3(l+1)+w3(l);
end

r2=radius*1000;
w4=zeros(1,Lmax+1);
b2=log(2)/(1-cos(r2/a));
w4(1)=1;
w4(2)=(1+exp(-2*b2))/(1-exp(-2*b2))-1/b2;
for l=1:(Lmax-1)
w4(l+2)=-(2*l+1)/b2*w4(l+1)+w4(l);
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
% sumx_grace=zeros(xx*yy,k);
sumg_gracex=zeros(xx,yy,k);
sumg_gracexx=zeros(xx,yy,k);
% sumx_gracex=zeros(xx*yy,k);
sumg_gracey=zeros(xx,yy,k);
sumg_graceyy=zeros(xx,yy,k);
% 未处理
for i=1:k
    for nn=1:xx
        %         for j=1:yy
        %             sumg_grace(nn,j,i)=a/3*pave*loveN*(Nlmx(:,:,nn).*DeltaGC(:,:,i)*cosdmf(:,j)+Nlmx(:,:,nn).*DeltaGS(:,:,i)*sindmf(:,j))/10;
        sumg_grace(nn,:,i)=a/3*pave*loveN*(Nlmx(:,:,nn).*DeltaGC(:,:,i)*cosdmf+Nlmx(:,:,nn).*DeltaGS(:,:,i)*sindmf)/10;
    end
end

% 高斯滤波300km
for i=1:k
    for nn=1:xx
%         for j=1:yy
            sumg_gracex(nn,:,i)=2*a/3*pave*loveN.*w1*(Nlmx(:,:,nn).*DeltaGC(:,:,i)*cosdmf+Nlmx(:,:,nn).*DeltaGS(:,:,i)*sindmf)/10;
    end
end

% Fa滤波300km
% for i=1:k
%     for nn=1:xx
%         for j=1:yy
%             sumg_gracey(nn,j,i)=a*pave/3*w2.*loveN*(Nlmx(:,:,nn).*DeltaGC(:,:,i)*(cosdmf(:,j).*w2')+Nlmx(:,:,nn).*DeltaGS(:,:,i)*(sindmf(:,j).*w2'))/10;
%         end
%     end
% end
% %高斯滤波400km
% for i=1:k
%     for nn=1:xx
%         for j=1:yy
%             sumg_gracexx(nn,j,i)=2*a/3*pave*loveN.*w3*(Nlmx(:,:,nn).*DeltaGC(:,:,i)*cosdmf(:,j)+Nlmx(:,:,nn).*DeltaGS(:,:,i)*sindmf(:,j))/10;
%         end
%     end
% end
% % Fa滤波400km
% for i=1:k
%     for nn=1:xx
%         for j=1:yy
%             sumg_graceyy(nn,j,i)=a*pave/3*w4.*loveN*(Nlmx(:,:,nn).*DeltaGC(:,:,i)*(cosdmf(:,j).*w4')+Nlmx(:,:,nn).*DeltaGS(:,:,i)*(sindmf(:,j).*w2'))/10;
%         end
%     end
% end


% lat=ones(yy,xx);
% for i=1:xx
%     lat(:,i)=lat(:,i)*x(i);
% end
% lon=zeros(yy,xx);
% for j=1:xx
%     lon(:,j)=lon(:,j)+y';
% end
% 
% latx=reshape(lat,xx*yy,1);
% lonx=reshape(lon,xx*yy,1);
% vary_f3=zeros(xx*yy,k,'double');
% vary_f=zeros(xx*yy,3,k,'double');
% vary_s3=zeros(xx*yy,k,'double');
% vary_s=zeros(xx*yy,3,k,'double');
% vary_fa3=zeros(xx*yy,k,'double');
% vary_fa=zeros(xx*yy,3,k,'double');
% vary_sa3=zeros(xx*yy,k,'double');
% vary_sa=zeros(xx*yy,3,k,'double');
% vary_fas3=zeros(xx*yy,k,'double');
% vary_fas=zeros(xx*yy,3,k,'double');
% % vary_chazhi3=zeros(xx*yy,k,'double');
% % vary_chazhi=zeros(xx*yy,3,k,'double');
% vary_chazhix3=zeros(xx*yy,k,'double');
% vary_chazhix=zeros(xx*yy,3,k,'double');
% for r=1:k
%     vary_f2=sumg_grace(:,:,r);
%     vary_f1=vary_f2';
%     vary_f3(:,r)=reshape(vary_f1,xx*yy,1);
%     vary_f(:,:,r)=[lonx,latx,vary_f3(:,r)];
%     vary_s2=sumg_gracex(:,:,r);
%     vary_s1=vary_s2';
%     vary_s3(:,r)=reshape(vary_s1,xx*yy,1);
%     vary_s(:,:,r)=[lonx,latx,vary_s3(:,r)];
%     vary_fa2=sumg_gracey(:,:,r);
%     vary_fa1=vary_fa2';
%     vary_fa3(:,r)=reshape(vary_fa1,xx*yy,1);
%     vary_fa(:,:,r)=[lonx,latx,vary_fa3(:,r)];
%     vary_sa2=sumg_gracexx(:,:,r);
%     vary_sa1=vary_sa2';
%     vary_sa3(:,r)=reshape(vary_sa1,xx*yy,1);
%     vary_sa(:,:,r)=[lonx,latx,vary_sa3(:,r)];
%     vary_fas2=sumg_graceyy(:,:,r);
%     vary_fas1=vary_fas2';
%     vary_fas3(:,r)=reshape(vary_fas1,xx*yy,1);
%     vary_fas(:,:,r)=[lonx,latx,vary_fas3(:,r)];
% %     vary_chazhi2=vary_fa2-vary_s2;
% %     vary_chazhi1=vary_chazhi2';
% %     vary_chazhi3(:,r)=reshape(vary_chazhi1,xx*yy,1);
% %     vary_chazhi(:,:,r)=[lonx,latx,vary_chazhi3(:,r)];
%     vary_chazhix2=vary_fa2-vary_s2;
%     vary_chazhix1=vary_chazhix2';
%     vary_chazhix3(:,r)=reshape(vary_chazhix1,xx*yy,1);
%     vary_chazhix(:,:,r)=[lonx,latx,vary_chazhix3(:,r)];
% 
% 
%     %     sta1=xlswrite('database.xls',vary_f4(:,:,r),output);
%     %     sta2=xlswrite('databasex.xls',vary_s4(:,:,r),output);
% end
% % for i=1:xx*yy
% %     sumx_grace(i,:)=sumg_grace(g(i),h(i),:);
% %     sumx_gracex(i,:)=sumg_gracex(g(i),h(i),:);
% % end