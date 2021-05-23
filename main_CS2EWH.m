clc;clear;
Output_data_address='C:\Results\Tibet\';%保存结果路径
% load('CS_model.mat');
% load('CS_CSR_mascon_GIA.mat');
% load('CS_original_163.mat');
% load('CS_test.mat');
% load('model_time.mat');
% GC=GC_model_140;%输入模型值
% GS=GS_model_140;
% load('CS_original_d96.mat');
% load('CS_GAD_158.mat');
% load('CS_p3m6_fan_140.mat');
% Data_time_select=Data_time_163;
% Data_time_select=[2005 12];
% Data_time_select=model_time;
% GC=GC_original_158(:,:,42)-mean(GC_original_158,3);%扣除所有月份的平均值
% GS=GS_original_158(:,:,42)-mean(GS_original_158,3);
% GC=GC_original_d96(:,:,42)-mean(GC_original_d96,3);%扣除所有月份的平均值
% GS=GS_original_d96(:,:,42)-mean(GS_original_d96,3);
% GC=GC_CSR_mascon_163;%输入模型值
% GS=GS_CSR_mascon_163;
% GC(1,1,:)=1;GC(2,1,:)=0;
% GC(2,2,:)=0;GS(2,2,:)=0;
% GC=GC_original_163(:,:,1:140);%输入原始数据
% GS=GS_original_163(:,:,1:140);
% GC=GC_GAD_158;%输入模型值
% GS=GS_GAD_158;
% Data_time_select=Data_time;
% GC=GC_p3m6_fan_140;%输入p3m6的结果
% GS=GS_p3m6_fan_140;
% Data_time_select=time_p3m6;
% GC=GC_original;=GS_original;

[~,Lmax,Data_number] = size(GC);
Lmax=Lmax-1;

%下面引入一阶项和C20项数据进行替换
% load('CS_deg1_all.mat');
load('C20_coef_all.mat');
% 挑选出所选月份的数据
% start_j=1;
% for i=1:Data_number
%     for j=start_j:length(CS_deg1_time(:,1))
%         if Data_time_select(i,1)==CS_deg1_time(j,1) && Data_time_select(i,2)==CS_deg1_time(j,2)
%             CS_deg1s(i,:)=CS_deg1(j,:);
%             start_j=j+1;
%             break;
%         end
%     end
% end
% start_j=1;
% for i=1:Data_number
%     for j=start_j:length(C20_time(:,1))
%         if Data_time_select(i,1)==C20_time(j,1) && Data_time_select(i,2)==C20_time(j,2)
%             C20_coefs(i,:)=C20_coef(j,:);
%             start_j=j+1;
%             break;
%         end
%     end
% end
% GC(2,1,:)=CS_deg1s(:,1);GC(2,2,:)=CS_deg1s(:,2);GS(2,2,:)=CS_deg1s(:,3);
% GC(3,1,:)=C20_coefs(:,1);

% %扣除所有月份的平均值
if Data_number>1
    GC_mean=mean(GC,3);
    GS_mean=mean(GS,3);
    for i=1:Data_number
        GC(:,:,i)=GC(:,:,i)-GC_mean;
        GS(:,:,i)=GS(:,:,i)-GS_mean;
    end
end

%Research region: Logitude:0~360;Latitude:-90~90
% minlon=0;maxlon=360;
% minlat=-90;maxlat=90;
minlon=0.5;maxlon=359.5;
minlat=-89.5;maxlat=89.5;
% minlon=0.25;maxlon=359.75;
% minlat=-89.75;maxlat=89.75;
% minlon=70;maxlon=105;
% minlat=20;maxlat=40;
%Resolution of longitude and latitude:(Unit:degree);
Res_lonlat=1;
%Make grid for research region
[ceta,fir,n_c,n_f,cetax,firx,nceta,nfir]=region_grid(minlat,maxlat,minlon,maxlon,Res_lonlat);
disp('Region Gridded is ready!');

%Legendre polynom
Pnm=Nlmx_v3(Lmax,ceta);
% Pnm_ref=Nlmx_v3(Lmax_ref,ceta);
disp('Legendre polynom is ready!');
%%
% % %去相关滤波
% De_P=3;De_M=6;
% [GC_de,GS_de]=Dec4grace(Lmax,GC,GS,Data_number,De_P,De_M);
% GC=GC_de;GS=GS_de;
% disp('Results with P3M6 done!')
% % % Fan滤波
% r1_fan=300;r2_fan=300;%两个滤波半径
% [GC_fan,GS_fan]=Fan4grace(r1_fan,r2_fan,Lmax,GC,GS,Data_number);
% GC=GC_fan;GS=GS_fan;
% disp('Results with Fan filter done!');
% % %高斯滤波
% r_guassian=300;
% [GC_g,GS_g]=Gau4grace(r_guassian,Lmax,GC,GS,Data_number);
% GC=GC_g;GS=GS_g;
% disp('Results with Gauss filter done!');
%%
%球谐系数CS转换成等效水高EWH，单位cm
[grace_ewh]=cs2ewh(Lmax,Data_number,Pnm,GC,GS,fir,n_c,n_f,nceta,nfir);
disp('Results with cs2ewh done!');
for j=1:Data_number
    filename = strcat(Output_data_address,'test_60_EWH_',num2str(Data_time_select(j,1)),num2str(Data_time_select(j,2),'%02d'));
    gra_outfile([firx,cetax,grace_ewh(:,j)],filename,'txt',Res_lonlat);
end
%转成180*360的格网保存起来
model_cs2grid_140=zeros(180,360,Data_number);
for kk=1:length(grace_ewh(:,1))
    model_cs2grid_140(cetax(kk)+90.5,firx(kk)+0.5,:)=grace_ewh(kk,:);
end
Data_time_158=Data_time_select;
save model_cs2grid_140 model_cs2grid_140 Data_time_158;

disp('Finish!');
