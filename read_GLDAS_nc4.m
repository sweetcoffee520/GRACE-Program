function [GLDAS_data,GLDAS_number]=read_GLDAS_nc4(Address)
%��ȡGLDASˮ��ģ�ͣ����ݸ�ʽ��.nc4
% ncdisp('E:\GLDAS\water\datefile\GLDAS_NOAH10_M.A200001.021.nc4')%�鿴��������ϸ��Ϣ
% data_name=gldas.Variables(ii).Name;
% ii=20:23;%����ˮ���¶ȶ��Ƿ��Ĳ�
% Address='E:\GLDAS\water\datafile\GSM2011\';%�ֱ���Ϊ1�������
% Address='D:\����\GLDASˮ��ģ��\GLDAS_NOAH025_M.2.1\all_data\';%�ֱ���Ϊ0.25�������
scale_factor=1;%�������ӣ��ֱ���Ϊ1���Ӧ��1��0.25���Ӧ��4
Output_data_address='E:\Results\Tibet\';
GFA=dir(fullfile(Address,'*.nc4'));
GLDAS_number=length(GFA);
filepath=Address;
 GLDAS_time=GrabTheDate_GLDAS_nc4(filepath);
for i=1:GLDAS_number
    file=[filepath,GFA(i).name];
%     yr_mon=str2double(file(length(file)-13:length(file)-8));
%     day=yr_mon2day(yr_mon);
    %%%%%%  ��ȡ����
    xx=ncinfo(file);
    lat=ncread(file,xx.Variables(1).Name);%γ�ȷ�Χ��-59.5~89.5
    lon=ncread(file,xx.Variables(2).Name);%����,��Χ��-179.5��179.5��ע������˳��
    %�Ѿ��ȵķ�Χת����0.5~359.5��ע������˳��������180.5��359.5�ȣ�Ȼ����0.5��179.5�ȣ���Ӧ�ĸ���ֵ����Ҳ�����˵���
    for ii=1:length(lon)
        if lon(ii)<0
            lon(ii)=lon(ii)+360;
        end
    end
    tmp_snow_rate=ncread(file,xx.Variables(10).Name);%��ѩ������
    soimoi1=ncread(file,xx.Variables(20).Name);%��һ�� 0-10cm ƽ������ˮ��,��λ��kg/m^2�����൱��mm
    soimoi1(isnan(soimoi1)==1) = 0;%NaN�ĵط���ֵΪ0
    soimoi2=ncread(file,xx.Variables(21).Name);%�ڶ��� 10-40cm ƽ������ˮ��,��λ��kg/m^2�����൱��mm
    soimoi2(isnan(soimoi2)==1) = 0;
    soimoi3=ncread(file,xx.Variables(22).Name);%������ 40-100cm ƽ������ˮ��,��λ��kg/m^2�����൱��mm
    soimoi3(isnan(soimoi3)==1) = 0;
    soimoi4=ncread(file,xx.Variables(23).Name);%���Ĳ� 100-200cm ƽ������ˮ��,��λ��kg/m^2�����൱��mm
    soimoi4(isnan(soimoi4)==1) = 0;
    tmp_soimoi=(soimoi1+soimoi2+soimoi3+soimoi4)/10;%�Ĳ��ۼ���������ת���ɵ�Чˮ�ߣ�����10����λ����cm
    
    tmp_surface_runoff=ncread(file,xx.Variables(13).Name)/10;%�ر���,��λ��kg/m/m��ת���ɵ�Чˮ�ߣ���λcm����ʱ����
    
    tmp_surface_snow=ncread(file,xx.Variables(18).Name)/10;%ѩˮ����,��λ��kg/m/m��ת���ɵ�Чˮ�ߣ���λcm
    tmp_surface_snow(isnan(tmp_surface_snow)==1) = 0;
    tmp_canopy_water=ncread(file,xx.Variables(33).Name)/10;%ֲ��ڲ㺬ˮ��,��λ��kg/m/m��ת���ɵ�Чˮ�ߣ���λcm
    tmp_canopy_water(isnan(tmp_canopy_water)==1) = 0;

  tmp_precipitation_rate=ncread(file,xx.Variables(35).Name)/10;%�ܽ�����,��λ��kg/m-2/s��ת���ɵ�Чˮ�ߣ���λcm
    tmp_sum_water=tmp_soimoi+tmp_surface_snow+tmp_canopy_water;%GLDASˮ��ģ���ܵ�ˮ����ע�ⵥλ���mm��Ҫ�������10
    
%     Ϊ��GRACE������������һ�£��Զ�ά���ת����һά
    for ii=1:150*scale_factor
        for jj=1:360*scale_factor
            sum_water((ii-1)*360*scale_factor+jj,i)=tmp_sum_water(jj,ii);
        end
    end
    %�Ѿ��ȵĵ�����˳�����Ϊ0.5~359.5
    tmp_sum_water1=[tmp_sum_water(180*scale_factor+1:360*scale_factor,:);tmp_sum_water(1:180*scale_factor,:)];
    GLDAS_data0(:,:,i)=tmp_sum_water1';
end
% %�۳����ڼ����ƽ��ֵ���õ�½��ˮ�������±仯��
 Delta_GLDAS=sum_water-repmat(mean(sum_water,2),1,GLDAS_number);
% Delta_GLDAS=sum_water;%���۳�ƽ��ֵ

tmp_latx=zeros(1,150*scale_factor*360*scale_factor);
tmp_lonx=zeros(1,150*scale_factor*360*scale_factor);
for ii=1:150*scale_factor
    for jj=1:360*scale_factor
        tmp_latx((ii-1)*360*scale_factor+jj)=lat(ii);
        tmp_lonx((ii-1)*360*scale_factor+jj)=lon(jj);
    end
end

lonx=tmp_lonx';
latx=tmp_latx';
Res_lonlat=1;%�ֱ���Ϊ1�ȣ������������ûɶ��
% save data fo all monthly GLDAS data
for i=1:GLDAS_number
    data_tmp=[lonx,latx,Delta_GLDAS(:,i)];
    A=isnan(data_tmp(:,3));
    ind=find(isnan(data_tmp(:,3)));
    data_tmp(ind,:)=[];
    GLDAS_data_158(:,:,i)=data_tmp;
    filename = strcat(Output_data_address,'GLDAS_NOAH10_M.A',num2str(GLDAS_time(i,1)),num2str(GLDAS_time(i,2),'%02d'));
%     filename = strcat(Output_data_address,'GLDAS_NOAH025_M.A',num2str(GLDAS_time(i,1)),num2str(GLDAS_time(i,2),'%02d'));
    gra_outfile([lonx,latx,Delta_GLDAS(:,i)],filename,'txt',Res_lonlat);
    GLDAS_data(:,:,i)=[data_tmp(:,1),data_tmp(:,2),data_tmp(:,3)];
end
save GLDAS_data GLDAS_data0 GLDAS_time;%���浽�����������û�м�ȥƽ��ֵ��
disp('Finish!');
