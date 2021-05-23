function [GLDAS_data,GLDAS_number]=read_GLDAS_nc4(Address)
%读取GLDAS水文模型，数据格式是.nc4
% ncdisp('E:\GLDAS\water\datefile\GLDAS_NOAH10_M.A200001.021.nc4')%查看变量的详细信息
% data_name=gldas.Variables(ii).Name;
% ii=20:23;%土壤水和温度都是分四层
% Address='E:\GLDAS\water\datafile\GSM2011\';%分辨率为1°的数据
% Address='D:\数据\GLDAS水文模型\GLDAS_NOAH025_M.2.1\all_data\';%分辨率为0.25°的数据
scale_factor=1;%比例因子，分辨率为1°对应于1，0.25°对应于4
Output_data_address='E:\Results\Tibet\';
GFA=dir(fullfile(Address,'*.nc4'));
GLDAS_number=length(GFA);
filepath=Address;
 GLDAS_time=GrabTheDate_GLDAS_nc4(filepath);
for i=1:GLDAS_number
    file=[filepath,GFA(i).name];
%     yr_mon=str2double(file(length(file)-13:length(file)-8));
%     day=yr_mon2day(yr_mon);
    %%%%%%  读取数据
    xx=ncinfo(file);
    lat=ncread(file,xx.Variables(1).Name);%纬度范围是-59.5~89.5
    lon=ncread(file,xx.Variables(2).Name);%经度,范围是-179.5到179.5，注意排列顺序
    %把经度的范围转换成0.5~359.5，注意排列顺序，先排列180.5到359.5度，然后是0.5到179.5度，对应的格网值后面也进行了调换
    for ii=1:length(lon)
        if lon(ii)<0
            lon(ii)=lon(ii)+360;
        end
    end
    tmp_snow_rate=ncread(file,xx.Variables(10).Name);%降雪的速率
    soimoi1=ncread(file,xx.Variables(20).Name);%第一层 0-10cm 平均土壤水分,单位是kg/m^2，可相当于mm
    soimoi1(isnan(soimoi1)==1) = 0;%NaN的地方赋值为0
    soimoi2=ncread(file,xx.Variables(21).Name);%第二层 10-40cm 平均土壤水分,单位是kg/m^2，可相当于mm
    soimoi2(isnan(soimoi2)==1) = 0;
    soimoi3=ncread(file,xx.Variables(22).Name);%第三层 40-100cm 平均土壤水分,单位是kg/m^2，可相当于mm
    soimoi3(isnan(soimoi3)==1) = 0;
    soimoi4=ncread(file,xx.Variables(23).Name);%第四层 100-200cm 平均土壤水分,单位是kg/m^2，可相当于mm
    soimoi4(isnan(soimoi4)==1) = 0;
    tmp_soimoi=(soimoi1+soimoi2+soimoi3+soimoi4)/10;%四层累加起来，并转化成等效水高，除以10，单位化成cm
    
    tmp_surface_runoff=ncread(file,xx.Variables(13).Name)/10;%地表径流,单位是kg/m/m，转化成等效水高，单位cm，暂时不用
    
    tmp_surface_snow=ncread(file,xx.Variables(18).Name)/10;%雪水当量,单位是kg/m/m，转化成等效水高，单位cm
    tmp_surface_snow(isnan(tmp_surface_snow)==1) = 0;
    tmp_canopy_water=ncread(file,xx.Variables(33).Name)/10;%植物冠层含水量,单位是kg/m/m，转化成等效水高，单位cm
    tmp_canopy_water(isnan(tmp_canopy_water)==1) = 0;

  tmp_precipitation_rate=ncread(file,xx.Variables(35).Name)/10;%总降雨量,单位是kg/m-2/s，转化成等效水高，单位cm
    tmp_sum_water=tmp_soimoi+tmp_surface_snow+tmp_canopy_water;%GLDAS水文模型总的水量，注意单位变成mm需要整体乘以10
    
%     为与GRACE的输出结果保持一致，对二维结果转换成一维
    for ii=1:150*scale_factor
        for jj=1:360*scale_factor
            sum_water((ii-1)*360*scale_factor+jj,i)=tmp_sum_water(jj,ii);
        end
    end
    %把经度的的排列顺序调整为0.5~359.5
    tmp_sum_water1=[tmp_sum_water(180*scale_factor+1:360*scale_factor,:);tmp_sum_water(1:180*scale_factor,:)];
    GLDAS_data0(:,:,i)=tmp_sum_water1';
end
% %扣除这期间的月平均值，得到陆地水储量逐月变化量
 Delta_GLDAS=sum_water-repmat(mean(sum_water,2),1,GLDAS_number);
% Delta_GLDAS=sum_water;%不扣除平均值

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
Res_lonlat=1;%分辨率为1度，这个参数好像没啥用
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
save GLDAS_data GLDAS_data0 GLDAS_time;%保存到矩阵的数据是没有减去平均值的
disp('Finish!');
