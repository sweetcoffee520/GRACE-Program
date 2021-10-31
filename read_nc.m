% 2002��4��-2021��4��TWSA����,ȥ��������ͼ�����
clc;
clear;
%��ȡnc���͵��ļ�
path = 'E:\dataset\altimeter\CMEMS\REP\dataset-duacs-rep-global-merged-allsat-phy-l4-monthly\2016';%���ݵñ���·��
file = dir(fullfile(path,'*.nc'));%��������
len = length(file);	
for i = 1:1
    filename=fullfile(path,file(i).name);    % ƴ��·�����ļ���������ʾ
    info=ncinfo(filename);                 % ��ȡnc�ļ���Ϣ��������Ϣ�����info��
%     lat=ncread(filename,'LATITUDE');            % ��ȡnc�ļ��б���
%     lon=ncread(filename,'LONGITUDE');
    time=ncread(filename,'depth_bnds');
%     lwe_thickness = ncread(filename,'lwe_thickness');
end

[a,~]=size(time);
year=zeros(a,1);
mon=zeros(a,1);
day=zeros(a,1);
Num=datenum('1950/01/01'); 
for i=1:a
num=time(i)+Num;                              % ��������������Ӧ������ֵ
str=datestr(double(num),'yyyy/mm/dd');  % ת���ַ���ʽ
Datecellstr = cellstr(str);               % ���
date = cell2mat(Datecellstr);
year(i)=str2num(date(1:4));
mon(i)=str2num(date(6:7));
day(i)=str2num(date(9:10));
end
year=year;
mon=mon;
day=day;

yy = [year,mon,day];