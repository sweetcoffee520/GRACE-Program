function info = read_nc()
%读取nc类型的文件
path = 'I:\86442\data\field\2002';%数据得保存路径
file = dir(fullfile(path,'*.nc'));%数据名称
len = length(file);
for i = 1:1
    filename=fullfile(path,file(i).name);    % 拼接路径和文件名，并显示
    info=ncinfo(filename);                 % 读取nc文件信息，变量信息存放在info中
    %     lat=ncread(filename,'LATITUDE');            % 提取nc文件中变量
    %     lon=ncread(filename,'LONGITUDE');
    time=ncread(filename,'time');
    depth = ncread(filename,'depth');
    %     lwe_thickness = ncread(filename,'lwe_thickness');
end

% [a,~]=size(time);
% year=zeros(a,1);
% mon=zeros(a,1);
% day=zeros(a,1);
% Num=datenum('1990/01/01');
% for i=1:a
%     num=time(i)+Num;                              % 计算后面的天数对应的日期值
%     str=datestr(double(num),'yyyy/mm/dd');  % 转成字符形式
%     Datecellstr = cellstr(str);               % 输出
%     date = cell2mat(Datecellstr);
%     year(i)=str2num(date(1:4));
%     mon(i)=str2num(date(6:7));
%     day(i)=str2num(date(9:10));
% end
% year=year;
% mon=mon;
% day=day;
% 
% yy = [year,mon,day];