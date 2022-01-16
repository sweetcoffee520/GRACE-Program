function [c,s]=calculate_CS(n,Res_lonlat,griddata,k)
% %计算叠加重力场位系数
% n阶，n不超过200
% input: filename 格网数据,分别是经度（度），纬度（度）,等效水高(mm）/面密度（kg/m2）
% output: C,S n阶的位系数
%%%%%常量参数
a=6.3781363000E+06; %椭球长半轴
pave=5517.0;      %地球平均密度
c=zeros(n+1,n+1,k);s=zeros(n+1,n+1,k);
% Res_lonlat=1;%格网间隔
%读取勒夫数，Wahr,1998,线性插值得到
filepath='D:\code\python_code\GRACE_py\data\loadLove.txt';
if exist(filepath,'file')==0
    error('不存在文件');
else
    kl=load(filepath);
end
kl(2)=0.021;
%
% %读取格网数据，第一行是经度，第二行是纬度
% fid=fopen(filename);
% if fid==-1%打开文件，若失败则返回空值
%     error('读取文件失败！');
% end;
% num = 0;griddata=[];%griddata 存储格网数据
% while 1
%     line = fgetl(fid);
%     if line==-1%没有数据就跳出循环
%         break;
%     end;
%     griddata(:,num+1)=str2num(line);
%     num = num + 1;
% end
% fclose('all');
for ii=1:k
    [num,~]=size(griddata(:,:,ii));
    %角度换弧度
    griddata(:,1,ii) = griddata(:,1,ii)/180*pi;
    griddata(:,2,ii) = griddata(:,2,ii)/180*pi;
    griddata(:,3,ii) = griddata(:,3,ii)*10;%注意这里的单位是cm还是mm,如果是cm需要乘以10
    
    % c=zeros(n+1,n+1);s=zeros(n+1,n+1);
    for i = 1:num  %共num个格网点数据
        p=legendre_non(pi/2-griddata(i,2,ii),n);% 勒让德系数
        %     save p%%保存变量
        for l = 0:n
            for m = 0:l
                c(l+1,m+1,ii)=c(l+1,m+1,ii)+(1/Res_lonlat)^2*griddata(i,3,ii)*p(l+1,m+1)*cos(m*griddata(i,1,ii))*cos(griddata(i,2,ii))/(180/pi)^2;
                s(l+1,m+1,ii)=s(l+1,m+1,ii)+(1/Res_lonlat)^2*griddata(i,3,ii)*p(l+1,m+1)*sin(m*griddata(i,1,ii))*cos(griddata(i,2,ii))/(180/pi)^2;
            end
        end
    end
    
    for l = 0:n
        for m = 0:l
            c(l+1,m+1,ii)=c(l+1,m+1,ii)*0.75*(1+kl(l+1))/(pi*a*pave*(2*l+1));%位系数C
            s(l+1,m+1,ii)=s(l+1,m+1,ii)*0.75*(1+kl(l+1))/(pi*a*pave*(2*l+1));%位系数S
        end
    end
end

