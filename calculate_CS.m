function [c,s]=calculate_CS(n,Res_lonlat,griddata,k)
% %�������������λϵ��
% n�ף�n������200
% input: filename ��������,�ֱ��Ǿ��ȣ��ȣ���γ�ȣ��ȣ�,��Чˮ��(mm��/���ܶȣ�kg/m2��
% output: C,S n�׵�λϵ��
%%%%%��������
a=6.3781363000E+06; %���򳤰���
pave=5517.0;      %����ƽ���ܶ�
c=zeros(n+1,n+1,k);s=zeros(n+1,n+1,k);
% Res_lonlat=1;%�������
%��ȡ�շ�����Wahr,1998,���Բ�ֵ�õ�
filepath='D:\code\python_code\GRACE_py\data\loadLove.txt';
if exist(filepath,'file')==0
    error('�������ļ�');
else
    kl=load(filepath);
end
kl(2)=0.021;
%
% %��ȡ�������ݣ���һ���Ǿ��ȣ��ڶ�����γ��
% fid=fopen(filename);
% if fid==-1%���ļ�����ʧ���򷵻ؿ�ֵ
%     error('��ȡ�ļ�ʧ�ܣ�');
% end;
% num = 0;griddata=[];%griddata �洢��������
% while 1
%     line = fgetl(fid);
%     if line==-1%û�����ݾ�����ѭ��
%         break;
%     end;
%     griddata(:,num+1)=str2num(line);
%     num = num + 1;
% end
% fclose('all');
for ii=1:k
    [num,~]=size(griddata(:,:,ii));
    %�ǶȻ�����
    griddata(:,1,ii) = griddata(:,1,ii)/180*pi;
    griddata(:,2,ii) = griddata(:,2,ii)/180*pi;
    griddata(:,3,ii) = griddata(:,3,ii)*10;%ע������ĵ�λ��cm����mm,�����cm��Ҫ����10
    
    % c=zeros(n+1,n+1);s=zeros(n+1,n+1);
    for i = 1:num  %��num������������
        p=legendre_non(pi/2-griddata(i,2,ii),n);% ���õ�ϵ��
        %     save p%%�������
        for l = 0:n
            for m = 0:l
                c(l+1,m+1,ii)=c(l+1,m+1,ii)+(1/Res_lonlat)^2*griddata(i,3,ii)*p(l+1,m+1)*cos(m*griddata(i,1,ii))*cos(griddata(i,2,ii))/(180/pi)^2;
                s(l+1,m+1,ii)=s(l+1,m+1,ii)+(1/Res_lonlat)^2*griddata(i,3,ii)*p(l+1,m+1)*sin(m*griddata(i,1,ii))*cos(griddata(i,2,ii))/(180/pi)^2;
            end
        end
    end
    
    for l = 0:n
        for m = 0:l
            c(l+1,m+1,ii)=c(l+1,m+1,ii)*0.75*(1+kl(l+1))/(pi*a*pave*(2*l+1));%λϵ��C
            s(l+1,m+1,ii)=s(l+1,m+1,ii)*0.75*(1+kl(l+1))/(pi*a*pave*(2*l+1));%λϵ��S
        end
    end
end

