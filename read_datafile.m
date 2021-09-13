function [DeltaGC,DeltaGS,DeltaRC,DeltaRS,DeltaGC1,DeltaGS1,DeltaRC1,DeltaRS1,...
    DeltaGC2,DeltaGS2,DeltaRC2,DeltaRS2,DeltaGC3,DeltaGS3,DeltaRC3,DeltaRS3,...
    DeltaGC4,DeltaGS4,DeltaRC4,DeltaRS4,DeltaGC5,DeltaGS5,DeltaRC5,DeltaRS5,...
    DeltaGC6,DeltaGS6,DeltaGC7,DeltaGS7]=read_datafile()
%%
%%%%%read GRACE-FO's data
Lmax=60;
Address='C:\GLDAS\water\数据\CSR-60阶\';
GFA = dir(fullfile(Address,'GSM*'));
k=length(GFA);
GC=zeros((Lmax+3)*Lmax/2+1,k);
GS=zeros((Lmax+3)*Lmax/2+1,k);
RC=zeros((Lmax+3)*Lmax/2+1,k);
RS=zeros((Lmax+3)*Lmax/2+1,k);
DeltaGC=zeros(Lmax+1,Lmax+1,k);
DeltaGS=zeros(Lmax+1,Lmax+1,k);
DeltaRC=zeros(Lmax+1,Lmax+1,k);
DeltaRS=zeros(Lmax+1,Lmax+1,k);
for i=1:k
    filepath=Address;
    file=[filepath,GFA(i).name];
    fid=fopen(file);
    if fid == -1
        'Error opening the file.';
    end 
    
    k_line=0;
    while 1

        nextline=fgetl(fid);      %读第一行
        if ~ischar(nextline)
            break,
        end %读到最后跳出
        if ~strcmp(nextline,'')
            if strcmp(nextline(1:6),'GRCOF2')
%         if length(nextline)==40      %笨方法，已知数据长度可以用
            k_line=k_line+1;
%            while k_line<=2||k_line==92
%                k_line=k_line+1;
%            end
           GC(k_line,i)=str2double(nextline(18:35));
           GS(k_line,i)=str2double(nextline(37:54));
           RC(k_line,i)=str2double(nextline(56:65));
           RS(k_line,i)=str2double(nextline(67:76));
            end
        end
    end
    sta=fclose(fid);
end
%%%%三次样条插值求出缺失项

C20_gra=read_SLR('C:\GLDAS\water\datafile\TN-11_C20_SLR_RL06.txt');
GC(3,1:9)=C20_gra(164:end);
for i=1:k
    for l=0:Lmax
        for m=0:l
            DeltaRC(l+1,m+1,i)=RC(l+1+(Lmax*2+1-m)*m/2,i);
            DeltaRS(l+1,m+1,i)=RS(l+1+(Lmax*2+1-m)*m/2,i);
        end
    end
end

[~,k]=size(GC);
Dgc=zeros((Lmax+3)*Lmax/2+1,k);
Dgs=zeros((Lmax+3)*Lmax/2+1,k);

GC0=0;
for i=1:k
    GC0=GC0+GC(:,i);
end
GC00=GC0/k;

GS0=0;
for i=1:k
    GS0=GS0+GS(:,i);
end
GS00=GS0/k;
for i=1:k
    Dgc(:,i)=GC(:,i)-GC00;
    Dgs(:,i)=GS(:,i)-GS00;
end
for i=1:k
    for l=0:Lmax
        for m=0:l
            DeltaGC(l+1,m+1,i)=Dgc(l+1+(Lmax*2+1-m)*m/2,i);
            DeltaGS(l+1,m+1,i)=Dgs(l+1+(Lmax*2+1-m)*m/2,i);
        end
    end
end

clearvars -except DeltaGC DeltaGS DeltaRC DeltaRS;
%%
%%%%read SWMA's data
% Lmax=70;
Lmax=40;
Address=['C:\GLDAS\water\数据\AIUB-70\A\';'C:\GLDAS\water\数据\AIUB-70\B\';...
    'C:\GLDAS\water\数据\AIUB-70\C\'];
GFA = dir(fullfile(Address(1,:),'SWM*'));
k=length(GFA);
ll=zeros((Lmax+3)*Lmax/2+1,k);
mm=zeros((Lmax+3)*Lmax/2+1,k);
GC=zeros((Lmax+3)*Lmax/2+1,k);
GS=zeros((Lmax+3)*Lmax/2+1,k);
RC=zeros((Lmax+3)*Lmax/2+1,k);
RS=zeros((Lmax+3)*Lmax/2+1,k);
DeltaGC1=zeros(Lmax+1,Lmax+1,k);
DeltaGS1=zeros(Lmax+1,Lmax+1,k);
DeltaRC1=zeros(Lmax+1,Lmax+1,k);
DeltaRS1=zeros(Lmax+1,Lmax+1,k);
for z=1:3
    GFA = dir(fullfile(Address(z,:),'SWM*'));
    for i=1:k
        filepath=Address(z,:);
        file=[filepath,GFA(i).name];
        fid=fopen(file);
        if fid == -1
            'Error opening the file.';
        end
        
        k_line=0;
        while 1
            
            nextline=fgetl(fid);      %读第一行
            if ~ischar(nextline)
                break,
            end %读到最后跳出
            if ~strcmp(nextline,'')
                if strcmp(nextline(1:4),'GEOC')
                    k_line1=str2double(nextline(9:11));
                    if  k_line1>40
                        break;
                    end
                    k_line=k_line+1;
                    ll(k_line,i,z)=str2double(nextline(9:11));
                    mm(k_line,i,z)=str2double(nextline(14:16));
                    GC(k_line,i,z)=str2double(nextline(18:35));
                    GS(k_line,i,z)=str2double(nextline(37:54));
                    RC(k_line,i,z)=str2double(nextline(56:73));
                    RS(k_line,i,z)=str2double(nextline(75:92));
                end
            end
        end
        sta=fclose(fid);
    end
end
numerator=1./(RC.*RC);
numerator(numerator==inf)=0;
denominator=sum(numerator,3);
denominator(denominator==inf)=0;
pc=numerator./denominator;
pc(isnan(pc)==1)=0;
gc=sum(GC.*pc,3);

numerator1=1./(RS.*RS);
numerator1(numerator1==inf)=0;
denominator1=sum(numerator1,3);
denominator1(denominator1==inf)=0;
ps=numerator1./denominator1;
ps(isnan(ps)==1)=0;
gs=sum(GS.*ps,3);

mc=sqrt(sum((RC.*RC).*(pc.*pc),3)./(sum(pc,3).*sum(pc,3)));
mc(isnan(mc)==1)=0;
ms=sqrt(sum((RS.*RS).*(ps.*ps),3)./(sum(ps,3).*sum(ps,3)));
ms(isnan(ms)==1)=0;

% Dgc=zeros((Lmax+3)*Lmax/2+1,k);
% Dgs=zeros((Lmax+3)*Lmax/2+1,k);
Dgc=zeros((Lmax+3)*Lmax/2+1,k);
Dgs=zeros((Lmax+3)*Lmax/2+1,k);

GC0=0;
for i=1:k
    GC0=GC0+gc(:,i);
end
GC00=GC0/k;

GS0=0;
for i=1:k
    GS0=GS0+gs(:,i);
end
GS00=GS0/k;
for i=1:k
    Dgc(:,i)=gc(:,i)-GC00;
    Dgs(:,i)=gs(:,i)-GS00;
end

for i=1:k
    for j=1:length(gc)
        DeltaGC1(ll(j,i,1)+1,mm(j,i,1)+1,i)=Dgc(j,i);
        DeltaGS1(ll(j,i,1)+1,mm(j,i,1)+1,i)=Dgs(j,i);
        DeltaRC1(ll(j,i,1)+1,mm(j,i,1)+1,i)=mc(j,i);
        DeltaRS1(ll(j,i,1)+1,mm(j,i,1)+1,i)=ms(j,i);
    end
end

clearvars -except DeltaGC DeltaGS DeltaRC DeltaRS DeltaGC1 DeltaGS1 DeltaRC1 DeltaRS1;
%%
%%%%%%read ASU-Swarm data
% Lmax=40;
Lmax=15;
Address='C:\GLDAS\water\数据\ASU-Swarm-40\';
GFA = dir(fullfile(Address,'GSW*'));
k=length(GFA);
GC=zeros((Lmax+3)*Lmax/2+1,k);
GS=zeros((Lmax+3)*Lmax/2+1,k);
RC=zeros((Lmax+3)*Lmax/2+1,k);
RS=zeros((Lmax+3)*Lmax/2+1,k);
DeltaGC2=zeros(Lmax+1,Lmax+1,k);
DeltaGS2=zeros(Lmax+1,Lmax+1,k);
DeltaRC2=zeros(Lmax+1,Lmax+1,k);
DeltaRS2=zeros(Lmax+1,Lmax+1,k);
for i=1:k
filepath=Address;
file=[filepath,GFA(i).name];
    fid=fopen(file);
if fid == -1
    'Error opening the file.';
end
k_line=0;
while 1
    
    nextline=fgetl(fid);      %读第一行
    if ~ischar(nextline)
        break,
    end %读到最后跳出
    if ~strcmp(nextline,'')
        if strcmp(nextline(1:3),'gfc')
            k_line1=str2double(nextline(6:7));
            if  k_line1>15
                break;
            end
            k_line=k_line+1;
            ll(k_line,i)=str2double(nextline(6:7));
            mm(k_line,i)=str2double(nextline(10:11));
            GC(k_line,i)=str2double(nextline(14:35));
            GS(k_line,i)=str2double(nextline(38:59));
            RC(k_line,i)=str2double(nextline(62:74));
            RS(k_line,i)=str2double(nextline(77:89));
        end
    end
end
sta=fclose(fid);
end
% C20_gra=read_SLR('C:\GLDAS\water\数据\TN-11_C20_SLR_RL06.txt');
% temp=importdata('C:\GLDAS\water\数据\TN-13_GEOC_CSR_RL06_2019-04-22.txt',' ',113);
% GC_Replace=temp.data;
% temp=[1:2,4:7,9:12,14:18,20:22,25:28,30:33,36:38,40:43];
% GC(3,temp)=C20_gra(131:163);
% k_line=0;
% k_line1=0;
% for i=1:size(GC_Replace,1)
%     if GC_Replace(i,2)==0
%         k_line=k_line+1;
%         GC10(k_line)=GC_Replace(i,3);
%         GS10(k_line)=GC_Replace(i,4);
%     elseif GC_Replace(i,2)==1
%         k_line1=k_line1+1;
%         GC11(k_line1)=GC_Replace(i,3);
%         GS11(k_line1)=GC_Replace(i,4);
%     end
% end
% GC(2,temp)=GC10(131:end);
% GS(2,temp)=GS10(131:end);
% GC(62,temp)=GC11(131:end);
% GS(62,temp)=GS11(131:end);
% Dgc=zeros((Lmax+3)*Lmax/2+1,k);
% Dgs=zeros((Lmax+3)*Lmax/2+1,k);
Dgc=zeros((Lmax+3)*Lmax/2+1,k);
Dgs=zeros((Lmax+3)*Lmax/2+1,k);

GC0=0;
for i=1:k
    GC0=GC0+GC(:,i);
end
GC00=GC0/k;

GS0=0;
for i=1:k
    GS0=GS0+GS(:,i);
end
GS00=GS0/k;
for i=1:k
    Dgc(:,i)=GC(:,i)-GC00;
    Dgs(:,i)=GS(:,i)-GS00;
end
for j=1:k
    for i=1:length(GC)
        DeltaGC2(ll(i)+1,mm(i)+1,j)=Dgc(i,j);
        DeltaGS2(ll(i)+1,mm(i)+1,j)=Dgs(i,j);
        DeltaRC2(ll(i)+1,mm(i)+1,j)=RC(i,j);
        DeltaRS2(ll(i)+1,mm(i)+1,j)=RS(i,j);
    end
end
clearvars -except DeltaGC DeltaGS DeltaRC DeltaRS DeltaGC1 DeltaGS1 DeltaRC1 DeltaRS1...
    DeltaGC2 DeltaGS2 DeltaRC2 DeltaRS2;
%%
%%%%read Grace data
Lmax=60;
Address='C:\GLDAS\water\数据\60阶\';
% GIA=dir(fullfile(FAddress,'xx*));  %%xx表示路径,FAddress表示父文件夹
% yr=length(GIA);
% for t=1:yr
% Ffilepath=FAddress;
% Adress=[Ffilepath,GIA(t).name);
GFA = dir(fullfile(Address,'GSM*'));
k=length(GFA);
GC=zeros((Lmax+3)*Lmax/2+1,k);
GS=zeros((Lmax+3)*Lmax/2+1,k);
RC=zeros((Lmax+3)*Lmax/2+1,k);
RS=zeros((Lmax+3)*Lmax/2+1,k);
DeltaGC3=zeros(Lmax+1,Lmax+1,k);
DeltaGS3=zeros(Lmax+1,Lmax+1,k);
DeltaRC3=zeros(Lmax+1,Lmax+1,k);
DeltaRS3=zeros(Lmax+1,Lmax+1,k);
for i=1:k
    filepath=Address;
    file=[filepath,GFA(i).name];
    fid=fopen(file);
    if fid == -1
        'Error opening the file.';
    end
    
    k_line=0;
    while 1
        
        nextline=fgetl(fid);      %读第一行
        if ~ischar(nextline)
            break,
        end %读到最后跳出
        if ~strcmp(nextline,'')
            if strcmp(nextline(1:6),'GRCOF2')
                %         if length(nextline)==40      %笨方法，已知数据长度可以用
                k_line=k_line+1;
                %            while k_line<=2||k_line==92
                %                k_line=k_line+1;
                %            end
                GC(k_line,i)=str2double(nextline(18:35));
                GS(k_line,i)=str2double(nextline(37:54));
                RC(k_line,i)=str2double(nextline(56:65));
                RS(k_line,i)=str2double(nextline(67:76));
            end
        end
    end
    sta=fclose(fid);
end
%%%%三次样条插值求出缺失项
%x=[1:5,7:96,98:101,103:112,114:117,119:122,124:127,129:132,134:137,139:142,144:148,150:152,155:158,160:163];
%x=[1:2,4:7,9:12,14:18,20:22,25:28,30:33,36:38,40:43];
% [GC10,GS10,GC11,GS11]=read_deg1('C:\GLDAS\water\数据\deg1_coef.txt');
C20_gra=read_SLR('C:\GLDAS\water\数据\TN-11_C20_SLR_RL06.txt');
GC(3,:)=C20_gra(1:163);
temp=importdata('C:\GLDAS\water\数据\TN-13_GEOC_CSR_RL06_2019-04-22.txt',' ',113);
GC_Replace=temp.data;
k_line=0;
k_line1=0;
for i=1:size(GC_Replace,1)
    if GC_Replace(i,2)==0
        k_line=k_line+1;
        GC10(k_line)=GC_Replace(i,3);
        GS10(k_line)=GC_Replace(i,4);
    elseif GC_Replace(i,2)==1
        k_line1=k_line1+1;
        GC11(k_line1)=GC_Replace(i,3);
        GS11(k_line1)=GC_Replace(i,4);
    end
end
GC(2,:)=GC10(1:end);
GS(2,:)=GS10(1:end);
GC(62,:)=GC11(1:end);
GS(62,:)=GS11(1:end);
%GS(62,[1:11,13:end])=GS11(6:153);

for i=1:k
    for l=0:Lmax
        for m=0:l
            DeltaRC3(l+1,m+1,i)=RC(l+1+(Lmax*2+1-m)*m/2,i);
            DeltaRS3(l+1,m+1,i)=RS(l+1+(Lmax*2+1-m)*m/2,i);
        end
    end
end

[~,k]=size(GC);
% Dgc=zeros((Lmax+3)*Lmax/2+1,k);
% Dgs=zeros((Lmax+3)*Lmax/2+1,k);
Dgc=zeros((Lmax+3)*Lmax/2+1,k);
Dgs=zeros((Lmax+3)*Lmax/2+1,k);

GC0=0;
for i=1:k
    GC0=GC0+GC(:,i);
end
GC00=GC0/k;

GS0=0;
for i=1:k
    GS0=GS0+GS(:,i);
end
GS00=GS0/k;
for i=1:k
    Dgc(:,i)=GC(:,i)-GC00;
    Dgs(:,i)=GS(:,i)-GS00;
end
for i=1:k
    for l=0:Lmax
        for m=0:l
            DeltaGC3(l+1,m+1,i)=Dgc(l+1+(Lmax*2+1-m)*m/2,i);
            DeltaGS3(l+1,m+1,i)=Dgs(l+1+(Lmax*2+1-m)*m/2,i);
        end
    end
end
clearvars -except DeltaGC DeltaGS DeltaRC DeltaRS DeltaGC1 DeltaGS1 DeltaRC1 DeltaRS1...
    DeltaGC2 DeltaGS2 DeltaRC2 DeltaRS2 DeltaGC3 DeltaGS3 DeltaRC3 DeltaRS3;

%%
%%%%%%read ITSG-60 data
% Lmax=60;
Lmax=40;
Address='C:\GLDAS\water\数据\ITSG-60\v2\';
GFA = dir(fullfile(Address,'coeff*'));
k=length(GFA);
GC=zeros((Lmax+3)*Lmax/2+1,k);
GS=zeros((Lmax+3)*Lmax/2+1,k);
RC=zeros((Lmax+3)*Lmax/2+1,k);
RS=zeros((Lmax+3)*Lmax/2+1,k);
DeltaGC4=zeros(Lmax+1,Lmax+1,k);
DeltaGS4=zeros(Lmax+1,Lmax+1,k);
DeltaRC4=zeros(Lmax+1,Lmax+1,k);
DeltaRS4=zeros(Lmax+1,Lmax+1,k);
for i=1:k
filepath=Address;
file=[filepath,GFA(i).name];
    fid=fopen(file);
if fid == -1
    'Error opening the file.';
end
k_line=0;
while 1
    
    nextline=fgetl(fid);      %读第一行
    if ~ischar(nextline)
        break,
    end %读到最后跳出
    if ~strcmp(nextline,'')
        if strcmp(nextline(1:3),'gfc')
            k_line1=str2double(nextline(8:9));
            if  k_line1>40
                break;
            end
            k_line=k_line+1;
                ll(k_line,i)=str2double(nextline(8:9));
                mm(k_line,i)=str2double(nextline(13:14));
                GC(k_line,i)=str2double(nextline(16:34));
                GS(k_line,i)=str2double(nextline(36:54));
                RC(k_line,i)=str2double(nextline(57:74));
                RS(k_line,i)=str2double(nextline(77:94));
            end
        end
    end
sta=fclose(fid);
end
% Dgc=zeros((Lmax+3)*Lmax/2+1,k);
% Dgs=zeros((Lmax+3)*Lmax/2+1,k);
Dgc=zeros((Lmax+3)*Lmax/2+1,k);
Dgs=zeros((Lmax+3)*Lmax/2+1,k);

GC0=0;
for i=1:k
    GC0=GC0+GC(:,i);
end
GC00=GC0/k;

GS0=0;
for i=1:k
    GS0=GS0+GS(:,i);
end
GS00=GS0/k;
for i=1:k
    Dgc(:,i)=GC(:,i)-GC00;
    Dgs(:,i)=GS(:,i)-GS00;
end
for j=1:k
    for i=1:size(GC)
        DeltaGC4(ll(i,1)+1,mm(i,1)+1,j)=Dgc(i,j);
        DeltaGS4(ll(i,1)+1,mm(i,1)+1,j)=Dgs(i,j);
        DeltaRC4(ll(i,1)+1,mm(i,1)+1,j)=RC(i,j);
        DeltaRS4(ll(i,1)+1,mm(i,1)+1,j)=RS(i,j);
    end
end

clearvars -except DeltaGC DeltaGS DeltaRC DeltaRS DeltaGC1 DeltaGS1 DeltaRC1 DeltaRS1...
    DeltaGC2 DeltaGS2 DeltaRC2 DeltaRS2 DeltaGC3 DeltaGS3 DeltaRC3 DeltaRS3...
    DeltaGC4 DeltaGS4 DeltaRC4 DeltaRS4;

%%
% Lmax=60;
Lmax=40;
Address='C:\GLDAS\water\数据\ULux_CHAMP\';
GFA = dir(fullfile(Address,'ULux*'));
k=length(GFA);
GC=zeros((Lmax+3)*Lmax/2+1,k);
GS=zeros((Lmax+3)*Lmax/2+1,k);
RC=zeros((Lmax+3)*Lmax/2+1,k);
RS=zeros((Lmax+3)*Lmax/2+1,k);
DeltaGC5=zeros(Lmax+1,Lmax+1,k);
DeltaGS5=zeros(Lmax+1,Lmax+1,k);
DeltaRC5=zeros(Lmax+1,Lmax+1,k);
DeltaRS5=zeros(Lmax+1,Lmax+1,k);
for i=1:k
filepath=Address;
file=[filepath,GFA(i).name];
    fid=fopen(file);
if fid == -1
    'Error opening the file.';
end
k_line=0;
while 1
    
    nextline=fgetl(fid);      %读第一行
    if ~ischar(nextline)
        break,
    end %读到最后跳出
    if ~strcmp(nextline,'')
        if strcmp(nextline(1:3),'gfc')
            k_line1=str2double(nextline(7:8));
            if  k_line1>40
                break;
            end
            k_line=k_line+1;
                ll(k_line,i)=str2double(nextline(7:8));
                mm(k_line,i)=str2double(nextline(12:13));
                GC(k_line,i)=str2double(nextline(15:33));
                GS(k_line,i)=str2double(nextline(35:53));
                RC(k_line,i)=str2double(nextline(55:73));
                RS(k_line,i)=str2double(nextline(75:93));
            end
        end
    end
sta=fclose(fid);
end
% Dgc=zeros((Lmax+3)*Lmax/2+1,k);
% Dgs=zeros((Lmax+3)*Lmax/2+1,k);
Dgc=zeros((Lmax+3)*Lmax/2+1,k);
Dgs=zeros((Lmax+3)*Lmax/2+1,k);

GC0=0;
for i=1:k
    GC0=GC0+GC(:,i);
end
GC00=GC0/k;

GS0=0;
for i=1:k
    GS0=GS0+GS(:,i);
end
GS00=GS0/k;
for i=1:k
    Dgc(:,i)=GC(:,i)-GC00;
    Dgs(:,i)=GS(:,i)-GS00;
end
for j=1:k
    for i=1:length(GC)
        DeltaGC5(ll(i,1)+1,mm(i,1)+1,j)=Dgc(i,j);
        DeltaGS5(ll(i,1)+1,mm(i,1)+1,j)=Dgs(i,j);
        DeltaRC5(ll(i,1)+1,mm(i,1)+1,j)=RC(i,j);
        DeltaRS5(ll(i,1)+1,mm(i,1)+1,j)=RS(i,j);
    end
end
clearvars -except DeltaGC DeltaGS DeltaRC DeltaRS DeltaGC1 DeltaGS1 DeltaRC1 DeltaRS1...
    DeltaGC2 DeltaGS2 DeltaRC2 DeltaRS2 DeltaGC3 DeltaGS3 DeltaRC3 DeltaRS3...
    DeltaGC4 DeltaGS4 DeltaRC4 DeltaRS4 DeltaGC5 DeltaGS5 DeltaRC5 DeltaRS5;

%%
% Lmax=40;
Lmax=15;
Address='C:\GLDAS\water\数据\COST-G\';
GFA = dir(fullfile(Address,'SW*'));
k=length(GFA);
GC=zeros((Lmax+3)*Lmax/2+1,k);
GS=zeros((Lmax+3)*Lmax/2+1,k);
DeltaGC6=zeros(Lmax+1,Lmax+1,k);
DeltaGS6=zeros(Lmax+1,Lmax+1,k);
DeltaRC6=zeros(Lmax+1,Lmax+1,k);
DeltaRS6=zeros(Lmax+1,Lmax+1,k);
for i=1:k
filepath=Address;
file=[filepath,GFA(i).name];
    fid=fopen(file);
if fid == -1
    'Error opening the file.';
end
k_line=0;
while 1
    
    nextline=fgetl(fid);      %读第一行
    if ~ischar(nextline)
        break,
    end %读到最后跳出
    if ~strcmp(nextline,'')
        if strcmp(nextline(1:3),'gfc')
            k_line1=str2double(nextline(6:8));
            if  k_line1>15
                break;
            end
            k_line=k_line+1;
            ll(k_line,i)=str2double(nextline(6:8));
            mm(k_line,i)=str2double(nextline(11:13));
            GC(k_line,i)=str2double(nextline(15:32));
            GS(k_line,i)=str2double(nextline(34:51));
        end
    end
end
sta=fclose(fid);
end
% C20_gra=read_SLR('C:\GLDAS\water\数据\TN-11_C20_SLR_RL06.txt');
% temp=importdata('C:\GLDAS\water\数据\TN-13_GEOC_CSR_RL06_2019-04-22.txt',' ',113);
% GC_Replace=temp.data;
% temp=[1:2,4:7,9:12,14:18,20:22,25:28,30:33,36:38,40:43];
% GC(3,temp)=C20_gra(131:163);
% k_line=0;
% k_line1=0;
% for i=1:size(GC_Replace,1)
%     if GC_Replace(i,2)==0
%         k_line=k_line+1;
%         GC10(k_line)=GC_Replace(i,3);
%         GS10(k_line)=GC_Replace(i,4);
%     elseif GC_Replace(i,2)==1
%         k_line1=k_line1+1;
%         GC11(k_line1)=GC_Replace(i,3);
%         GS11(k_line1)=GC_Replace(i,4);
%     end
% end
% GC(2,temp)=GC10(131:end);
% GS(2,temp)=GS10(131:end);
% GC(62,temp)=GC11(131:end);
% GS(62,temp)=GS11(131:end);
% Dgc=zeros((Lmax+3)*Lmax/2+1,k);
% Dgs=zeros((Lmax+3)*Lmax/2+1,k);
Dgc=zeros((Lmax+3)*Lmax/2+1,k);
Dgs=zeros((Lmax+3)*Lmax/2+1,k);

GC0=0;
for i=1:k
    GC0=GC0+GC(:,i);
end
GC00=GC0/k;

GS0=0;
for i=1:k
    GS0=GS0+GS(:,i);
end
GS00=GS0/k;
for i=1:k
    Dgc(:,i)=GC(:,i)-GC00;
    Dgs(:,i)=GS(:,i)-GS00;
end
for j=1:k
    for i=1:length(GC)
        DeltaGC6(ll(i,1)+1,mm(i,1)+1,j)=Dgc(i,j);
        DeltaGS6(ll(i,1)+1,mm(i,1)+1,j)=Dgs(i,j);
    end
end
clearvars -except DeltaGC DeltaGS DeltaRC DeltaRS DeltaGC1 DeltaGS1 DeltaRC1 DeltaRS1...
    DeltaGC2 DeltaGS2 DeltaRC2 DeltaRS2 DeltaGC3 DeltaGS3 DeltaRC3 DeltaRS3...
    DeltaGC4 DeltaGS4 DeltaRC4 DeltaRS4 DeltaGC5 DeltaGS5 DeltaRC5 DeltaRS5...
    DeltaGC6 DeltaGS6 DeltaRC6 DeltaRS6;


%%
Lmax=60;
Address='C:\GLDAS\water\数据\CSR-60阶\';
GFA = dir(fullfile(Address,'GSM*'));
k=length(GFA);
GC=zeros((Lmax+3)*Lmax/2+1,k);
GS=zeros((Lmax+3)*Lmax/2+1,k);
RC=zeros((Lmax+3)*Lmax/2+1,k);
RS=zeros((Lmax+3)*Lmax/2+1,k);
for i=1:k
    filepath=Address;
    file=[filepath,GFA(i).name];
    fid=fopen(file);
    if fid == -1
        'Error opening the file.';
    end 
    
    k_line=0;
    while 1

        nextline=fgetl(fid);      %读第一行
        if ~ischar(nextline)
            break,
        end %读到最后跳出
        if ~strcmp(nextline,'')
            if strcmp(nextline(1:6),'GRCOF2')
%         if length(nextline)==40      %笨方法，已知数据长度可以用
            k_line=k_line+1;
%            while k_line<=2||k_line==92
%                k_line=k_line+1;
%            end
           GC(k_line,i)=str2double(nextline(18:35));
           GS(k_line,i)=str2double(nextline(37:54));
           RC(k_line,i)=str2double(nextline(56:65));
           RS(k_line,i)=str2double(nextline(67:76));
            end
        end
    end
    sta=fclose(fid);
end
C20_gra=read_SLR('C:\GLDAS\water\datafile\TN-11_C20_SLR_RL06.txt');
GC(3,1:9)=C20_gra(164:end);

Lmax=60;
Address='C:\GLDAS\water\数据\60阶\';
% GIA=dir(fullfile(FAddress,'xx*));  %%xx表示路径,FAddress表示父文件夹
% yr=length(GIA);
% for t=1:yr
% Ffilepath=FAddress;
% Adress=[Ffilepath,GIA(t).name);
GFA = dir(fullfile(Address,'GSM*'));
k=length(GFA);
GC1=zeros((Lmax+3)*Lmax/2+1,k);
GS1=zeros((Lmax+3)*Lmax/2+1,k);
RC1=zeros((Lmax+3)*Lmax/2+1,k);
RS1=zeros((Lmax+3)*Lmax/2+1,k);
for i=1:k
    filepath=Address;
    file=[filepath,GFA(i).name];
    fid=fopen(file);
    if fid == -1
        'Error opening the file.';
    end
    
    k_line=0;
    while 1
        
        nextline=fgetl(fid);      %读第一行
        if ~ischar(nextline)
            break,
        end %读到最后跳出
        if ~strcmp(nextline,'')
            if strcmp(nextline(1:6),'GRCOF2')
                %         if length(nextline)==40      %笨方法，已知数据长度可以用
                k_line=k_line+1;
                %            while k_line<=2||k_line==92
                %                k_line=k_line+1;
                %            end
                GC1(k_line,i)=str2double(nextline(18:35));
                GS1(k_line,i)=str2double(nextline(37:54));
                RC1(k_line,i)=str2double(nextline(56:65));
                RS1(k_line,i)=str2double(nextline(67:76));
            end
        end
    end
    sta=fclose(fid);
end
%%%%三次样条插值求出缺失项
%x=[1:5,7:96,98:101,103:112,114:117,119:122,124:127,129:132,134:137,139:142,144:148,150:152,155:158,160:163];
%x=[1:2,4:7,9:12,14:18,20:22,25:28,30:33,36:38,40:43];
% [GC10,GS10,GC11,GS11]=read_deg1('C:\GLDAS\water\数据\deg1_coef.txt');
C20_gra=read_SLR('C:\GLDAS\water\数据\TN-11_C20_SLR_RL06.txt');
GC1(3,:)=C20_gra(1:163);
temp=importdata('C:\GLDAS\water\数据\TN-13_GEOC_CSR_RL06_2019-04-22.txt',' ',113);
GC_Replace=temp.data;
%%%%三次样条插值求出缺失项
k_line=0;
k_line1=0;
for i=1:size(GC_Replace,1)
    if GC_Replace(i,2)==0
        k_line=k_line+1;
        GC10(k_line)=GC_Replace(i,3);
        GS10(k_line)=GC_Replace(i,4);
    elseif GC_Replace(i,2)==1
        k_line1=k_line1+1;
        GC11(k_line1)=GC_Replace(i,3);
        GS11(k_line1)=GC_Replace(i,4);
    end
end
GC1(2,:)=GC10(1:end);
GS1(2,:)=GS10(1:end);
GC1(62,:)=GC11(1:end);
GS1(62,:)=GS11(1:end);

GCSUM=cat(2,GC,GC1);
GSSUM=cat(2,GS,GS1);
RCSUM=cat(2,RC,RC1);
RSSUM=cat(2,RS,RS1);

k=size(GCSUM,2);

DeltaGC7=zeros(Lmax+1,Lmax+1,k);
DeltaGS7=zeros(Lmax+1,Lmax+1,k);
DeltaRC7=zeros(Lmax+1,Lmax+1,k);
DeltaRS7=zeros(Lmax+1,Lmax+1,k);

for i=1:k
    for l=0:Lmax
        for m=0:l
            DeltaRC7(l+1,m+1,i)=RCSUM(l+1+(Lmax*2+1-m)*m/2,i);
            DeltaRS7(l+1,m+1,i)=RSSUM(l+1+(Lmax*2+1-m)*m/2,i);
        end
    end
end

Dgc=zeros((Lmax+3)*Lmax/2+1,k);
Dgs=zeros((Lmax+3)*Lmax/2+1,k);

GC0=0;
for i=1:k
    GC0=GC0+GCSUM(:,i);
end
GC00=GC0/k;

GS0=0;
for i=1:k
    GS0=GS0+GSSUM(:,i);
end
GS00=GS0/k;
for i=1:k
    Dgc(:,i)=GCSUM(:,i)-GC00;
    Dgs(:,i)=GSSUM(:,i)-GS00;
end
for i=1:k
    for l=0:Lmax
        for m=0:l
            DeltaGC7(l+1,m+1,i)=Dgc(l+1+(Lmax*2+1-m)*m/2,i);
            DeltaGS7(l+1,m+1,i)=Dgs(l+1+(Lmax*2+1-m)*m/2,i);
        end
    end
end

clearvars -except DeltaGC DeltaGS DeltaRC DeltaRS DeltaGC1 DeltaGS1 DeltaRC1 DeltaRS1...
    DeltaGC2 DeltaGS2 DeltaRC2 DeltaRS2 DeltaGC3 DeltaGS3 DeltaRC3 DeltaRS3...
    DeltaGC4 DeltaGS4 DeltaRC4 DeltaRS4 DeltaGC5 DeltaGS5 DeltaRC5 DeltaRS5...
    DeltaGC6 DeltaGS6 DeltaRC6 DeltaRS6 DeltaGC7 DeltaGS7 DeltaRC7 DeltaRS7;