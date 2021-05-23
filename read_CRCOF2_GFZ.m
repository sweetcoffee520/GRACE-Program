function [DeltaGC,DeltaGS,k]=read_CRCOF2_GFZ(Address,Lmax)
Lmax=60;
Address='C:\GLDAS\water\datafile\GSM\';
% GIA=dir(fullfile(FAddress,'xx*));  %%xx表示路径,FAddress表示父文件夹
% yr=length(GIA);
% for t=1:yr
% Ffilepath=FAddress;
% Adress=[Ffilepath,GIA(t).name);
GFA = dir(fullfile(Address,'GSM*'));
k=length(GFA);
GC=zeros((Lmax+3)*Lmax/2+1,k);
GS=zeros((Lmax+3)*Lmax/2+1,k);
DeltaGC=zeros(Lmax+1,Lmax+1,k);
DeltaGS=zeros(Lmax+1,Lmax+1,k);
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
           RC(k_line,i)=str2num(nextline(56:65));
           RS(k_line,i)=str2num(nextline(67:76));
            end
        end
    end
    sta=fclose(fid);
end

C20_gra=read_SLR('C:\GLDAS\water\datafile\C20_RL06.txt');
[GC10,GS10,GC11,GS11]=read_deg1('C:\GLDAS\water\datafile\deg1_coef.txt');
GC(3,:)=C20_gra([13:17,19:108,110:113,115:124,126:129,131:134,136:139,141:144,146:149,151:154,156:160,162:164,167:170,172:175]);
% C20_gra([25:37,55:85]) %提取数据
% GC(3,:)=C20_gra(109:120,:);
GC(2,[1:11,13:end])=GC10(6:153);
GS(2,[1:11,13:end])=GS10(6:153);
GC(62,[1:11,13:end])=GC11(6:153);
GS(62,[1:11,13:end])=GS11(6:153);
% GS(2,:)=GS10(88:99);
% GC(62,:)=GC11(88:99);
% GS(62,:)=GS11(88:99);

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
            DeltaRC(l+1,m+1,i)=RC(l+1+(Lmax*2+1-m)*m/2,i);
            DeltaRS(l+1,m+1,i)=RS(l+1+(Lmax*2+1-m)*m/2,i);
        end
    end
end

L=1./(2*(1:Lmax+1)+1);
G=sum(DeltaRC.*DeltaRC+DeltaRS.*DeltaRS,2);
for i=1:k
sgm(i,:)=sqrt(L.*G(:,1,i)');
end
