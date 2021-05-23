% function [DeltaGC,DeltaGS,k]=read_CRCOF2_CSR(Address,Lmax)
Lmax=96;
Address='C:\GLDAS\water\datafile\CSR-96��\';
% GIA=dir(fullfile(FAddress,'xx*));  %%xx��ʾ·��,FAddress��ʾ���ļ���
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

        nextline=fgetl(fid);      %����һ��
        if ~ischar(nextline)
            break,
        end %�����������
        if ~strcmp(nextline,'')
            if strcmp(nextline(1:6),'GRCOF2')
%         if length(nextline)==40      %����������֪���ݳ��ȿ�����
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

C20_gra=read_SLR('C:\GLDAS\water\datafile\TN-11_C20_SLR_RL06.txt');
GC(3,:)=C20_gra(164:end);

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