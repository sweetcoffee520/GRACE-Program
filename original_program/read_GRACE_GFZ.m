%-----------------------------------------------
%Read GRACE data from CSR or GFZ;
%Input:
%        Address of the data
%Output:
%        DeltaC and DeltaS
%------------------------------------------------
function [DeltaGC,DeltaGS,k]=read_GRACE_GFZ(Address,Lmax)
%Lmax=120;
%Address='G:\GRACE\EIGEN_GFZ\EIGEN_GFZ_0\'
GFA = dir(fullfile(Address,'GSM*'));
k=length(GFA);
%k=1

GC=zeros((Lmax+3)*Lmax/2+1,k);
GS=zeros((Lmax+3)*Lmax/2+1,k);


Dgc=zeros((Lmax+3)*Lmax/2+1,k);
Dgs=zeros((Lmax+3)*Lmax/2+1,k);
DeltaGC=zeros(Lmax+1,Lmax+1,k);
DeltaGS=zeros(Lmax+1,Lmax+1,k);

for i=1:k
    filepath=Address;
    file=[filepath,GFA(i).name];
    fid=fopen(file);
    if fid == -1 
        ('Error opening the file.');
    end
    
    k_line=0;
    while 1
        nextline=fgetl(fid);      %读第一行
        if ~ischar(nextline)
            break,
        end %读到最后跳出
        
        
        if strcmp(nextline(1:6),'GRCOF2')
           k_line=k_line+1;
%            ll(k_line,i)=str2num(nextline(9:11));
%            mm(k_line,i)=str2num(nextline(14:16));
           GC(k_line,i)=str2double(nextline(18:35));
           GS(k_line,i)=str2double(nextline(37:54));
%            RC(k_line,i)=str2num(nextline(56:65));
%            RS(k_line,i)=str2num(nextline(67:76));
%            time_begin(k_line,i)=str2num(nextline(78:90));
%            time_end(k_line,i)=str2num(nextline(92:104));
        end
    end
end

%C20 replaced by C20 from SLR
%C20=textread('E:\data\C20.txt');
%GC(3,:)=C20';

%mean coefficients
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

%Caculate deltaC and deltaD
for i=1:k
    Dgc(:,i)=GC(:,i)-GC00;
    Dgs(:,i)=GS(:,i)-GS00;
end


%Change the Dgc(:,i) into Detagc(:,:,i)
% % for i=1:k
% % for l=0:Lmax
% %     for m=0:l
% %     DeltaGC(l+1,m+1,i)=Dgc(l+1+(Lmax*2+1-m)*m/2,i);
% %     DeltaGS(l+1,m+1,i)=Dgs(l+1+(Lmax*2+1-m)*m/2,i);
% %     end
% % end
% % end

for i=1:k
for l=0:60
    for m=0:l
    DeltaGC(l+1,m+1,i)=Dgc(l+1+(Lmax*2+1-m)*m/2,i);
    DeltaGS(l+1,m+1,i)=Dgs(l+1+(Lmax*2+1-m)*m/2,i);
    end
end
end