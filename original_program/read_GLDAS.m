%-----------------------------------------------
%Read GLDAS data ;
%Input:
%        Address of the data
%Output:
%        DeltaC and DeltaS
%------------------------------------------------
function [DeltaGC,DeltaGS,k]=read_GLDAS(Address,Lmax)

GFA = dir(fullfile(Address,'*.txt'));
k=length(GFA);
for i=1:k
        GFA(i).name=strcat(Address,GFA(i).name);
    [ll(:,i) mm(:,i) GC(:,i) GS(:,i)]=textread(GFA(i).name,'%f %f %f %f');
end

%C20 replaced by C20 from SLR
%C20=textread('E:\data\C20.txt');
%GC(3,:)=C20';

%mean coefficients
% GC0=0;
% for i=1:k
%     GC0=GC0+GC(:,i);
% end
% GC00=GC0/k;
% 
% GS0=0;
% for i=1:k
%     GS0=GS0+GS(:,i);
% end
% GS00=GS0/k;
% 
% %Caculate deltaC and deltaD
% for i=1:k
%     Dgc(:,i)=GC(:,i)-GC00;
%     Dgs(:,i)=GS(:,i)-GS00;
% end


%Change the Dgc(:,i) into Detagc(:,:,i)
% % for i=1:k
% % for l=0:Lmax
% %     for m=0:l
% %     DeltaGC(l+1,m+1,i)=Dgc(l+1+(Lmax*2+1-m)*m/2,i);
% %     DeltaGS(l+1,m+1,i)=Dgs(l+1+(Lmax*2+1-m)*m/2,i);
% %     DeltaGC(l+1,m+1,i)=GC(l+1+(Lmax*2+1-m)*m/2,i);
% %     DeltaGS(l+1,m+1,i)=GS(l+1+(Lmax*2+1-m)*m/2,i);
% %     end
% % end
% % end

for i=1:k
    for j=1:length(ll(:,i))
        DeltaGC(ll(j)+1,mm(j)+1)=GC(j);
        DeltaGS(ll(j)+1,mm(j)+1)=GS(j);
    end
end

        
