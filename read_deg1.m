function [GC10,GS10,GC11,GS11]=read_deg1(Fileaddress)
% Fileaddress='E:\GLDAS\water\datafile\deg1_coef.txt';
fid=fopen(Fileaddress);
    if fid == -1 
        ('Error opening the file'); 
    end
    
k_line=0;
k_line1=0;
while 1
    nextline=fgetl(fid);
    if ~ischar(nextline)
        break,
    end %读到最后跳出
    
    if ~strcmp(nextline,'')
    if strcmp(nextline(1:1),'2')
        if strcmp(nextline(16:16),'0')
            k_line=k_line+1;
        GC10(k_line,1)=str2double(nextline(17:35));
        GS10(k_line,1)=str2double(nextline(36:54));
        elseif strcmp(nextline(16:16),'1')
            k_line1=k_line1+1;
        GC11(k_line1,1)=str2double(nextline(17:35));
        GS11(k_line1,1)=str2double(nextline(36:54));
        end
    end 
    end
end
fclose(fid);

% TimeC20=zeros(length(C20),3);
% for i=1:length(C20)
%     [TimeC20(i,1),TimeC20(i,2),TimeC20(i,3)]=MJD2YMD(MJD_C20(i,1));
% end

% C20_gra=zeros(length(Data_time(:,1)),1);
% for i=1:length(C20)
%     for j=1:length(Data_time(:,1))
%         if TimeC20(i,1)==Data_time(j,1)&&TimeC20(i,2)==Data_time(j,2);
%             C20_gra(j,1)=C20(i,1);
%         end
%     end
% end

