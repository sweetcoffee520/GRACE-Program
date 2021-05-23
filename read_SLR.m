% %Read SLR C20
function [C20_gra]=read_SLR(Fileaddress)
% Fileaddress='C:\GLDAS\water\datafile\C20_RL06.txt';
fid=fopen(Fileaddress);
    if fid == -1 
        ('Error opening the file'); 
    end
    
k_line=0;

while 1
    nextline=fgetl(fid);
    if ~ischar(nextline)
        break,
    end %读到最后跳出
    
        if ~strcmp(nextline,'#')
    if ~strcmp(nextline,'')
%             if strcmp(nextline(2:2),'2')
        if strcmp(nextline(1:1),'5')
            k_line=k_line+1;
            %         time(k_line,1)=str2double(nextline(2:10));
            C20_gra(k_line,1)=str2double(nextline(20:39));
        end
    end
        end
end
fclose(fid);

