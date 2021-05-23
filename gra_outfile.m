function gra_outfile(data,filename,format,outsize)
%
% gra_outfile(data,filename,format,outsize)
% 
% +Purpose:
%    write out truncated simulated data in given format
%
% Created by Feng,W.P., @ GU, 2012-04-09
%
%
if nargin < 1
    disp('Usage: gra_outfile(data,filename,format,outsize)');
    disp('  data,     a matrix of n*3, n is the number of points and 3 for lon,lat and gravity');
    disp('  filename, the root name of output file, special postfices have been given in the code ');
    disp('            .ascii for ascii file, .bin for binary foramt and .gid for grided output');
    disp('  format,   3 types of formats can be supported: ascii, binary and grided file with rsc header');
    disp('  outsize,  if format is grided binary, a outsize will be requested.');
    disp(' +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++');
    disp('  Created by Feng, W.P., @ GU, 2012-09-04');
    disp(' +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++');
    return
end
%
switch upper(format)
    case 'ASCII'
        outfile = [filename,'.ascii'];
        fid = fopen(outfile,'w');
        fprintf(fid,'%f %f %f\n',data');
        fclose(fid);
    case 'BINARY'
        outfile = [filename,'.bin'];
        fid     = fopen(outfile,'w');
        fwrite(fid,data','float');
        fclose(fid);
    case 'GRIDED'
        if numel(outsize)<2
            outsize(2) = outsize(1);
        end
        %
        outfile = [filename,'.grid'];
        minlon  = min(data(:,1));
        maxlon  = max(data(:,1));
        minlat  = min(data(:,2));
        maxlat  = max(data(:,2));
        %[lon,lat] = meshgrid(minlon:outsize(1):maxlon,minlat:outsize(1):maxlat);
        len       = fix((maxlon-minlon)/outsize(1))+1;
        wid       = fix((maxlat-minlat)/outsize(2))+1;
        outimage  = zeros(wid,len);
        clon      = fix((data(:,1)-minlon)/outsize(1))+1;
        clat      = fix((data(:,2)-minlat)/outsize(2))+1;
        ind           = sub2ind(size(outimage),clat(:),clon(:));
        outimage(ind) = data(:,3);
        %
        fid = fopen(outfile,'w');
        fwrite(fid,fliplr(outimage'),'float');
        fclose(fid);
        info         = sim_roirsc();
        info.x_first = minlon;
        info.y_first = maxlat;
        info.x_step  = outsize(1);
        info.y_step  = -1*outsize(2);
        info.width   = len;
        info.file_length = wid;
        %
        outrsc = [outfile,'.rsc'];
        sim_croirsc(outrsc,info);
        %
        
end