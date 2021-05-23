function Time_GLDAS=GrabTheDate_GLDAS_nc4(Address)

list = dir(fullfile(Address,'*.nc4'));
k=length(list);
Time_GLDAS=zeros(k,2);
for i=1:k
    liste=char(list(i).name);
    line_length=length(liste);
    Time_GLDAS(i,1)=str2double(liste(line_length-13:line_length-10));
    Time_GLDAS(i,2)=str2double(liste(line_length-9:line_length-8));
end