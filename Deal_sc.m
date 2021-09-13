function [WSD_fa,WSD_gri,ave_WSD_fa,ave_WSD_gri,WSDI_fa,WSDI_gri,ave_WSDI_fa,ave_WSDI_gri]=Deal_sc(ave_vary_fa_a,vary_fa_a,ave_griddata_sc,griddatasc)
[n,o,m]=size(griddatasc);
temp=zeros(1,163);
%用相邻两个月份求平均填充缺失的等效水高
temp([1:5,7:96,98:101,103:112,114:117,119:122,124:127,129:132,134:137,139:142,144:148,150:152,155:158,160:163])=ave_vary_fa_a;
temp([153,154])=ave_griddata_sc([153,154]);   %由于中间连续缺失两个月，所以用水文数据替代
temp([6,97,102,113,118,123,128,133,138,143,149,159])=(temp([6,97,102,113,118,123,128,133,138,143,149,159]-1)+temp([6,97,102,113,118,123,128,133,138,143,149,159]+1))/2;
ave_vary_fa_a=temp;
%用相邻两个月份求平均填充缺失的月份经纬值
temp1=zeros(n,o,m);
temp1(:,:,[1:5,7:96,98:101,103:112,114:117,119:122,124:127,129:132,134:137,139:142,144:148,150:152,155:158,160:163])=vary_fa_a;
temp1(:,:,[153,154])=griddatasc(:,:,[153,154]); %由于中间连续缺失两个月，所以用水文数据替代
temp1(:,:,[6,97,102,113,118,123,128,133,138,143,149,159])=(temp1(:,:,[6,97,102,113,118,123,128,133,138,143,149,159]-1)+temp1(:,:,[6,97,102,113,118,123,128,133,138,143,149,159]+1))/2;
vary_fa_a=temp1;

for i=1:12
    x=i:12:156;
    WSD_fa(:,:,x)=vary_fa_a(:,3,x)-sum(vary_fa_a(:,3,x),3)/length(x);
    WSD_fax(:,:,x)=[vary_fa_a(:,1:2,x),WSD_fa(:,1,x)];
    WSD_gri(:,:,x)=griddatasc(:,3,x)-sum(griddatasc(:,3,x),3)/length(x);
    WSD_grix(:,:,x)=[griddatasc(:,1:2,x),WSD_gri(:,1,x)];
    ave_WSD_fa(x)=ave_vary_fa_a(x)-sum(ave_vary_fa_a(x))/length(x);
    ave_WSD_gri(x)=ave_griddata_sc(x)-sum(ave_griddata_sc(x))/length(x);
end


%计算WSDI
WSDI_fa=(WSD_fa(:,1,:)-mean(WSD_fa(:,1,:),3))./std(WSD_fa(:,1,:),1,3);
WSDI_gri=(WSD_gri(:,1,:)-mean(WSD_gri(:,1,:),3))./std(WSD_gri(:,1,:),1,3);
ave_WSDI_fa=(ave_WSD_fa-mean(ave_WSD_fa))/std(ave_WSD_fa);
ave_WSDI_gri=(ave_WSD_gri-mean(ave_WSD_gri))/std(ave_WSD_gri);

for r=1:12
    M=int2str(r);
    N='.txt';
    L='Datasc';
    output=[L,M,N];
    fid=fopen(['C:\GLDAS\water\WSD\',output],'w');
    if fid == -1
        'Error opening the file.';
    end
    
    fprintf(fid,'%-5.1f %-4.1f %8.6f\n',WSD_fax(:,:,83+r)');
    %          fprintf(fid,'%010.5f\n',vary_f3(:,r));
    
    sta=fclose(fid);
end    