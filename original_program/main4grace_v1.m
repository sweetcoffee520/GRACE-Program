% Parameter_address='G:\GRACE\RL05\result\parametre_grace_v4.txt';
% % Parameter_address         =   'G:\Parametre\China.txt';
% % Parameter_address         =   'G:\Parametre\MekongRiver.txt';

% % Input_GRACE_data_address  =   'E:\Programs and data\GLDAS data\GEO_GLDAS0\';
% % Output_data_address       =   'G:\Results\GLDAS\';

Parameter_address         =   'G:\Parametre\China_v1.txt';
% % Input_GRACE_data_address  =   'G:\GRACE\CSR_RL05_2004_2010\';
% % Output_data_address       =   'G:\Results\RL04 vs RL05\';

C20_address               =   'G:\C20\TN-07_C20_SLR.txt';
C20_Flag                  =   1; %1-default replace by SLR C20

% % Input_GRACE_data_address  =   'G:\DATA\GRACE\EIGEN_GFZ\EIGEN_GFZ_120\';
% % Output_data_address       =   'G:\Results\';


%%Load the parameter data
getparameter_1

%Make grid for research region
[ceta,fir,n_c,n_f,cetax,firx,nceta,nfir]=region_grid(minlat,maxlat,minlon,maxlon,Res_lonlat);
disp('Region Gridded is ready!')

%Legendre polynom
Pnm=Nlmx_v3(Lmax,ceta);
disp('Legendre polynom is ready!')


% % Read GRACE or GLDAS data
[DeltaGC,DeltaGS,Data_number,Data_time]=readdata(Data_type,Input_GRACE_data_address,Lmax,C20_address,C20_Flag);



switch upper(Result_type)
    case 'GRAVITY'
        [GRACE_all]=CS2Grav(Data_type,DeltaGC,DeltaGS,Data_number,Pnm,De_filter,Filter_index,Lmax,De_P,De_M,Gaussian_r,Fan_r1,Fan_r2,fir,n_c,n_f,nceta,nfir);
    case 'EWH'
        [GRACE_all]=CS2EWH(Data_type,DeltaGC,DeltaGS,Data_number,Pnm,De_filter,Filter_index,Lmax,De_P,De_M,Gaussian_r,Fan_r1,Fan_r2,fir,n_c,n_f,nceta,nfir);
    case 'DISPLACEMENT'
        [GRACE_all]=CS2Disp(Data_type,DeltaGC,DeltaGS,Data_number,Pnm,De_filter,Filter_index,Lmax,De_P,De_M,Gaussian_r,Fan_r1,Fan_r2,fir,n_c,n_f,nceta,nfir);
end

% save data fo all monthly GRACE data
if isequal(Save_monthly_data,'0')
    for i=1:6
        if isequal(Filter_index{i},'1')
        for j=1:Data_number
                filename = strcat(Output_data_address,Keyword,'_',Data_type,'_',Filter_type{i},'_',Result_type,'_',num2str(Data_time(j,1)),num2str(Data_time(j,2),'%02d'));
                gra_outfile([firx,cetax,GRACE_all(:,j,i)],filename,'ascii',Res_lonlat)
        end
        end
    end
end  