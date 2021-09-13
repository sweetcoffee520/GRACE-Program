function [GRACE_all]=CS2Grav(Data_type,DeltaGC,DeltaGS,Data_number,Nlmx,De_filter,Filter_index,Lmax,De_P,De_M,r_guassian,r1_fan,r2_fan,fir,n_c,n_f,nceta,nfir)

k=Data_number;
%%---------------------------------------------------------------------
if isequal(De_filter,'1')
    [DeltaGC_0,DeltaGS_0]=Dec4grace(Lmax,DeltaGC,DeltaGS,k,De_P,De_M);
    disp('Decorrelated coefficients  produced!')
else 
    Filter_index(4:6)={'0'};
    disp('Error. Without Decorrelated coefficients!')
end

%%-------------------------------------------------------------------------
GRACE_all(length(nceta),k,6)=0;

if isequal(Filter_index{1},'1')
    [grace_non]=grace2gravity(Lmax,k,Nlmx,DeltaGC,DeltaGS,fir,n_c,n_f,nceta,nfir);
    disp('Results without filter done!')
    GRACE_all(:,:,1)=grace_non;
end
%%---------------------------------
if isequal(Filter_index{2},'1')
    [DeltaGC_g,DeltaGS_g]=Gau4grace(r_guassian,Lmax,DeltaGC,DeltaGS,k);
    [grace_gau]=grace2gravity(Lmax,k,Nlmx,DeltaGC_g,DeltaGS_g,fir,n_c,n_f,nceta,nfir);
    disp('Results with Gaussian filter done!')
    GRACE_all(:,:,2)=grace_gau;
end
%%---------------------------------
if isequal(Filter_index{3},'1')
    [DeltaGC_f,DeltaGS_f]=Fan4grace(r1_fan,r2_fan,Lmax,DeltaGC,DeltaGS,k);
    [grace_fan]=grace2gravity(Lmax,k,Nlmx,DeltaGC_f,DeltaGS_f,fir,n_c,n_f,nceta,nfir);
    disp('Results with Fan filter done!')
    GRACE_all(:,:,3)=grace_fan;
end
%%--------------------------------------------%%
if isequal(Filter_index{4},'1')
    [grace_p3m6_non]=grace2gravity(Lmax,k,Nlmx,DeltaGC_0,DeltaGS_0,fir,n_c,n_f,nceta,nfir);
    disp('Results with decorrlation P3M6 done!')
    GRACE_all(:,:,4)=grace_p3m6_non;
end
%%---------------------------------
if isequal(Filter_index{5},'1')
    [DeltaGC_g,DeltaGS_g]=Gau4grace(r_guassian,Lmax,DeltaGC_0,DeltaGS_0,k);
    [grace_p3m6_gau]=grace2gravity(Lmax,k,Nlmx,DeltaGC_g,DeltaGS_g,fir,n_c,n_f,nceta,nfir);
    disp('Results with P3M6 and Gaussian filter done!')
    GRACE_all(:,:,5)=grace_p3m6_gau;
end
%%------------------------------------------%%
if isequal(Filter_index{6},'1')
    [DeltaGC_f,DeltaGS_f]=Fan4grace(r1_fan,r2_fan,Lmax,DeltaGC_0,DeltaGS_0,k);
    [grace_p3m6_fan]=grace2gravity(Lmax,k,Nlmx,DeltaGC_f,DeltaGS_f,fir,n_c,n_f,nceta,nfir);
    disp('Results with P3M6 and Fan filter done!')
    GRACE_all(:,:,6)=grace_p3m6_fan;
end

if isequal(Data_type,'GLDAS')
    GRACE_all=GRACE_all*10;
end