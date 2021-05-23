function  [a,b,c,d,e,f,g,h]=region_grid(c1,c2,f1,f2,res_lonlat)
ceta=c1:res_lonlat:c2;n_c=length(ceta);
fir=f1:res_lonlat:f2;n_f=length(fir);%China
for i=1:n_c
    for j=1:n_f
        firx(j+(i-1)*n_f,1)=fir(j);
        nfir(j+(i-1)*n_f,1)=j;
        cetax(j+(i-1)*n_f,1)=ceta(i);
        nceta(j+(i-1)*n_f,1)=i;
    end
end 

a=ceta; %%纬度矩阵
b=fir;  %%经度矩阵
c=n_c;  %%纬度矩阵长度
d=n_f;  %%经度矩阵长度
e=cetax;
f=firx;
g=nceta;
h=nfir;