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

a=ceta;
b=fir;
c=n_c;
d=n_f;
e=cetax;
f=firx;
g=nceta;
h=nfir;