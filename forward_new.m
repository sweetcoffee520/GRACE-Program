function grid=forward_new(m_grace)
load('lg_msk_land_025_new.mat');
land_index=lg_msk_land;
load('lg_msk_ocean_025_new.mat');
ocean_index=lg_msk_ocean;

%load('region_ocean.mat');
%region_index=region_data;
%region_index=read_msk_25('bohai.bln');
%other_index=replace_bln(land_index);
m0=m_grace;

for i=1:45
       
  cs=gmt_grid2cs(m0,60);
%cs(1,1)=0;
%cs(1,2)=0;
%cs(2,1)=0;
%cs(2,2)=0;

   grid2=gmt_cs2grid(cs,200,0.25);
   %deta=m_grace-grid2;
  land_deta=(m_grace-grid2).*land_index;
  % other_deta=grid2.*other_index;
  a=sum(sum(land_deta));
  b=-a/686581;  %海陆质量平衡686581为海洋格网数
ocean_deta=(m_grace-grid2).*ocean_index;
  ocean_new=ocean_deta+b*ocean_index;
   m0=(m0+land_deta+ocean_new);

end
  tmp=m0;
 grid=tmp;
end
%grid=m0;
