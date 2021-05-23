function [Dc_w2,Ds_w2]=Fan4grace(radius1,radius2,Lmax,Dc,Ds,data_number)



a=6.378136460E+06;% µ¥Î»m

w1=zeros(1,Lmax+1);
r1=radius1*1000;%%Unit:km
b1=log(2)/(1-cos(r1/a));
w1(1)=1;
w1(2)=(1+exp(-2*b1))/(1-exp(-2*b1))-1/b1;
for l=1:(Lmax-1)
w1(l+2)=-(2*l+1)/b1*w1(l+1)+w1(l);
end

r2=radius2*1000;
w2=zeros(1,Lmax+1);
b2=log(2)/(1-cos(r2/a));
w2(1)=1;
w2(2)=(1+exp(-2*b2))/(1-exp(-2*b2))-1/b2;
for l=1:(Lmax-1)
w2(l+2)=-(2*l+1)/b2*w2(l+1)+w2(l);
end

% Fan_w=w1*w2';

Dc_w1=zeros(Lmax+1,Lmax+1,data_number);
Ds_w1=zeros(Lmax+1,Lmax+1,data_number);
Dc_w2=zeros(Lmax+1,Lmax+1,data_number);
Ds_w2=zeros(Lmax+1,Lmax+1,data_number);
for i=1:data_number
    for j=1:Lmax+1
        Dc_w1(:,j,i)=Dc(:,j,i).*w1';
        Ds_w1(:,j,i)=Ds(:,j,i).*w1';
    end
end

for i=1:data_number
    for j=1:Lmax+1
        Dc_w2(j,:,i)=Dc_w1(j,:,i).*w2;
        Ds_w2(j,:,i)=Ds_w1(j,:,i).*w2;
    end
end
        