function [Dc_w,Ds_w]=Gau4grace(radius,Lmax,Dc,Ds,data_number)
%Gaussian Filter-------------------------------------------------------
% % r=300;
% % Lmax=60;
w1=zeros(1,Lmax+1);

a=6.378136460E+06;% µ¥Î»m
r1=radius*1000;%%Unit:km
b1=log(2)/(1-cos(r1/a));
w1(1)=1;

w1(2)=(1+exp(-2*b1))/(1-exp(-2*b1))-1/b1;
for l=1:(Lmax-1)
w1(l+2)=-(2*l+1)/b1*w1(l+1)+w1(l);
end


Dc_w=zeros(Lmax+1,Lmax+1,data_number);
Ds_w=zeros(Lmax+1,Lmax+1,data_number);
for i=1:data_number
    for j=1:Lmax+1
        Dc_w(:,j,i)=Dc(:,j,i).*w1';
        Ds_w(:,j,i)=Ds(:,j,i).*w1';
    end
end

        