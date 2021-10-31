% function rP=legendre_n(x,Nmax)
% tested correct!
%{
# input: 
x : cos(theta)
Nmax : the max degree to compute

# output:
rP(1:Nmax+1, 1:Nmax+1) will return

# Recurrence formula:
P(n,m) = A*p(n-1,m) + B*p(n-2,m)
%}
x = cosd(90+88);
Nmax=60;
rP(1:Nmax+1,1:Nmax+1)=0;

if abs(x)>1.0 
    error 'ERROR in legendre: x should <= 1.0'
end

%m=n=0
rP(1,1)=1.0;
%n=1,m=0         
rP(2,1)=sqrt(3.0)*x;
%n=m=1         
rP(2,2)=sqrt(3.0)*sqrt(1.0-x*x);
%m=n
for in=3:Nmax+1
    rP(in,in)=sqrt( (1.0-x*x) * (2*in-1) / (2*in-2) )*rP(in-1,in-1);
end
%n-m=1
for in=3:Nmax+1
    rP(in,in-1)=x*sqrt(2*in-1)*rP(in-1,in-1);
end
%others
for im=1:Nmax-1
    for in=im+2:Nmax+1
        rP(in,im)=x*sqrt((2.0*in-3)*(2.0*in-1)/((in-im)*(in+im-2)))*rP(in-1,im) ...
            -sqrt((2.0*in-1)*(in+im-3.0)*(in-im-1.0)/((2.0*in-5)*(in+im-2)*(in-im)))*rP(in-2,im);
    end
end