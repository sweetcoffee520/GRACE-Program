 function mNLegendre = Nlmx_v3(Lmax,x)
%
%
% Lmax=60;x=53;
%
% Developed by Qiong Li, University of Wuhan, 08-2012
% mNLegendre is a matrix for saving the coefficients of Legendre
% function...
% Improved by Feng,W.P. @ GU, 2012-08-09
%   -> support matrices computation ...
%   -> speed has been improved over 20 times...
%
%
 kk    = length(x);
Theta = 90-x;
%//
% mNLegendre = zeros(Lmax,Lmax,kk);
mNLegendre = zeros(Lmax,Lmax,kk);
mx         = zeros(Lmax+1,1);       % the length should be greater or equal to m+1.
mx(1,1)    = 1;
for i = 1:kk
for l = 0:Lmax
    %
    cy = legendre(l,cosd(Theta),'norm');
    %
    for m=0:l
        mNLegendre(l+1,m+1,:)= Nlegendre_v3(l,m,cy,mx);
%           mNLegendre= Nlegendre_v3(l,m,cy,mx);
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = Nlegendre_v3(n,m,cy,mx)
% 
% y = Nlegendre(n,m,cy,model)
y  = cy'*mzeros_v3(n+1,m+1)*sqrt(2*(2-mx(m+1,1)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = mzeros_v3(n,m)
% Developed by Li Qiong, which will be used in Nlegendre()
y      = zeros(n,1);
y(m,1) = 1;