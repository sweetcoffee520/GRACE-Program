function [D_Deltegc,D_Deltegs] = Dec4grace(max_deg,GC,GS,k,p,startorder,startdegree,endorder)
% max_deg=60;
% startdegree=6;
% p=4;
% m=6;
% startorder=6;
% endorder=50;
%%
% 在Q.Li. 2009.12.8基础上修改
%************************
% max_deg, the maximum of degree for legendre function
% GC, the coefficients
% GS, the coefficients
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% max_deg, the maximum of degree
% k,       version, default 1
% p,       order for polynominal fitting
% m,       the order for fitting
%% 
if nargin < 6 || isempty(startorder)
    startorder = 6;
end

if nargin < 7 || isempty(startdegree)
    startdegree = 11;
end

if nargin <8 || isempty(endorder)
    endorder = 50;
end
%
D_Deltegc = GC;
D_Deltegs = GS;
% D_Deltegc1=zeros(max_deg+1,max_deg+1,k);
% D_Deltegs1=zeros(max_deg+1,max_deg+1,k);

%% 
% Improved by Feng,W.P., @ GU, 2012-09-3
% -> startorder, startdegree and endorder work in the latest version.
if startdegree >= max_deg
    startdegree = max_deg - 1;
end
if startorder > startdegree
    startorder = startdegree - 1;
end
if endorder <= startorder
    endorder = startorder + 1;
end
%%
lown  = startdegree;  % start degree, 11 in default
lowm  = startorder;   % start order, 6
upm   = endorder;     % end order, 50
order = p ;           % order for fitting function
%
for j=1:k
    D_gc(:,:) = GC(:,:,j);
    D_gs(:,:) = GS(:,:,j);
    for m = lowm:upm
        if(m<=lown)
            ms = lown;
            if(mod(lown,2)~=0)
                fitptn = (max_deg-lown+1)/2;
                cs_x   = zeros(1,fitptn);
                bias   = 1;
                Coe    = coe_fit(m,ms,fitptn,D_gc,D_gs,bias,order);
                for i=1:fitptn
                    cs_x(i)                    = ms+2*i-bias;
                    D_Deltegc(cs_x(i),m+1,j) = Coe(cs_x(i),m+1,1);
                    D_Deltegs(cs_x(i),m+1,j) = Coe(cs_x(i),m+1,2);
                end
                fitptn = (max_deg-lown+1)/2;
                cs_x   = zeros(1,fitptn);
                bias   = 0;
                Coe    = coe_fit(m,ms,fitptn,D_gc,D_gs,bias,order);
                for i=1:fitptn
                    cs_x(i)                    = ms+2*i-bias;
                    D_Deltegc(cs_x(i),m+1,j) = Coe(cs_x(i),m+1,1);
                    D_Deltegs(cs_x(i),m+1,j) = Coe(cs_x(i),m+1,2);
                end
            end
            %
            if(mod(lown,2)==0)
                fitptn = (max_deg-lown)/2+1;
                cs_x   = zeros(1,fitptn);
                bias   = 1;
                Coe    = coe_fit(m,ms,fitptn,D_gc,D_gs,bias,order);
                for i=1:fitptn
                    cs_x(i)                    = ms+2*i-bias;
                    D_Deltegc(cs_x(i),m+1,j) = Coe(cs_x(i),m+1,1);
                    D_Deltegs(cs_x(i),m+1,j) = Coe(cs_x(i),m+1,2);
                end
                fitptn = (max_deg-lown)/2;
                cs_x   = zeros(1,fitptn);
                bias   = 0;
                Coe    = coe_fit(m,ms,fitptn,D_gc,D_gs,bias,order);
                for i=1:fitptn
                    cs_x(i)                    = ms+2*i-bias;
                    D_Deltegc(cs_x(i),m+1,j) = Coe(cs_x(i),m+1,1);
                    D_Deltegs(cs_x(i),m+1,j) = Coe(cs_x(i),m+1,2);
                end
            end
        else
            ms = m;
            if(mod(m,2)==0)
                fitptn = (max_deg-m)/2+1;
                cs_x   = zeros(1,fitptn);
                bias   = 1;
                Coe    = coe_fit(m,ms,fitptn,D_gc,D_gs,bias,order);
                for i=1:fitptn
                    cs_x(i)                    = ms+2*i-bias;
                    D_Deltegc(cs_x(i),m+1,j) = Coe(cs_x(i),m+1,1);
                    D_Deltegs(cs_x(i),m+1,j) = Coe(cs_x(i),m+1,2);
                end
                fitptn = (max_deg-m)/2;
                bias   = 0;
                Coe    = coe_fit(m,ms,fitptn,D_gc,D_gs,bias,order);
                for i=1:fitptn
                    cs_x(i)                    = ms+2*i-bias;
                    D_Deltegc(cs_x(i),m+1,j) = Coe(cs_x(i),m+1,1);
                    D_Deltegs(cs_x(i),m+1,j) = Coe(cs_x(i),m+1,2);
                end
            else
                fitptn = (max_deg-m+1)/2;
                cs_x   = zeros(1,fitptn);
                bias   = 1;
                Coe    = coe_fit(m,ms,fitptn,D_gc,D_gs,bias,order);
                for i=1:fitptn
                    cs_x(i)                    = ms+2*i-bias;
                    D_Deltegc(cs_x(i),m+1,j) = Coe(cs_x(i),m+1,1);
                    D_Deltegs(cs_x(i),m+1,j) = Coe(cs_x(i),m+1,2);
                end
                fitptn = (max_deg-m+1)/2;
                bias   = 0;
                Coe    = coe_fit(m,ms,fitptn,D_gc,D_gs,bias,order);
                for i=1:fitptn
                    cs_x(i)                    = ms+2*i-bias;
                    D_Deltegc(cs_x(i),m+1,j) = Coe(cs_x(i),m+1,1);
                    D_Deltegs(cs_x(i),m+1,j) = Coe(cs_x(i),m+1,2);
                end
            end
        end
    end
end
%%

function Coe = coe_fit(m,ms,fitptn,Detagc,Detags,bias,order)
%****************************
%
% 在Q.Li. 2009.12.8基础上修改
%****************************
cs_x = zeros(fitptn,1);
cy   = cs_x;
sy   = cs_x;
Cec  = cs_x;
Ces  = cs_x;
Coe  = cs_x;
%
for i=1:fitptn
    cs_x(i)=ms+2*i-bias;
    cy(i)=Detagc(cs_x(i),m+1);
    sy(i)=Detags(cs_x(i),m+1);
end
%
pc = polyfit(cs_x,cy,order);
ps = polyfit(cs_x,sy,order);
%
for i=1:fitptn
    Cec(i)               = polyval(pc,cs_x(i));
    Ces(i)               = polyval(ps,cs_x(i));
    Coe(cs_x(i),m+1,1) = Detagc(cs_x(i),m+1)-Cec(i);
    Coe(cs_x(i),m+1,2) = Detags(cs_x(i),m+1)-Ces(i);
end
%%