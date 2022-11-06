function [Tn,dTndx,d2Tndx2]=f_EvalChebyshevPolynomials(n,x)
%
% function [Tn,dTndx,d2Tndx2]=f_EvalChebyshevPolynomials(n,x)
% Version 1.0
%
% This function generates the Chebyshev polynomials of degree n in points
% x. Both n and x are intended to be row vectors, with different
% dimensions, where the "usual" normalization is adopted.
%
% The function is based on the goniometric definitions of Chebyshev
% polynomials, the first derivative is computed through the second kind
% polynomials, and the second derivative is evaluated with by using the
% differential equation. The values of the function and
% of the derivatives at the [-1,1] interval endpoints are computed
% analytically.
%
% Alberto Tibaldi, 19/01/2016

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finding internal points, x=+1 and x=-1 points (for analytical limits)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
threshold=1e-11; % threshold to identify internal points 
ind=find(abs(x-1)>threshold & abs(x-(-1))>threshold); % internal points
ind_p1=find(abs(x-1)<threshold); % x=+1
ind_m1=find(abs(x-(-1))<threshold); % x=-1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Synthesis of Chebyshev polynomials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tn=ones(length(n),length(x)); % P_0(x) = 1
if(not(isempty(ind)))
    thx=acos(x);
    [THx,dum]=meshgrid(thx,n);
    [X,N]=meshgrid(x,n);
    Tn=cos(N.*THx);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Synthesis of Chebyshev polynomials first derivative
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dTndx=zeros(length(n),length(x));
if(not(isempty(ind)))
    Un=sin((N+1).*THx)./sin(THx);
    dTndx=zeros(size(Un));
    dTndx(2:end,:)=N(2:end,:).*Un(1:end-1,:);
end
dTndx_p1=(n.^2).';
dTndx_m1=((-1).^(n-1).*n.^2).';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Synthesis of Chebyshev polynomials second derivative
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d2Tndx2=zeros(length(n),length(x));
if(not(isempty(ind)))
    d2Tndx2(2:end,ind)=(X(2:end,ind).*dTndx(2:end,ind)-N(2:end,ind).^2.*Tn(2:end,ind))./(1-X(2:end,ind).^2);
end
d2Tndx2_p1=(n.^2.*(-1+n.^2)/3).'; % x=+1 limit
d2Tndx2_m1=((-1).^(n).*(n.^2.*(-1+n.^2)/3)).'; % x=-1 limit

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setting limit values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(not(isempty(ind_p1)))
    Tn(:,ind_p1)=1;
    dTndx(:,ind_p1)=dTndx_p1;
    d2Tndx2(:,ind_p1)=d2Tndx2_p1;
end    
if(not(isempty(ind_m1)))
    Tn(:,ind_m1)=((-1).^n).';
    dTndx(:,ind_m1)=dTndx_m1;
    d2Tndx2(:,ind_m1)=d2Tndx2_m1;
end

return