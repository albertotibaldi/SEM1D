function[xvet,wvet]=quadad(a,N)

m=0:1:N-2;
alpha=(m+1)./(2*m+1).*sqrt((2*m+1)./(2*m+3));
m=0:1:N-1;
beta=0*m;

T=diag(alpha,1)+diag(beta,0)+diag(alpha,-1);
[U,D]=eig(T);
xvet=diag(D).';
wvet=U(1,:).^2;
wvet=wvet*2;

xvet=xvet*a;
wvet=wvet*a;