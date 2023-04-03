function [B,C,P,T,U,R,R2X,R2Y] = pls(X,Y,A)
% Aim:
% centralized partial least squares regression based on simpls
% Input: 
% X, matrix (n,p), predictor matrix (assumed to be centered)
% Y, matrix (n,m), predictand (assumed to be centered)
% A, scalar, number of PLS factors
% ------------------------------------------------------------------------
% Output: 
% B, matrix (p+1,m), regression coefficients with intercept in the 1st row
% C, matrix (m,h), Y loadings
% P, matrix (p,h), X loadings
% T, matrix (n,h), X scores (standardized) 
% U, matrix (n,h), Y scores
% R, matrix (p,h), X weights
% R2X, vecor (1,h), X-variance
% R2Y, vecor (1,h), Y-variance
% Copyright Zhang Jin (zhangjin@mail.nankai.edu.cn).
mx = mean(X);
my = mean(Y);
Xnew = X - mx;
Ynew = Y - my;
[B,C,P,T,U,R,R2X,R2Y] = plssim(Xnew,Ynew,A);
B = [my-mx*B; B];
end

function [B,C,P,T,U,R,R2X,R2Y]=plssim(X,Y,A,S,XtX)
% ------------------------------------------------------------------------
% Function: [B,C,P,T,U,R,R2X,R2Y]=plssim(X,Y,h,S,XtX)
% ------------------------------------------------------------------------
% Aim:
% Partial Least Squares for tall X matrices, SIM-PLS 
% ------------------------------------------------------------------------
% Input: 
% X, matrix (n,p), predictor matrix (assumed to be centered)
% Y, matrix (n,m), predictand (assumed to be centered)
% h, scalar, number of PLS factors
% S, matrix (n,m), S=X'*Y
% XtX, matrix (n,n), XtX=X'*X (boosts speed for tall X matrices n>>p)
% ------------------------------------------------------------------------
% Output: 
% B, matrix (p,m), regression coefficients
% C, matrix (m,h), Y loadings
% P, matrix (p,h), X loadings
% T, matrix (n,h), X scores (standardized) 
% U, matrix (n,h), Y scores
% R, matrix (p,h), X weights
% R2X, vecor (1,h), X-variance
% R2Y, vecor (1,h), Y-variance
% ------------------------------------------------------------------------
% Example: 
% 1/ for non tall matrix: 
% [B]=plssim(X,Y,10,[],X'*Y)
% 2/ for tall matrix:     
% [B]=plssim(X,Y,10,X'*Y,X'*X)
% ------------------------------------------------------------------------
% The above routine is included into the toolbox with personal agreement 
% of its author Sijmen de Jong
% ------------------------------------------------------------------------
% Reference:
% S. de Jong, SIMPLS: An alternative approach to partial least squares 
% regression, Chemometrics and Intelligent Laboratory Systems, 
% 18 (1993) 251-263
[n,px]=size(X); 
[n,m]=size(Y);   				

if nargin<5, 
    S=[]; 
end

if isempty(S) 
    S=(Y'*X)'; 
end		            % if XtX not inputted, S=[]; always when S=[] then S=(Y'*X)'

if nargin<4 
    XtX=[]; 
end					% if S is not inputted, XtX=[];

if isempty(XtX) & n>3*px 
    XtX=X'*X; 
end			        % when XtX=[] and X is very "tall", the booster XtX is calculated

if nargin<3 
    A=10; 
end

A=min([A px n-1]);			% if A is not inputted, then the defaul A is min[10 px n-1]
T=zeros(n,A); 
U=T;						% initialization of variables
R=zeros(px,A); 
P=R; 
V=R;
C=zeros(m,A); 
R2Y=zeros(1,A);
z=zeros(m,1); 
v=zeros(px,1);

if n>px 
    S0 = S; 
end
StS=S'*S;				    % SIMPLS algorithm
nm1=n-1;
tol=0;

for a=1:A
    StS=StS-z*z'; 
    [Q,LAMBDA]=eig(StS); 
    [lambda,j]=max(diag(LAMBDA)); 
    q=Q(:,j(1));
    r=S*q;
    t=X*r;
    
    if isempty(XtX)
        p=(t'*X)'; 
    else
        p=XtX*r;
    end
    
    if n>px, 
        d=sqrt(r'*p/nm1); 
    else 
        d=sqrt(t'*t/nm1); 
    end
    
    v=p-V(:,1:max(1,a-1))*(p'*V(:,1:max(1,a-1)))'; 
    v=v/sqrt(v'*v); 
    z=(v'*S)'; 
    S=S-v*z'; 
    V(:,a)=v;
    R(:,a)=r/d; 						    % X weights
    P(:,a)=p/(d*nm1); 						% X loadings
    T(:,a)=t/d;							    % X scores
    U(:,a)=Y*q;							    % Y scores
    C(:,a)=q*(lambda(1)/(nm1*d)); 			% Y loadings
    R2Y(1,a)=lambda(1)/d;					% Y-variance accounted for
    B(:,a)=R*C';					        % B-coefficients of the regression Y on X
   
end

clear StS V LAMBDA Q p q r t v z;

if d<tol,
    A=a-1; 
    a=A; 
    T=T(:,1:A); 
    U=U(:,1:A); 
    R=R(:,1:A); 
    P=P(:,1:A); 
    C=C(:,1:A);
end

while a>1
    U(:,a) = U(:,a)-T(:,1:a-1)*(U(:,a)'*T(:,1:a-1)/nm1)'; 
    a=a-1; 
end

if isempty(XtX),
    sumX2=sum(X.^2);
else 
    sumX2=sum(diag(XtX)); 
end

R2X=100*nm1/sum(sumX2)*(sum(P.^2)); 
R2Y=100/nm1/sum(sum(Y.^2))*(R2Y(1:A).^2);
end