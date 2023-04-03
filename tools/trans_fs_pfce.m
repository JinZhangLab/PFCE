function [ b_s ,cost] = trans_fs_pfce( Xm,Xs,ys,b_m,const,const_thres)
% -------------------------------------------
% Calibration enhancement by using full-supervised Parameter-Free Framework
% for Calibration Enhancemen (fs-PEFC).
% Input:
%       Xm (M¡ÁN) The spectra of master used for calibration enhancment. The
%                 M and N are the numer of samples and variables of master spectra,
%                 respectively.
%       Xs (M¡ÁN) The spectra of master used for calibration enhancment.
%       ys (M¡Á1) The reference value of master used for calibration enhancment.
%                 The Xm, Xs and ys should obtained from the same
%                 samples and correspond one by one.
%       b_m (N+1)¡Á1 The linear model coefficient of master. The model can be
%                     build by any linear fit method such as least
%                     regression, ridge regression and partial least squares
%                     regression, etc. The The intercept should be included
%                     and placed at the first element of the vector.
%       const    1,2,and 3 represent the Corr, L1, and L2 constraint, respectively
%       const_thres, contraint threshold, default 0.98. 
% Output:
%       b_s (N+1)¡Á1 The spectra of master used for calibration enhancment.
% References:
%    [1] A Parameter-Free Framework for Calibration Enhancement of Near-Infrared Spectroscopy Based on Correlation Constraint
%   Copyright Zhang Jin (zhangjin@mail.nankai.edu.cn).
% -------------------------------------------
if nargin <6
    const_thres = [0.98, 0.9, 0.01];
end
if nargin <5
    const = 1;
end
options = optimoptions('fmincon','Display','none','Algorithm','sqp','UseParallel',true);
%fun = @(b_s)(((ys-[ones(size(Xs,1),1) Xs]*b_s)'*(ys-[ones(size(Xs,1),1) Xs]*b_s))+...
%    (([ones(size(Xm,1),1) Xm]*b_m-[ones(size(Xs,1),1) Xs]*b_s)'*([ones(size(Xm,1),1) Xm]*b_m-[ones(size(Xs,1),1) Xs]*b_s))); %Objective function

A = [];
b = [];
Aeq = [];
beq = [];
lb = [];
ub = [];
ntask = 2;
if const == 0
    constfun = [];
elseif const ==1
    constfun = @(b_s)const1([b_s;b_m],ntask,const_thres);
elseif const ==2
    constfun = @(b_s)const2([b_s;b_m],ntask,const_thres);
elseif const ==3
    constfun = @(b_s)const3([b_s;b_m],ntask,const_thres);
else
    error("The selection of constraint is error!")
end

x0 = b_m;
[b_s,cost] = fmincon(@(b_s)CostFun_fs_pfce(b_s,b_m, Xm,Xs,ys),x0,A,b,Aeq,beq,lb,ub,constfun,options);%Optimization
end


