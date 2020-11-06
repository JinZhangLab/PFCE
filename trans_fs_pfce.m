function [ b_s ] = trans_fs_pfce( Xm,Xs,ys,b_m)
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
%                     and placed at the element of the vector.
% Output:
%       b_s (S+1)¡Á1 The spectra of master used for calibration enhancment.
% References:
%    [1] A Parameter-Free Framework for Calibration Enhancement of Near-Infrared Spectroscopy Based on Correlation Constraint
%   Copyright Zhang Jin (zhangjin@mail.nankai.edu.cn).
% -------------------------------------------
options = optimoptions('fmincon','Display','none','Algorithm','sqp');
fun = @(b_s)(((ys-[ones(size(Xs,1),1) Xs]*b_s)'*(ys-[ones(size(Xs,1),1) Xs]*b_s))+...
    (([ones(size(Xm,1),1) Xm]*b_m-[ones(size(Xs,1),1) Xs]*b_s)'*([ones(size(Xm,1),1) Xm]*b_m-[ones(size(Xs,1),1) Xs]*b_s))); %Objective function

A = [];
b = [];
Aeq = [];
beq = [];
lb = [];
ub = [];
x0 = b_m;
b_s = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,@(b_s)mycon(b_s,b_m),options); %Optimization
end

function [c,ceq]=mycon(b_s,b_m)
r = corrcoef(b_s(2:end),b_m(2:end));
c = 0.98-r(1,2);
ceq = [];
end