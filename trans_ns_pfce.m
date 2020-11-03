function [ b_s ] = trans_ns_pfce( Xm,Xs,b_m)
% -------------------------------------------
% Calibration enhancement by using non-supervised Parameter-Free Framework
% for Calibration Enhancemen (ns-PEFC). 
% Input:
%       Xm (M¡ÁN) The spectra of master used for calibration enhancment. The
%                 M and N are the numer of samples and variables of master spectra,
%                 respectively.
%       Xs (M¡ÁS) The spectra of master used for calibration enhancment. The
%                 S is the numer of variables of slave spectra. The Xm and Xs
%                 should obtained from the same
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
options = optimoptions('fmincon','Display','iter','Algorithm','sqp','UseParallel',false);
fun = @(b_s)(([ones(size(Xm,1),1) Xm]*b_m-[ones(size(Xs,1),1) Xs]*b_s)'*([ones(size(Xm,1),1) Xm]*b_m-[ones(size(Xs,1),1) Xs]*b_s));%Objective function

A = [];
b = [];
Aeq = [];
beq = [];
lb = [];
ub = [];
x0 = b_m;
b_s = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,@(b_s)mycon(b_s,b_m),options);%Optimization
end


function [c,ceq]=mycon(b_s,b_m)
r = corr(b_s(2:end),b_m(2:end));
c = 0.98-r;
ceq = [];
end