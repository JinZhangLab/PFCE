function [ b_s ] = trans_ss_pfce( Xs,ys,b_mm)
% -------------------------------------------
% Calibration enhancement by using semi-supervised Parameter-Free Framework
% for Calibration Enhancemen (ss-PEFC). 
% Input:
%       Xm (M×N) The spectra of master used for calibration enhancment. The
%                 M and N are the numer of samples and variables of master spectra,
%                 respectively.
%       ys (M×1) The reference value of master used for calibration enhancment.
%                 The Xm, Xs and ys should obtained from the same
%                 samples and correspond one by one.
%       b_mm (N+1)×1 The linear model coefficient of master. The model can be 
%                     build by any linear fit method such as least
%                     regression, ridge regression and partial least squares
%                     regression, etc. The The intercept should be included 
%                     and placed at the element of the vector.
% Output:
%       b_s (S+1)×1 The spectra of master used for calibration enhancment.
% References:
%    [1] A Parameter-Free Framework for Calibration Enhancement of Near-Infrared Spectroscopy Based on Correlation Constraint
%   Copyright Zhang Jin (zhangjin@mail.nankai.edu.cn).
% -------------------------------------------
global b_m
b_m = b_mm;
options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
fun = @(b_s)((ys-[ones(size(Xs,1),1) Xs]*b_s)'*(ys-[ones(size(Xs,1),1) Xs]*b_s));%目标函数

A = [];
b = [];
Aeq = [];
beq = [];
lb = [];
ub = [];
x0 = b_m;
b_s = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,@mycon,options);%线性
end


function [c,ceq]=mycon(b_s)
global b_m
r = corrcoef(b_s(2:end),b_m(2:end));
c = 0.98-r(1,2);
ceq = [];
end