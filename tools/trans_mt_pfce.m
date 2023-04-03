function [Ball,cost] = trans_mt_pfce( X,y,const,const_thres,optLV)
% -------------------------------------------
% Multitask learning by using Parameter-Free Framework for Calibration Enhancemen (mt-PEFC).
% Input:
% X, (ntask), cell array, The NIR spectra of n tasks, for all task, the row
% number (sample number) could be differ, but column numbers should be equal
% y, (ntask), cell array, The reference value of n tasks corressponding to the X.
% const,    1,2,and 3 represent the Corr, L1, and L2 constraint, respectively
% const_thres, contraint threshold, default 0.98. 
% Output:
% Ball,  (N+1)-by-ntask matrix, Model coefficients obtained for all n taks.
% cost, scalar, return of cost function with optimized models.
% References:
%  A Parameter-Free Framework for Calibration Enhancement of Near-Infrared Spectroscopy Based on Correlation Constraint
%  Copyright J Zhang (zhangjin@mail.nankai.edu.cn).
% -------------------------------------------
if length(X)~=length(y)
    error("Error, paired spectra and references value from multi source shold be provied!!!")
end
ntask = length(X);
rows = nan(ntask,1); cols = nan(ntask,1);
Xaug = [];yaug =[];
for i = 1:length(X)
    [rows(i),cols(i)] = size(X{i});
    Xaug = [Xaug;X{i}];
    yaug = [yaug;y{i}];
end


if nargin<5
    optLV = plscv(Xaug,yaug,20);
end
if nargin <4
    const_thres = 0.98;
end
if nargin <3
    const = 1;
end
options = optimoptions('fmincon','Display','none','Algorithm','sqp','UseParallel',true);

% Initize multitask learning parameters by using global model obtained by
% PLS with optimal lantent variables determined by 10 fold-cross validation
[B] = pls(Xaug,yaug,optLV);
bGlobal = B(:,end);
% To ensure convergence speed, regresssion coefficients by PLS was used.

b0 = repmat(bGlobal,ntask,1);
A = [];
b = [];
Aeq = [];
beq = [];
lb = [];
ub = [];

if const == 0
    constrfun = [];
elseif const ==1
    constrfun = @(bglo)const1(bglo,ntask,const_thres);
elseif const ==2
    constrfun = @(bglo)const2(bglo,ntask,const_thres);
elseif const ==3
    constrfun = @(bglo)const3(bglo,ntask,const_thres);
else
    error("The selection of constraint is error!")
end


[bGlobal,cost] = fmincon(@(bglo)CostFun_mt_pfce(bglo, X,y),b0,A,b,Aeq,beq,lb,ub,constrfun,options);%Optimization
Ball = reshape(bGlobal,[],ntask);
end



