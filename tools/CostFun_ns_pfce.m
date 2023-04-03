function [cost] = CostFun_ns_pfce(b_s,b_m, Xm,Xs)
% Cost function used in null supervised PFCE
%   Copyright J Zhang (zhangjin@mail.nankai.edu.cn).
yerr = [ones(size(Xm,1),1) Xm]*b_m-[ones(size(Xs,1),1) Xs]*b_s;
cost = yerr(:)'*yerr(:)/length(yerr);
end