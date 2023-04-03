function [cost] = CostFun_ss_pfce(b_s, Xs, ys)
% Cost function used in semi supervised PFCE
%   Copyright J Zhang (zhangjin@mail.nankai.edu.cn).
yerr = [ones(size(Xs,1),1) Xs]*b_s- ys;
cost = yerr(:)'*yerr(:);
end