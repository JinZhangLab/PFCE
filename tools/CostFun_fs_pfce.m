function [cost] = CostFun_fs_pfce(b_s,b_m, Xm,Xs,ys)
% Cost function used in full supervised PFCE
%   Copyright J Zhang (zhangjin@mail.nankai.edu.cn).
yerr1 = [ones(size(Xm,1),1) Xm]*b_m-[ones(size(Xs,1),1) Xs]*b_s;
yerr2 = ys - [ones(size(Xs,1),1) Xs]*b_s;
cost = yerr1(:)'*yerr1(:)/length(yerr1)+yerr2(:)'*yerr2(:)/length(yerr2);
end
