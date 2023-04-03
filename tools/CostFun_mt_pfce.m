function [cost] = CostFun_mt_pfce(bglo, X, y)
% Cost function used in multi task PFCE
%   Copyright J Zhang (zhangjin@mail.nankai.edu.cn).
ntask = length(X);
B = reshape(bglo,[],ntask);
cost = 0;
for i = 1:ntask
    yerri = y{i}- [ones(size(X{i},1),1) X{i}]*B(:,i);
    cost = cost + yerri(:)'*yerri(:)/length(yerri);
end
end