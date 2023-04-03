function [c,ceq]=const1(bglo,ntask,thes)
% Correlation constraint (Corr) used in PFCE
%   [c,ceq]=const1(bglo,ntask,thes) returns the unequal and equal
%   constraint (c, ceq) using the linear regression coefficients (bglo) of n tasks
%   (ntask) with a specified threshould (thes)
%   bglo is vector with (1+n)*ntask elements. For all tasks, the corresponding 
%   coefficients were aranged according to their orders. For each task, the first and
%   rest figures of model coefficients are intercept and coefficients. The constraints
%   are only imposed on coefficients, but not on intercept.
%   Copyright J Zhang (zhangjin@mail.nankai.edu.cn).
B = reshape(bglo,[],ntask);
B = B(2:end,:);
R= corr(B);
R_triu = triu(R,1);
r = sum(R_triu(:))/sum(ntask-1:-1:1);
c = thes-r;
ceq = [];
end