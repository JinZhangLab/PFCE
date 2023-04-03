function [c,ceq]=const3(bglo,ntask,thes)
% Normalized L1 constraint (L1) used in PFCE
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
norm1s = sum(abs(B));
norm1Err = 0;
for i = 1:ntask
    Err = B-B(:,i);
    Erri = Err;
    Erri(:,i) = [];
    norm1si = norm1s;
    norm1si(i) = [];
    norm1Err = norm1Err + sum(abs(Erri))./sqrt(norm1s(i)*norm1si);
end
dist = norm1Err/(ntask*(ntask-1));
c = thes-(1-dist);
ceq = [];
end
