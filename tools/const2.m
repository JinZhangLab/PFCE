function [c,ceq]=const2(bglo,ntask,thes)
% Normalized L2 constraint (L2) used in PFCE
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
Vars = diag(B'*B);
VarErr = 0;
for i =1:ntask
    Err = B-B(:,i);
    Erri = Err;
    Erri(:,i) = [];
    Varsi = Vars;
    Varsi(i) = [];
    VarErr = VarErr + diag(Erri'*Erri)./sqrt(Vars(i)*Varsi);
end
dist = VarErr/(ntask*(ntask-1));
c = thes-(1-dist);
ceq = [];
end
