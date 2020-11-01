function [R,rms]=rmse(y,yhat)
%rmse Root-mean-square error.
%   [R,rms] = rmse(y,yhat) computes the correlation coefficient (R) and
%   root mean square error (rms) betwwen the reference value (y) and predicted
%   value (yhat). R and rms are both scales.
%   Copyright Zhang Jin (zhangjin@mail.nankai.edu.cn).
R=corr(y(:),yhat(:));
rms=sqrt(sum((y(:)-yhat(:)).^2)/length(y));
end