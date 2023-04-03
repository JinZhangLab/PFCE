function [optLV, RMSECV, Rcv, Ycv_hat, RMSEC, Rc, Y_hat] = plscv(X,y,maxLV, isPlot)
% Aim:
% Cross validation for partial least squares regression
% Input: 
% X,     matrix (n,p),  predictor matrix (assumed to be centered)
% Y,     matrix (n,m),  predictand (assumed to be centered)
% maxLV, scalar,        maximum number of PLS factors in cross validation
% isPlot,Logical value, Whether or not to draw the results 
% ------------------------------------------------------------------------
% Output: 
% optLV,  scaler,         optimal number of PLS factor from cross validtion
% -----------------------------------------------------------------------
% The following metrics were calculated based on prediction from cross validtion
% -----------------------------------------------------------------------
% RMSECV, vector (maxLV),   RMSE (root mean squared error) of cross validation
% Rcv,    vector (maxLV),   correlation coefficient (R) of cross validtion
% Ycv_hat,matrix (n,maxLV), prediction from cross validtion
% -----------------------------------------------------------------------
% The following metrics were calculated based on prediction form the PLS model with optLV
% -----------------------------------------------------------------------
% RMSEC,  vector (maxLV),   RMSE form the PLS model with optLV
% Rc,     vector (maxLV),   R between reference and prediction form the PLS model with optLV
% Y_hat,  matrix (n,maxLV), prediction form the PLS model with optLV
% Copyright J Zhang (zhangjin@mail.nankai.edu.cn).
[nrow, ncol]= size(X);
if nargin < 4
    isPlot = false;
end
if nargin <3
    maxLV = max([10,min([nrow-1, ncol-1, rank(X)])]);
end


Ycv_hat = zeros(nrow,maxLV);
for i =1:10% 10 fold cross validtion
    ival = mod(1:nrow,10) ==i-1;
    ical = ~ival;
    [B] = pls(X(ical,:),y(ical),maxLV);
    Ycv_hat(ival,:) = [ones(sum(ival),1), X(ival,:)] * B;
end
ErrorCV = Ycv_hat - y*ones(1,maxLV);
RMSECV = sqrt(sum((ErrorCV).^2)./(nrow));
Rcv = corr(Ycv_hat, y);
% stde = std(Eaug);
[~ ,optLV] = min(RMSECV);


[B] = pls(X,y,maxLV);
Y_hat = [ones(nrow,1), X] * B;
ErrorCV = Y_hat - y*ones(1,maxLV);
RMSEC = sqrt(sum((ErrorCV).^2)./(nrow));
Rc = corr(Y_hat, y);

if isPlot
    figure; 
    subplot(211);hold on;
    plot(1:maxLV, RMSEC, 'r-*');
    plot(1:maxLV, RMSECV,'g-^');
    legend(["RMSEC";"RMSECV";"optimal LV"])
    scatter([optLV, optLV], [RMSEC(optLV), RMSECV(optLV)], 50, 'ko');
    subplot(212); hold on;
    plot(1:maxLV, Rc, 'r-*');
    plot(1:maxLV, Rcv, 'g-^');
    scatter([optLV, optLV], [Rc(optLV), Rcv(optLV)], 50,'ko');
    legend(["RMSEC";"Rcv";"optimal LV"])
end
end

