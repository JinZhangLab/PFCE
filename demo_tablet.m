% -------------------------------------------
% Demo of Parameter-Free Framework for Calibration Enhancemen (PEFC)
% including non-supervised PEFC (ns-pefc), semi-supervised PEFC (ss-PEFC)
% and full-supervised PEFC (fs-PEFC).
% References:
%    [1] A Parameter-Free Framework for Calibration Enhancement of Near-Infrared Spectroscopy Based on Correlation Constraint
%   Copyright Zhang Jin (zhangjin@mail.nankai.edu.cn).
% -------------------------------------------
clear;clc; close all;
addpath(genpath(pwd))
load('Tablet.mat') % Load tablet dataset
% -------------------------------------------
%           Calibration
% -------------------------------------------
% Cross validation was performed, and the optimal number of component (n_comp) in
% PLS modeling was set to 3.
n_comp = 3;
B_m = pls(Xc,yc,3); % PLS modeling for master 
b_m = B_m(:,end); 
yhat_m  = [ones(size(Xm_t,1),1) Xm_t] * b_m; % Prediction of master
yhat_s  = [ones(size(Xs_t,1),1) Xs_t] * b_m; % Prediction of slave
[R_m, rmsep_m]   = rmse(yhat_m, yt); % R and RMSEP of master
[R_s, rmsep_s]   = rmse(yhat_s, yt); % R and RMSEP of slave
%-------------------------------------------
%           ns-pfce calibration maintenance
%-------------------------------------------
% Genearating the slave model by using ns-pfce
[ b_s_ns_pfce ] = trans_ns_pfce( Xm,Xs,b_m);
% Prediction of slave by using the model of ns-pfce
yhat_s_ns_pfce  = [ones(size(Xs_t,1),1) Xs_t]*b_s_ns_pfce; 
% R and RMSEP of slave aftet calibration maintenacne by ns-pfce
[R_s_ns_pfce, rmsep_s_ns_pfce]    = rmse(yhat_s_ns_pfce,yt);     
%-------------------------------------------
%           ss-pfce calibration enhancement
%-------------------------------------------
% Genearating the slave model by using ss-pfce
[ b_s_ss_pfce ]              = trans_ss_pfce( Xs,ys,b_m);
% Prediction of slave by using the model of ss-pfce
yhat_s_ss_pfce                = [ones(size(Xs_t,1),1) Xs_t]*b_s_ss_pfce;
% R and RMSEP of slave aftet calibration maintenacne by ns-pfce
[R_s_ss_pfce, rmsep_s_ss_pfce]    = rmse(yhat_s_ss_pfce,yt);
%-------------------------------------------
%           fs-pfce calibration enhancement
%-------------------------------------------
% Genearating the slave model by using fs-pfce
[ b_s_fs_pfce ]              = trans_fs_pfce( Xm,Xs,ys,b_m);
% Prediction of slave by using the model of ss-pfce
yhat_s_fs_pfce                = [ones(size(Xs_t,1),1) Xs_t]*b_s_fs_pfce;
% R and RMSEP of slave aftet calibration maintenacne by ns-pfce
[R_s_fs_pfce, rmsep_s_fs_pfce]    = rmse(yhat_s_fs_pfce,yt);
disp('------------Test end----------------')
