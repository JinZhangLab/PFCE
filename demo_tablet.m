% -------------------------------------------
% Demo of Parameter-Free Framework for Calibration Enhancemen (PEFC)
% including non-supervised PEFC (ns-pefc), semi-supervised PEFC (ss-PEFC)
% and full-supervised PEFC (fs-PEFC).
% References:
% [1] J Zhang, BY Li, Y Hu, et. al. A Parameter-Free Framework for Calibration Enhancement of Near-Infrared Spectroscopy Based on Correlation Constraint, analytica chimica acta, 2021, 1142: 169-178. https://doi.org/10.1016/j.aca.2020.11.006
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
% Calibration set includes Xc and yc, which are spectra measured on master instrument and the correspondig reference values, respectively.
B_m = pls(Xc,yc,3); % PLS modeling for master model 
% Master model (b_m) are given in the form of (N+1)*1 array, and the first one element is intercept, and the last N elements are coefficients. 
b_m = B_m(:,end); 

% Test set includes Xm_t and Xs_t, which are spectra measured on master and slave instrument, respectively.
yhat_m  = [ones(size(Xm_t,1),1) Xm_t] * b_m; % Prediction of test spectra of master with master model
yhat_s  = [ones(size(Xs_t,1),1) Xs_t] * b_m; % Prediction of test spectra of slave with master model
[R_m, rmsep_m]   = rmse(yhat_m, yt); % R and RMSEP of spectral prediction of master with master model
[R_s, rmsep_s]   = rmse(yhat_s, yt); % R and RMSEP of spectral prediction of slave with master model
%-------------------------------------------
%           ns-pfce calibration maintenance
%-------------------------------------------
% Genearating the enhanced model for slave instrument by using ns-pfce
% Xm and Xs are paired spectra measured on master and slave instruments from the same batch of standard samples.
[ b_s_ns_pfce ] = trans_ns_pfce( Xm,Xs,b_m);
% Prediction of slave by using the model of ns-pfce
yhat_s_ns_pfce  = [ones(size(Xs_t,1),1) Xs_t]*b_s_ns_pfce; % Prediction of test spectra of slave with the enhanced model by ns-pfce
% R and RMSEP of slave aftet calibration maintenacne by ns-pfce
[R_s_ns_pfce, rmsep_s_ns_pfce]    = rmse(yhat_s_ns_pfce,yt); % R and RMSEP of spectral prediction of slave with enhanced model by ns-pfce
%-------------------------------------------
%           ss-pfce calibration enhancement
%-------------------------------------------
% Genearating the enhanced model for slave instrument by using ss-pfce
% Xs are spectra measured only on slave instruments from a few samples, and ys are the corresponding reference values.
[ b_s_ss_pfce ]              = trans_ss_pfce( Xs,ys,b_m);
% Prediction of slave by using the model of ss-pfce
yhat_s_ss_pfce                = [ones(size(Xs_t,1),1) Xs_t]*b_s_ss_pfce; % Prediction of test spectra of slave with the enhanced model by ss-pfce
% R and RMSEP of slave aftet calibration maintenacne by ns-pfce
[R_s_ss_pfce, rmsep_s_ss_pfce]    = rmse(yhat_s_ss_pfce,yt);% R and RMSEP of spectral prediction of slave with enhanced model by ss-pfce
%-------------------------------------------
%           fs-pfce calibration enhancement
%-------------------------------------------
% Genearating the enhanced model for slave instrument by using fs-pfce
% Xm and Xs are paired spectra measured on master and slave instruments from the same batch of standard samples,
% and ys are the corresponding reference values.
[ b_s_fs_pfce ]              = trans_fs_pfce( Xm,Xs,ys,b_m);
% Prediction of slave by using the model of ss-pfce
yhat_s_fs_pfce                = [ones(size(Xs_t,1),1) Xs_t]*b_s_fs_pfce; % Prediction of test spectra of slave with the enhanced model by fs-pfce
% R and RMSEP of slave aftet calibration maintenacne by ns-pfce
[R_s_fs_pfce, rmsep_s_fs_pfce]    = rmse(yhat_s_fs_pfce,yt); % R and RMSEP of test spectral prediction of slave with enhanced model by fs-pfce
disp('------------Test end----------------')
