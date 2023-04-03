% -------------------------------------------
% Demonstration of Parameter-Free Framework for Calibration Enhancemen
% (PEFC) with tablet dataset.
% References:
% [1]	Zhang J., Li B. Y., Hu Y., et al. A parameter-free framework for calibration enhancement of near-infrared spectroscopy based on correlation constraint [J]. Analytica Chimica Acta, 2021, 1142: 169-178.
%  Copyright J Zhang (zhangjin@mail.nankai.edu.cn).
% -------------------------------------------
% Clear workspace, command window, and close all figures
clear;clc; close all;

% Add the current directory to the MATLAB path
addpath(genpath(pwd))

% Load example tablet dataset
load('Tablet.mat') 

% The dataset contains:
% wv                 the wavelength of NIR spectra
% Xcal1, Xcal2       the NIR spectra in calibration set of master and slave,respectively
% Xtrans1, Xtrans2   the NIR spectra in standard set of the master and slave,respectively
% Xtest1 Xtest2      the NIR spectra in test set of the master and slave, respectively
% ycal ytrans ytest  the reference value of calibraiton, standard, test set, respectively 

%% -------------------------------------------
%           Calibration
% -------------------------------------------
% Number of PLS components
n_comp = 3;

% PLS modeling for master 
B_m = pls(Xcal1,ycal,n_comp); 
b_m = B_m(:,end); 

% Prediction of spectra in test set for master instrument with the master model (b_m)
yhat_m  = [ones(size(Xtest1,1),1) Xtest1] * b_m; 

% Prediction of spectra in test set for slave instrument directly with the master model (b_m)
yhat_s  = [ones(size(Xtest2,1),1) Xtest2] * b_m; 

% R and RMSEP of master model applied to the spectra of master and slave
[R_m, rmsep_m]   = rmse(yhat_m, ytest); 
[R_s, rmsep_s]   = rmse(yhat_s, ytest); 
%% ------------------------------------------
%         Global setting for PFCE 
%--------------------------------------------
% Three Constraints were introduced in PFCE2. 1, 2, and 3 for Corr, L2 and L1 constraints, respectively 
constraint_type = 1;

% Predefined constraint threshold. You can try more values for  improved results. 
constraint_thres = 0.98;
%% -------------------------------------------
%           ns-pfce 
%-------------------------------------------
% Genearating the slave model by using ns-pfce with Corr constraint
constraint_type = 1;
[ b_s_ns_pfce1 ] = trans_ns_pfce(Xtrans1,Xtrans2,b_m,constraint_type,constraint_thres);
% Prediction of slave by using the model enhanced by ns-pfce (b_s_ns_pfce1)
yhat_s_ns_pfce1  = [ones(size(Xtest2,1),1) Xtest2]*b_s_ns_pfce1; 
% R and RMSEP of slave obtained with enhanced model by ns-pfce (b_s_ns_pfce1)
[R_s_ns_pfce1, rmsep_s_ns_pfce1]    = rmse(yhat_s_ns_pfce1,ytest); 


% Repeat with L2 constraint
constraint_type = 2;
[ b_s_ns_pfce2 ] = trans_ns_pfce(Xtrans1,Xtrans2,b_m,constraint_type,constraint_thres);
yhat_s_ns_pfce2  = [ones(size(Xtest2,1),1) Xtest2]*b_s_ns_pfce2; 
[R_s_ns_pfce2, rmsep_s_ns_pfce2]    = rmse(yhat_s_ns_pfce2,ytest); 


% Repeat with L1 constraint
constraint_type = 3;
[ b_s_ns_pfce3 ] = trans_ns_pfce(Xtrans1,Xtrans2,b_m,constraint_type,constraint_thres);
yhat_s_ns_pfce3  = [ones(size(Xtest2,1),1) Xtest2]*b_s_ns_pfce3; 
[R_s_ns_pfce3, rmsep_s_ns_pfce3]    = rmse(yhat_s_ns_pfce3,ytest); 
%% -------------------------------------------
%           ss-pfce 
%-------------------------------------------
% Genearating the slave model by using ss-pfce
constraint_type = 1;
[ b_s_ss_pfce1 ]              = trans_ss_pfce(Xtrans2,ytrans,b_m, constraint_type,constraint_thres);
% Prediction of slave by using the model of ss-pfce
yhat_s_ss_pfce1                = [ones(size(Xtest2,1),1) Xtest2]*b_s_ss_pfce1;
% R and RMSEP of slave aftet calibration maintenacne by ss-pfce
[R_s_ss_pfce1, rmsep_s_ss_pfce1]    = rmse(yhat_s_ss_pfce1,ytest);


% Repeat with L2 constraint
constraint_type = 2;
[ b_s_ss_pfce2 ] = trans_ss_pfce( Xtrans2,ytrans,b_m,constraint_type,constraint_thres);
yhat_s_ss_pfce2  = [ones(size(Xtest2,1),1) Xtest2]*b_s_ss_pfce2; 
[R_s_ss_pfce2, rmsep_s_ss_pfce2]    = rmse(yhat_s_ss_pfce2,ytest); 


% Repeat with L1 constraint
constraint_type = 3;
[ b_s_ss_pfce3 ] = trans_ss_pfce( Xtrans2,ytrans,b_m,constraint_type,constraint_thres);
yhat_s_ss_pfce3  = [ones(size(Xtest2,1),1) Xtest2]*b_s_ss_pfce3; 
[R_s_ss_pfce3, rmsep_s_ss_pfce3]    = rmse(yhat_s_ss_pfce3,ytest); 
%% -------------------------------------------
%           fs-pfce 
%-------------------------------------------
% Genearating the slave model by using fs-pfce
constraint_type = 1;
[ b_s_fs_pfce1 ]              = trans_fs_pfce(Xtrans1,Xtrans2,ytrans,b_m,constraint_type,constraint_thres);
% Prediction of slave by using the model of fs-pfce
yhat_s_fs_pfce1                = [ones(size(Xtest2,1),1) Xtest2]*b_s_fs_pfce1;
% R and RMSEP of slave aftet calibration maintenacne by fs-pfce
[R_s_fs_pfce1, rmsep_s_fs_pfce1]    = rmse(yhat_s_fs_pfce1,ytest);


% Repeat with L2 constraint
constraint_type= 2;
[ b_s_fs_pfce2 ] = trans_fs_pfce( Xtrans1,Xtrans2,ytrans,b_m,constraint_type,constraint_thres);
yhat_s_fs_pfce2  = [ones(size(Xtest2,1),1) Xtest2]*b_s_fs_pfce2; 
[R_s_fs_pfce2, rmsep_s_fs_pfce2]    = rmse(yhat_s_fs_pfce2,ytest); 


% Repeat with L1 constraint
constraint_type = 3;
[ b_s_fs_pfce3 ] = trans_fs_pfce( Xtrans1,Xtrans2,ytrans,b_m,constraint_type,constraint_thres);
yhat_s_fs_pfce3  = [ones(size(Xtest2,1),1) Xtest2]*b_s_fs_pfce3; 
[R_s_fs_pfce3, rmsep_s_fs_pfce3]    = rmse(yhat_s_fs_pfce3,ytest);
%% -------------------------------------------
%          mt-PFCE
%-------------------------------------------
% Enhancing the master and slave model by using mt-pfce with Corr constraint
constraint_type = 1;
[B_mt_pfce1]              = trans_mt_pfce({Xcal1,Xtrans2},{ycal,ytrans},constraint_type,constraint_thres,n_comp);
% Prediction of slave by using the model of fs-pfce
yhat_m_mt_pfce1                = [ones(size(Xtest1,1),1) Xtest1]*B_mt_pfce1(:,1);
yhat_s_mt_pfce1                = [ones(size(Xtest2,1),1) Xtest2]*B_mt_pfce1(:,2);
% R and RMSEP of slave aftet calibration maintenacne by fs-pfce
[R_m_mt_pfce1, rmsep_m_mt_pfce1]    = rmse(yhat_m_mt_pfce1,ytest);
[R_s_mt_pfce1, rmsep_s_mt_pfce1]    = rmse(yhat_s_mt_pfce1,ytest);


% Repeat with L2 constraint
constraint_type = 2;
[ B_mt_pfce2]              = trans_mt_pfce({Xcal1,Xtrans2},{ycal,ytrans},constraint_type,constraint_thres,n_comp);
yhat_m_mt_pfce2                = [ones(size(Xtest1,1),1) Xtest1]*B_mt_pfce2(:,1);
yhat_s_mt_pfce2                = [ones(size(Xtest2,1),1) Xtest2]*B_mt_pfce2(:,2);
[R_m_mt_pfce2, rmsep_m_mt_pfce2]    = rmse(yhat_m_mt_pfce2,ytest);
[R_s_mt_pfce2, rmsep_s_mt_pfce2]    = rmse(yhat_s_mt_pfce2,ytest);


% Repeat with L1 constraint
constraint_type = 3;
[ B_mt_pfce3]              = trans_mt_pfce({Xcal1,Xtrans2},{ycal,ytrans},constraint_type,constraint_thres,n_comp);
yhat_m_mt_pfce3                = [ones(size(Xtest1,1),1) Xtest1]*B_mt_pfce3(:,1);
yhat_s_mt_pfce3                = [ones(size(Xtest2,1),1) Xtest2]*B_mt_pfce3(:,2);
[R_m_mt_pfce3, rmsep_m_mt_pfce3]    = rmse(yhat_m_mt_pfce3,ytest);
[R_s_mt_pfce3, rmsep_s_mt_pfce3]    = rmse(yhat_s_mt_pfce3,ytest);
disp('------------Test end----------------')
%%
% plot the NS-PFCE results
figure;

subplot(2,3,1); plot(ytest,yhat_m,'o'); 
hold on; plot([min(ycal) max(ycal)],[min(ycal) max(ycal)],'-');
axis equal; xlabel('Reference Value'); ylabel('Prediction'); 
title(sprintf('%s \n R = %.2f, RMSEP = %.2f', "Master", R_m, rmsep_m));

subplot(2,3,2); plot(ytest,yhat_s,'o');
hold on; plot([min(ycal) max(ycal)],[min(ycal) max(ycal)],'-');
axis equal; xlabel('Reference Value'); ylabel('Prediction'); 
title(sprintf('%s \n R = %.2f, RMSEP = %.2f',"Slave", R_s, rmsep_s));

subplot(2,3,4); plot(ytest,yhat_s_ns_pfce1,'o');
hold on; plot([min(ycal) max(ycal)],[min(ycal) max(ycal)],'-');
axis equal; xlabel('Reference Value'); ylabel('Prediction');
title(sprintf('%s \n R = %.2f, RMSEP = %.2f',"NS-PFCE (Corr)", R_s_ns_pfce1, rmsep_s_ns_pfce1));

subplot(2,3,5); plot(ytest,yhat_s_ns_pfce2,'o'); 
hold on; plot([min(ycal) max(ycal)],[min(ycal) max(ycal)],'-');
axis equal; xlabel('Reference Value'); ylabel('Prediction');
title(sprintf('%s \n R = %.2f, RMSEP = %.2f', "NS-PFCE (L2)", R_s_ns_pfce2, rmsep_s_ns_pfce2));

subplot(2,3,6); plot(ytest,yhat_s_ns_pfce3,'o');
hold on; plot([min(ycal) max(ycal)],[min(ycal) max(ycal)],'-');
axis equal; xlabel('Reference Value'); ylabel('Prediction');
title(sprintf('%s \n R = %.2f, RMSEP = %.2f', "NS-PFCE (L1)", R_s_ns_pfce3, rmsep_s_ns_pfce3));

sgtitle('Calibration enhancement with NS-PFCE for Tablet dataset');


% plot the SS-PFCE results
figure;

subplot(2,3,1); plot(ytest,yhat_m,'o'); 
hold on; plot([min(ycal) max(ycal)],[min(ycal) max(ycal)],'-');
axis equal; xlabel('Reference Value'); ylabel('Prediction'); 
title(sprintf('%s \n R = %.2f, RMSEP = %.2f', "Master", R_m, rmsep_m));

subplot(2,3,2); plot(ytest,yhat_s,'o');
hold on; plot([min(ycal) max(ycal)],[min(ycal) max(ycal)],'-');
axis equal; xlabel('Reference Value'); ylabel('Prediction'); 
title(sprintf('%s \n R = %.2f, RMSEP = %.2f',"Slave", R_s, rmsep_s));

subplot(2,3,4); plot(ytest,yhat_s_ss_pfce1,'o');
hold on; plot([min(ycal) max(ycal)],[min(ycal) max(ycal)],'-');
axis equal; xlabel('Reference Value'); ylabel('Prediction');
title(sprintf('%s \n R = %.2f, RMSEP = %.2f',"SS-PFCE (Corr)", R_s_ss_pfce1, rmsep_s_ss_pfce1));

subplot(2,3,5); plot(ytest,yhat_s_ss_pfce2,'o'); 
hold on; plot([min(ycal) max(ycal)],[min(ycal) max(ycal)],'-');
axis equal; xlabel('Reference Value'); ylabel('Prediction');
title(sprintf('%s \n R = %.2f, RMSEP = %.2f', "SS-PFCE (L2)", R_s_ss_pfce2, rmsep_s_ss_pfce2));

subplot(2,3,6); plot(ytest,yhat_s_ss_pfce3,'o');
hold on; plot([min(ycal) max(ycal)],[min(ycal) max(ycal)],'-');
axis equal; xlabel('Reference Value'); ylabel('Prediction');
title(sprintf('%s \n R = %.2f, RMSEP = %.2f', "SS-PFCE (L1)", R_s_ss_pfce3, rmsep_s_ss_pfce3));

sgtitle('Calibration enhancement with SS-PFCE for Tablet dataset');


% plot the FS-PFCE results
figure;

subplot(2,3,1); plot(ytest,yhat_m,'o'); 
hold on; plot([min(ycal) max(ycal)],[min(ycal) max(ycal)],'-');
axis equal; xlabel('Reference Value'); ylabel('Prediction'); 
title(sprintf('%s \n R = %.2f, RMSEP = %.2f', "Master", R_m, rmsep_m));

subplot(2,3,2); plot(ytest,yhat_s,'o');
hold on; plot([min(ycal) max(ycal)],[min(ycal) max(ycal)],'-');
axis equal; xlabel('Reference Value'); ylabel('Prediction'); 
title(sprintf('%s \n R = %.2f, RMSEP = %.2f',"Slave", R_s, rmsep_s));

subplot(2,3,4); plot(ytest,yhat_s_fs_pfce1,'o');
hold on; plot([min(ycal) max(ycal)],[min(ycal) max(ycal)],'-');
axis equal; xlabel('Reference Value'); ylabel('Prediction');
title(sprintf('%s \n R = %.2f, RMSEP = %.2f',"FS-PFCE (Corr)", R_s_fs_pfce1, rmsep_s_fs_pfce1));

subplot(2,3,5); plot(ytest,yhat_s_fs_pfce2,'o'); 
hold on; plot([min(ycal) max(ycal)],[min(ycal) max(ycal)],'-');
axis equal; xlabel('Reference Value'); ylabel('Prediction');
title(sprintf('%s \n R = %.2f, RMSEP = %.2f', "FS-PFCE (L2)", R_s_fs_pfce2, rmsep_s_fs_pfce2));

subplot(2,3,6); plot(ytest,yhat_s_fs_pfce3,'o');
hold on; plot([min(ycal) max(ycal)],[min(ycal) max(ycal)],'-');
axis equal; xlabel('Reference Value'); ylabel('Prediction');
title(sprintf('%s \n R = %.2f, RMSEP = %.2f', "FS-PFCE (L1)", R_s_fs_pfce3, rmsep_s_fs_pfce3));

sgtitle('Calibration enhancement with FS-PFCE for Tablet dataset');


% plot the MT-PFCE results
figure;

subplot(3,3,1); plot(ytest,yhat_m,'o'); 
hold on; plot([min(ycal) max(ycal)],[min(ycal) max(ycal)],'-');
axis equal; xlabel('Reference Value'); ylabel('Prediction'); 
title(sprintf('%s \n R = %.2f, RMSEP = %.2f', "Master", R_m, rmsep_m));

subplot(3,3,2); plot(ytest,yhat_s,'o');
hold on; plot([min(ycal) max(ycal)],[min(ycal) max(ycal)],'-');
axis equal; xlabel('Reference Value'); ylabel('Prediction'); 
title(sprintf('%s \n R = %.2f, RMSEP = %.2f',"Slave", R_s, rmsep_s));

subplot(3,3,4); plot(ytest,yhat_m_mt_pfce1,'o');
hold on; plot([min(ycal) max(ycal)],[min(ycal) max(ycal)],'-');
axis equal; xlabel('Reference Value'); ylabel('Prediction');
title(sprintf('%s \n R = %.2f, RMSEP = %.2f',"MT-PFCE (Corr), Master", R_m_mt_pfce1, rmsep_m_mt_pfce1));

subplot(3,3,5); plot(ytest,yhat_m_mt_pfce2,'o'); 
hold on; plot([min(ycal) max(ycal)],[min(ycal) max(ycal)],'-');
axis equal; xlabel('Reference Value'); ylabel('Prediction');
title(sprintf('%s \n R = %.2f, RMSEP = %.2f', "MT-PFCE (L2), Master", R_m_mt_pfce2, rmsep_m_mt_pfce2));

subplot(3,3,6); plot(ytest,yhat_m_mt_pfce3,'o');
hold on; plot([min(ycal) max(ycal)],[min(ycal) max(ycal)],'-');
axis equal; xlabel('Reference Value'); ylabel('Prediction');
title(sprintf('%s \n R = %.2f, RMSEP = %.2f', "MT-PFCE (L1), Master", R_m_mt_pfce3, rmsep_m_mt_pfce3));

subplot(3,3,7); plot(ytest,yhat_s_mt_pfce1,'o');
hold on; plot([min(ycal) max(ycal)],[min(ycal) max(ycal)],'-');
axis equal; xlabel('Reference Value'); ylabel('Prediction');
title(sprintf('%s \n R = %.2f, RMSEP = %.2f',"MT-PFCE (Corr), Slave", R_s_mt_pfce1, rmsep_s_mt_pfce1));

subplot(3,3,8); plot(ytest,yhat_s_mt_pfce2,'o'); 
hold on; plot([min(ycal) max(ycal)],[min(ycal) max(ycal)],'-');
axis equal; xlabel('Reference Value'); ylabel('Prediction');
title(sprintf('%s \n R = %.2f, RMSEP = %.2f', "MT-PFCE (L2), Slave", R_s_mt_pfce2, rmsep_s_mt_pfce2));

subplot(3,3,9); plot(ytest,yhat_s_mt_pfce3,'o');
hold on; plot([min(ycal) max(ycal)],[min(ycal) max(ycal)],'-');
axis equal; xlabel('Reference Value'); ylabel('Prediction');
title(sprintf('%s \n R = %.2f, RMSEP = %.2f', "MT-PFCE (L1), Slave", R_s_mt_pfce3, rmsep_s_mt_pfce3));

sgtitle('Calibration enhancement with FS-PFCE for Tablet dataset');
%%
Results0 = [rmsep_m, R_m; 
    rmsep_s, R_s];

Results1 = [rmsep_s_ns_pfce1, R_s_ns_pfce1,  rmsep_s_ns_pfce2, R_s_ns_pfce2,  rmsep_s_ns_pfce3, R_s_ns_pfce3;
            rmsep_s_ss_pfce1, R_s_ss_pfce1,  rmsep_s_ss_pfce2, R_s_ss_pfce2,  rmsep_s_ss_pfce3, R_s_ss_pfce3;
            rmsep_s_fs_pfce1, R_s_fs_pfce1,  rmsep_s_fs_pfce2, R_s_fs_pfce2,  rmsep_s_fs_pfce3, R_s_fs_pfce3;
            rmsep_m_mt_pfce1, R_m_mt_pfce1,  rmsep_m_mt_pfce2, R_m_mt_pfce2,  rmsep_m_mt_pfce3, R_m_mt_pfce3;
            rmsep_s_mt_pfce1, R_s_mt_pfce1,  rmsep_s_mt_pfce2, R_s_mt_pfce2,  rmsep_s_mt_pfce3, R_s_mt_pfce3;];
