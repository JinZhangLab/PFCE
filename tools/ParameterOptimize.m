function [optThres,RMSECal,RMSECVal] = ParameterOptimize(method,inputData,thresRange,ConstType,nFold)
% parameter optimization for PFCE to find out optimal contraint threshould with nfold
% cross validtion.
% method      The method for which the parameter are to be optimized.
%             (NS_PFCE,SS_PFCE,FS_PFCE,MT_PFCE)
% inpuptData  The input data for the correspoindind method, shold be cell
%             format.
% ThresRnage  The constraint range to be searched, default 0.00:0.02:1.
% ConstType   Constraint type.
% nFold       The number of fold in cross validtion.
% copyright J Zhang (zhangjin@mail.nankai.edu.cn)

if nargin < 5
    nFold = min([10, size(inputData{1})]);
end
if nargin <4
    ConstType = 1;
end
if nargin <3
    thresRange = 0.00:0.02:1;
end
if method == "NS_PFCE" || method == "SS_PFCE" || method == "FS_PFCE"
    RMSECal = zeros(nFold,length(thresRange));
    RMSECVal = zeros(nFold,length(thresRange));
elseif method == "MT_PFCE"
    RMSECal = zeros(nFold,length(thresRange),length(inputData{1}));
    RMSECVal = zeros(nFold,length(thresRange),length(inputData{1}));
else
    error("The inputted method in not supportted.")
end
for i =1:nFold% N fold cross validtion
    for j = 1:length(thresRange)
        thresi = thresRange(j);
        if method =="NS_PFCE"
            Xtrans1 = inputData{1};
            Xtrans2 = inputData{2};
            ytrans = inputData{3};
            b_m = inputData{4};
            nsample = size(Xtrans1,1);
            ival = mod(1:nsample,nFold) ==i-1;
            ical = ~ival;
            [b_s, Cost_cali ] = trans_ns_pfce(Xtrans1(ical,:),Xtrans2(ical,:),b_m,ConstType,thresi);
            yhat_c_si  = [ones(sum(ical),1) Xtrans2(ical,:)] * b_s;
            [Rc_i, rmsec_i]   = rmse(yhat_c_si, ytrans(ical));
            yhat_si  = [ones(sum(ival),1) Xtrans2(ival,:)] * b_s;
            [R_i, rmsep_i]   = rmse(yhat_si, ytrans(ival));
        elseif method =="SS_PFCE"
            Xtrans2 = inputData{1};
            ytrans = inputData{2};
            b_m = inputData{3};
            nsample = length(ytrans);
            ival = mod(1:nsample,nFold) ==i-1;
            ical = ~ival;
            [b_s, Cost_cali] = trans_ss_pfce(Xtrans2(ical,:),ytrans(ical,:),b_m, ConstType,thresi);
            yhat_c_si  = [ones(sum(ical),1) Xtrans2(ical,:)] * b_s;
            [Rc_i, rmsec_i]   = rmse(yhat_c_si, ytrans(ical));
            yhat_si  = [ones(sum(ival),1) Xtrans2(ival,:)] * b_s;
            [R_i, rmsep_i]   = rmse(yhat_si, ytrans(ival));
        elseif method =="FS_PFCE"
            Xtrans1 = inputData{1};
            Xtrans2 = inputData{2};
            ytrans = inputData{3};
            b_m = inputData{4};
            nsample = length(ytrans);
            ival = mod(1:nsample,nFold) ==i-1;
            ical = ~ival;
            [b_s, Cost_cali] = trans_fs_pfce(Xtrans1(ical,:),Xtrans2(ical,:),ytrans(ical),b_m,ConstType,thresi);
            yhat_c_si  = [ones(sum(ical),1) Xtrans2(ical,:)] * b_s;
            [Rc_i, rmsec_i]   = rmse(yhat_c_si, ytrans(ical));
            yhat_si  = [ones(sum(ival),1) Xtrans2(ival,:)] * b_s;
            [R_i, rmsep_i]   = rmse(yhat_si, ytrans(ival));
        elseif method =="MT_PFCE"
            X = inputData{1};
            y = inputData{2};
            Xcal = cell(0);
            ycal = cell(0);
            Xval = cell(0);
            yval = cell(0);
            for k = 1:length(X)
                nsample = length(y{k});
                ival = mod(1:nsample,nFold) ==i-1;
                ical = ~ival;
                Xk = X{k};
                yk = y{k};
                Xcal{k} = Xk(ical,:);
                ycal{k} = yk(ical);
                Xval{k} = Xk(ival,:);
                yval{k} = yk(ival);
            end
            [bglo, Cost_cali] = trans_mt_pfce(Xcal,ycal,ConstType,thresi);
            rmsec_i = [];
            rmsep_i = [];
            for k = 1:length(X)
                nsample = length(y{k});
                ival = mod(1:nsample,nFold) ==i-1;
                ical = ~ival;
                Xk = X{k};
                yk = y{k};
                yhat_c_si  = [ones(sum(ical),1) Xk(ical,:)] * bglo(:,k);
                [Rc_i, rmsec_k_i]   = rmse(yhat_c_si, yk(ical));
                yhat_si  = [ones(sum(ival),1) Xk(ival,:)] * bglo(:,k);
                [R_i, rmsep_k_i]   = rmse(yhat_si, yk(ival));
                rmsec_i = [rmsec_i rmsec_k_i];
                rmsep_i = [rmsep_i rmsep_k_i];
            end
        else
            error('No such method!!!')
        end
        if method == "NS_PFCE" || method == "SS_PFCE" || method == "FS_PFCE"
            RMSECal(i,j)=rmsec_i;
            RMSECVal(i,j)=rmsep_i;
        elseif method == "MT_PFCE"
            RMSECal(i,j,:)=rmsec_i;
            RMSECVal(i,j,:)=rmsep_i;
        else
            error("The inputted method in not supportted.")
        end
    end
end
[minCost,idxMinCost] = min(squeeze(mean(RMSECVal)));
optThres = thresRange(idxMinCost);
end
