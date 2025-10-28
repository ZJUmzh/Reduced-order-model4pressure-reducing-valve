clc
rng(0)
%load('svmmdl.mat')
for i = 1:size(An,1)
    close all
svmmdl(i).model = fitrsvm(ph,An(:,i) ,'OptimizeHyperparameters', 'none', ...
    'HyperparameterOptimizationOptions', struct('AcquisitionFunctionName', 'expected-improvement-plus', ...
        'MaxTime', 300),'KernelFunction','linear','Standardize',true,'Epsilon',svmmdl(i).model.Epsilon*0.1,...
        'BoxConstraint',svmmdl(i).model.BoxConstraints(1),'KernelScale',svmmdl(i).model.KernelParameters.Scale);
end
%svmmdl(i).model.Epsilon

save('E:\study materials\2024\2024autumn\PRV\matlab_code\newdata\fitting\svmmdl_new.mat','svmmdl');