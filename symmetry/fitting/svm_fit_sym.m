function An_new = svm_fit_sym(ph_new,ph,An,nmodes)
% load('E:\study materials\2024\2024autumn\PRV\matlab_code\newdata\fitting\svmmdl_new.mat')

for i = 1:nmodes
    close all
    if i==1
    svmmdl(i).model = fitrsvm(ph,An(:,i) ,'OptimizeHyperparameters', 'auto', ...
        'HyperparameterOptimizationOptions', struct('AcquisitionFunctionName', 'expected-improvement-plus', ...
            'MaxTime', 300),'KernelFunction','linear','Standardize',true);
    else 
        params = hyperparameters('fitrsvm',ph,An(:,i));
        params(1).Range = [1e-20,1e20];
        % svmmdl(i).model = fitrsvm(ph,An(:,i) ,'OptimizeHyperparameters', {'Epsilon','BoxConstraint','PolynomialOrder','KernelScale','Standardize'}, ...
        svmmdl(i).model = fitrsvm(ph,An(:,i) ,'OptimizeHyperparameters', params, ...
        'HyperparameterOptimizationOptions', struct('AcquisitionFunctionName', 'expected-improvement-plus', ...
            'MaxTime', 300),'KernelFunction','polynomial','Standardize',true);
    end

% ypred = resubPredict(gprmdl);

An_new(:,i) = predict(svmmdl(i).model,ph_new);
end

% for i = 1:nmodes
%     close all
%     if i==1
%         params = hyperparameters('fitrsvm',ph,An(:,i));
%         params(1).Range = [1e-20,1e20];
%     svmmdl(i).model = fitrsvm(ph,An(:,i) ,'OptimizeHyperparameters', params, ...
%         'HyperparameterOptimizationOptions', struct('AcquisitionFunctionName', 'expected-improvement-plus', ...
%             'MaxTime', 300),'KernelFunction','linear','Standardize',true);
%     else 
%         params = hyperparameters('fitrsvm',ph,An(:,i));
%         params(1).Range = [1e-20,1e20];
%         % svmmdl(i).model = fitrsvm(ph,An(:,i) ,'OptimizeHyperparameters', {'Epsilon','BoxConstraint','PolynomialOrder','KernelScale','Standardize'}, ...
%         svmmdl(i).model = fitrsvm(ph,An(:,i) ,'OptimizeHyperparameters', params, ...
%         'HyperparameterOptimizationOptions', struct('AcquisitionFunctionName', 'expected-improvement-plus', ...
%             'MaxTime', 300),'KernelFunction','rbf','Standardize',true);
%     end
% 
% % ypred = resubPredict(gprmdl);
% 
% An_new(:,i) = predict(svmmdl(i).model,ph_new);
% end

% for i = 1:nmodes
%     close all
%     if i==1
%         params = hyperparameters('fitrsvm',ph,An(:,i));
%         params(1).Range = [1e-20,1e20];
%     svmmdl(i).model = fitrsvm(ph,An(:,i) ,'OptimizeHyperparameters', params, ...
%         'HyperparameterOptimizationOptions', struct('AcquisitionFunctionName', 'expected-improvement-plus', ...
%             'MaxTime', 300),'KernelFunction','polynomial','Standardize',true);
%     else 
%         params = hyperparameters('fitrsvm',ph,An(:,i));
%         params(1).Range = [1e-20,1e20];
%         params(2).Range = [1e-4,1e4];
%         % svmmdl(i).model = fitrsvm(ph,An(:,i) ,'OptimizeHyperparameters', {'Epsilon','BoxConstraint','PolynomialOrder','KernelScale','Standardize'}, ...
%         svmmdl(i).model = fitrsvm(ph,An(:,i) ,'OptimizeHyperparameters', params, ...
%         'HyperparameterOptimizationOptions', struct('AcquisitionFunctionName', 'expected-improvement-plus', ...
%             'MaxTime', 300),'KernelFunction','rbf','Standardize',true);
%     end
% 
% % ypred = resubPredict(gprmdl);
% 
% An_new(:,i) = predict(svmmdl(i).model,ph_new);
% end

save('svmmodel.mat','svmmdl');
end