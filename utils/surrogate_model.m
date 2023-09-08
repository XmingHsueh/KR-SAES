function [func_surrogate,model_type] = surrogate_model(database,surrogate_name)

X = database(:,1:end-1);
y = database(:,end);

switch(surrogate_name)
    case 'GPR'
        gprModel = fitrgp(X,y,'Basis','linear','FitMethod','exact','PredictMethod','exact','Standardize',1);
        func_surrogate = @(x)predict(gprModel,x);
        model_type = 'regression';
    case 'GLM'
        glmModel = fitglm(X,y,'purequadratic');
        func_surrogate = @(x)predict(glmModel,x);
        model_type = 'regression';
    case 'EnTr'
        ensembleModel = fitrensemble(X,y,'NumLearningCycles',50);
        func_surrogate = @(x)predict(ensembleModel,x);
        model_type = 'regression';
end