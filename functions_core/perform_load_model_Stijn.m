function [model] = perform_load_model_Stijn(modelPath)
% Load model

model_names = dir(modelPath);

fprintf('Choose the desired model from the models folder:\n');
for i = 1:size(model_names,1)
    if model_names(i).isdir == 0
        fprintf('[%i] %s\n',i,model_names(i).name);
    end
end

prompt_model = 'Enter model number: ';
model_nr = input(prompt_model);

modelName = strcat(model_names(model_nr).name);
load(modelName);