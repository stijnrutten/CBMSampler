function [converged,c,x_max] = perform_convergence_CSM_Stijn(accuracy_range,string,model,x,x_string_names,i,converged,c_threshold,c,x_max,nSteps,nSamples)

% Get data from histogram fit function
figure
for j = 1:size(model.rxns,1)
    [x_fit.(string)(j,:),y_fit.(string)(j,:)] = CustomHistfit(x.(string)(j,:),30,'ev');
end
close figure 1

% Find maximum probability value of histogram function (y-axis = probability, x-axis = flux)
[y_max.(string),y_max.(string)(:,2)] = max(y_fit.(string),[],2);

% Find maximum flux value of histogram function (y-axis = probability, x-axis = flux)
for k = 1:size(x_fit.(string),1)
    x_max.(x_string_names{i})(k,1) = x_fit.(x_string_names{i})(k,y_max.(x_string_names{i})(k,2));
end

% Check if the flux value at maximum probability at the current
% calculation deviates less than a given percentage (2*0.025 = 5%,
% since data could be positive or negative. If the flux value lies
% within the given percentage, allocate 1 to the 'converged' vector of
% that reation.

if i > 1
    for l = 1:size(x_fit.(string),1)
        if converged(l,1) == 0
            if x_max.(x_string_names{i-1})(l) ~= 0
                old = (accuracy_range/2)*abs(x_max.(x_string_names{i-1})(l));
            elseif x_max.(x_string_names{i-1})(l) == 0
                old = (accuracy_range/2)*abs(max(x_fit.(string)(l,:)) - min(x_fit.(string)(l,:)));
            end
            range_min = x_max.(x_string_names{i-1})(l)-old;
            range_max = x_max.(x_string_names{i-1})(l)+old;
            if x_max.(x_string_names{i})(l) > range_min && x_max.(x_string_names{i})(l) < range_max && converged(l,1) == 0
                converged(l,1) = 1;
                converged(l,2) = nSamples;
                converged(l,3) = nSteps;
                fprintf('Convergence for reaction %i achieved at nSteps = %i and nSamples = %i.\n',l,nSteps,nSamples);
            end
        end
    end
end

% Print summary of converged reactions
fprintf('%i reactions converged with %i as number of target reactions (%i total reactions).\n',sum(converged(:,1)),ceil(c_threshold*size(x_fit.(string),1)),size(model.rxns,1));

% Change the convergence tracker (c) to 1 if a sufficient amount of
% reactions are converged.
if sum(converged(:,1)) >= ceil(c_threshold*size(x_fit.(string),1))
    fprintf('Convergence achieved for at least %i percent of the reactions.\n',c_threshold*100);
    c = 1;
    return
end