function [converged,c] = perform_convergence_STM_Stijn(string,x,i,converged,nSteps,nSamples,model,c_threshold,c)

% Get data from histogram function    
    if i > 2
        for l = 1:size(x.(string),1)
            x_string_names = fieldnames(x);            
            x1 = x.(x_string_names{i-2})(l,:);
            x2 = x.(x_string_names{i-1})(l,:);
            x3 = x.(x_string_names{i})(l,:);
            
            hist = histogram(x1,'NumBins',30,'Normalization','probability');
            hist_data.x1.Values = hist.Values;
            hist_data.x1.BinEdges = hist.BinEdges;
            
            hist = histogram(x2,'NumBins',30,'Normalization','probability');
            hist_data.x2.Values = hist.Values;
            hist_data.x2.BinEdges = hist.BinEdges;
            
            hist = histogram(x3,'NumBins',30,'Normalization','probability');
            hist_data.x3.Values = hist.Values;
            hist_data.x3.BinEdges = hist.BinEdges;
            
            [h_Values_1_3,p_Values_1_3] = kstest2(hist_data.x1.Values,hist_data.x3.Values);
            [h_Values_2_3,p_Values_2_3] = kstest2(hist_data.x2.Values,hist_data.x3.Values);
            [h_BinEdges_2_3,p_BinEdges_2_3] = kstest2(hist_data.x2.BinEdges,hist_data.x3.BinEdges);
            
            if p_Values_1_3 < 0.95 && p_Values_2_3 > 0.95 && p_BinEdges_2_3 > 0.95 && converged(l,1) == 0
                converged(l,1) = 1;
                converged(l,2) = p_Values_2_3;
                converged(l,3) = p_BinEdges_2_3;
                converged(l,4) = nSamples;
                converged(l,5) = nSteps;
                fprintf('Convergence for reaction %i achieved at nSteps = %i and nSamples = %i.\n',l,nSteps,nSamples);
            end
        end
    end
    
    close all
    
    % Print summary of converged reactions
    fprintf('%i reactions converged with %i as number of target reactions (%i total reactions).\n',sum(converged(:,1)),ceil(c_threshold*size(x.(string),1)),size(model.rxns,1));
    
    % Change the convergence tracker (c) to 1 if a sufficient amount of
    % reactions are converged.
    if sum(converged(:,1)) >= ceil(c_threshold*size(x.(string),1))
        fprintf('Convergence achieved for at least %i percent of the reactions.\n',c_threshold*100);
        c = 1;
        return
    end