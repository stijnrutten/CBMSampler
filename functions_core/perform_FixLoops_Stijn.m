function [sModel] = perform_FixLoops_Stijn(model,sModel)
% Fix thermodynamically infeasible loops

fprintf('Performing FixLoops...\n');

tic;
[sModel.points_fixed,sModel.fixed] = FixLoops(model,sModel.points);
TSampling = toc;

fprintf('Done! It took %f seconds.\n', TSampling);