function [convergence_method] = perform_select_convergence_method_Stijn()
% select convergence method: Statistical test Method (STM) or Common Sense
% Method (CSM)

prompt = 'Which convergence method do you want to use? STM/CSM: ';
str = input(prompt,'s');
if strcmp(str,'STM') == 1 || strcmp(str,'stm') == 1
    convergence_method = 1;
elseif strcmp(str,'CSM') == 1 || strcmp(str,'csm') == 1
    convergence_method = 2;
elseif strcmp(str,'STM') ~= 1 || strcmp(str,'stm') ~= 1 || strcmp(str,'CSM') ~= 1 || strcmp(str,'csm') ~= 1
    convergence_method = 3;
end