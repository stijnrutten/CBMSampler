function [nSamples,nSteps] = perform_increase_parameters_Stijn(c,nSamples,nSteps)
% Increase the amount of samples for the next loop

if c == 0
    nSamples = 100*ceil((1.1*nSamples)/100); % increase nSamples for next loop
    if nSamples > 500
        nSteps = 25*ceil((0.05*nSamples)/25); % increase nSteps for next loop
    end
end