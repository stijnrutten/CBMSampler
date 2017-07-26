function gurobiOptions = GurobiOptions(varargin)

%% Parameters

doubleParams = {...
%   Parameter         min   max
    'feasibilityTol', 1e-9, 1e-2
    'optimalityTol' , 1e-9, 1e-2
    'timeLimit'     , 0   , realmax
    };

intParams = {...
%   Parameter     min  max
    'scaleFlag' , 0  , 2
    'method'    , -1 , 4
    'outputFlag', 0  , 1
    };

%% Intialize

numDoubleParams = size(doubleParams, 1);

params = lower([doubleParams(:, 1); intParams(:, 1)]);

%% Create options struct

gurobiOptions = struct();

for iArg = 1:2:(nargin - 1)
    field = varargin{iArg};
    
    if ischar(field)
        iParam = find(ismember(params, lower(field)));
        
        if ~isempty(iParam) % Valid parameter
            if iParam <= numDoubleParams % Double parameter
                field = doubleParams{iParam, 1};
                lb = doubleParams{iParam, 2};
                ub = doubleParams{iParam, 3};

                val = varargin{iArg + 1};
            else % Int parameter
                iParam = iParam - numDoubleParams;
                
                field = intParams{iParam, 1};
                lb = int32(intParams{iParam, 2});
                ub = int32(intParams{iParam, 3});

                val = int32(varargin{iArg + 1});
            end

            if val > ub
                warning('Inputted value of parameter ''%s'' exceeds its upper bound.', ...
                    field);
            elseif val < lb
                warning('Inputted value of parameter ''%s'' exceeds its lower bound.', ...
                    field);
            else % Valid parameter value
                gurobiOptions.(field) = val;
            end
        else % Invalid parameter
            warning('Unrecognized parameter ''%s''.', field);
        end
    end
end

end