% More information about building c++ programs with Gurobi:
% https://www.gurobi.com/documentation/6.5/quickstart_mac/cpp_building_and_running_t.html

clearvars; clc;

%% Input

armaDir = 'C:\armadillo-7.950.1'; % Armadillo directory
gurobiDir = 'C:\gurobi751'; % Gurobi directory

msvcVersion = '2015'; % Version of msvc with which the Gurobi libraries are compiled
gurobiVersion = '75';

%% MEX command

% Path to the current directory:
dirCur = fileparts(mfilename('fullpath'));

% C++ files:
cppFiles = {'CBMSamplerMEX'; 'Utility'; 'SamplerFVA'; 'Sampler'};
cppFiles = strjoin(cellfun(@(x) sprintf(' "%s\\src\\%s.cpp"', dirCur, x), ...
	cppFiles, 'UniformOutput', false), '');

% C++ include directory:
includeDir = sprintf(' -I"%s\\include"', dirCur);

% Armadillo include directories:
armaIncludeDir = sprintf(' -I"%s\\include"', armaDir);

% Gurobi include directories:
gurobiIncludeDir = sprintf(' -I"%s\\win64\\include"', gurobiDir);

% Gurobi library:
gurobiLibDir = sprintf(' -L"%s\\win64\\lib"', gurobiDir);
gurobiLibs = sprintf(' -lgurobi_c++mt%s -lgurobi%s', msvcVersion, gurobiVersion);

% Blas/Lapack libraries:
blasLapackLibDir = sprintf(' -L"%s\\lib_win64"', dirCur);
blasLapackLibs = ' -lblas_win64_MT -llapack_win64_MT';

% Matlab libraries:
matlabLibDir = sprintf(' -L"%s"', fullfile(matlabroot, 'bin', computer('arch')));
matlabLibs = ' -lmwblas -lmwlapack -llibut'; % These libraries are shipped with MATLAB
% TODO: Perform svd and matrix multiplication using Armadillo. In that
% case, we don't need the functions "svd" and "MatrixMultiply" defined in
% "Sampler.cpp", which depend on the above MATLAB libraries. We can then
% remove these dependencies.

% Output directory:
dirOutput = sprintf('%s\\output', dirCur);

% Command:
cmd = ['mex -v -largeArrayDims', ...
	includeDir, armaIncludeDir, gurobiIncludeDir, ...
	cppFiles, ...
	gurobiLibDir, gurobiLibs, ... % Gurobi libraries
	blasLapackLibDir, blasLapackLibs, ... % Blas and Lapack libraries
	matlabLibDir, matlabLibs, ... % MATLAB libraries
	' COMPFLAGS="$COMPFLAGS /MT /openmp"', ...
	' -outdir "', dirOutput, '"']; % Output

%% Compile MEX

eval(cmd); % Execute command