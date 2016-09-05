function [update] = generateSolver(varargin)
    p = inputParser;
    p.addRequired('model', @isstruct);
    p.addRequired('N', @(x)(isnumeric(x) && length(x) == 1 && x > 0 && mod(x,1) == 0));
    p.addRequired('H');
    p.addRequired('hlin');
    p.addRequired('xi');
    p.addParameter('solverName', 'hybridMPC', @ischar); % Name of generated solver
    p.addParameter('gamma', 0.5);
    p.addParameter('maxIt', 1000, @(x)(isnumeric(x) && length(x) == 1 && x > 0 && mod(x,1) == 0)); % Maximum number of iterations
    p.addParameter('ctol', 1e-3, @(x)(isnumeric(x) && length(x) == 1 && x >= eps)); % Consensus tolerance
    p.addParameter('stoppingCriterion', 'consensus', @(s)(any(strcmp(s,{'consensus', 'maxIt'})))); % Stopping criterion
    p.addParameter('overwriteSolver', false);
    p.addParameter('lpSolver', 'yalmip', @(s)(any(strcmp(s,{'yalmip', 'cplex', 'gurobi'})))); % LP solver to use in code generation
    p.addParameter('gendir', './gen', @ischar);
    p.addParameter('verbose', 0, @(x)(isnumeric(x) && length(x) == 1 && x >= 0 && mod(x,1) == 1));
    p.parse(varargin{:});
    options = p.Results;
    
    %% Problem data
    model = options.model;
    N = options.N;
    H = options.H;
    h = options.hlin;
    xi = options.xi;
    
    %% Generate iteration data
    [M, ~, W, update] = pcg.getIterationData(H, h, N, model.dims, xi);
    
    %% Generate individual projections via MPT
    if ~options.overwriteSolver && exist(options.gendir, 'dir')
        yn = input(['The directory "' options.gendir '" already exists.\nDo you want to overwrite and generate a new solver "' options.solverName '"? [yes|no] '], 's');
        while ~any(strcmp(yn,{'yes', 'no', 'abort'}))
            yn = input('Please type [yes|no|abort] in order to proceed. ', 's');
        end
    else
        yn = 'yes';
    end
    switch yn
        case 'abort',
            error('Aborting.');
        case 'yes',
            warning('off', 'MATLAB:MKDIR:DirectoryExists');
            mkdir(options.gendir);
            warning('on', 'MATLAB:MKDIR:DirectoryExists');
            pcg.generateProjections(model, N, 'verbose', options.verbose, 'gendir', options.gendir, 'lpSolver', options.lpSolver);
            % Generate problem-specific embedded code
            pcg.generateProxCode(M, W, model.dims.nr, 'solverName', options.solverName, ...
                                                      'verbose', options.verbose, ...
                                                      'stepSize', options.gamma, ...
                                                      'ctol', options.ctol, ...
                                                      'maxIt', options.maxIt, ...
                                                      'stoppingCriterion', options.stoppingCriterion);
        case 'no',
    end

end