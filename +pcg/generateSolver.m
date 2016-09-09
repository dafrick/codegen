%% Automatic code generation for hybrid model predictive control with piecewise affine (PWA) dynamics
%
% [update, xi] = GENERATESOLVER(model, N, Hx, hx, Hu, hu, options)
%
% min sum_{k=0}^N 0.5*-x_{k+1}'*Hx*x_{k+1} + hx'*x_{k+1} + 0.5*u_{k}'*Hu*u_{k} + hu'*u_{k}
%  s.t. x_{k+1} = dyn.A{i}*x_{k} + dyn.B{i}*u{k} + dyn.c{i},
%                  if (x{k},u{k}) satisfy dyn.lb{i} <= [x{k}; u{k}] <= dyn.ub{i},
%                                         dyn.Aineq{i}*[x{k}; u{k}] <= dyn.bineq{i}, and
%                                         dyn.Aeq{i}[x{k}; u{k}] = dyn.beq{i}.
%       const.lb <= [x{k}; u{k}] <= const.ub,
%       const.Aineq*[x{k}; u{k}] <= const.bineq, and
%       const.Aeq[x{k}; u{k}] = const.beq.
%
% Generates code for a hybrid MPC problem
% INPUT:
%  model - System & controller model. Needs to contain structs:
%           - dims, with the problem dimensions
%              - nx, number of states
%              - nu, number of inputs
%              - nw, number of auxiliaries - needs to be equal to nx
%              - nr, number of regions of the PWA dyamics
%              - n,  number of decision variables per stage (= nx+nu+nw),
%              - nn, number of decision variables (= n*N)
%           - dyn, with the description of the model dynamics and the regions on which they are defined
%              - A, a cell of length nr with a matrix of dimension (nx)x(nx) and
%              - B, a cell of length nr with a matrix of dimension (nx)x(nu) and
%              - c, a cell of length nr with a vector of length nx containing the
%                 region dynamics x(k+1) = Ax(k) + Bu(k) +c
%              - lb, a cell of length nr containing the region lower bounds [xmin; umin]
%              - ub, a cell of length nr containing the region upper bounds [xmax; umax]
%              - Aineq, a cell of length nr with a matrix of dimension (m)*(nx+nu) and
%              - bineq, a cell of length nr with a vector of length m containing the
%                 region inequality constraints Aineq*[x; u] <= bineq
%              - Aeq, a cell of length nr with a matrix of dimension (m)x(nx+nu) and
%              - beq, a cell of length nr with a vector of length m containing the
%                 region equality constraints Aeq*[x; u] = beq
%            - const, with the description of the operating constraints
%              - lb, a vector of lower bounds [xmin; umin]
%              - ub, a vector of upper bounds [xmax; umax]
%              - Aineq, matrix of dimension (m)x(nx+nu) and
%              - bineq, a vector of length m containing the 
%                 inequality constraints Aineq*[x; u] <= bineq
%              - Aeq, a matrix of dimension (m)x(nx+nu) and
%              - beq, a vector of length m containing the
%                 equality constraints Aeq*[x; u] = beq
%  N - Control horizon
%  Hx - Quadratic state-cost matrix
%  hx - Linear state-cost
%  Hu - Quadratic input-cost matrix
%  hu - Linear input-cost
% OUTPUT:
%  update
% INPUT (optional):
%  'solverName' - Name of the generated solver (Default: hybridMPC)
%  'N'          - Control horizon (Default: 10)
%  'xi'         - Proximal scaling (Default: 10)
%  'x0'         - Initial state (Default: [1;1])
%  'stoppingCriterion' - Criterion used to terminate solver, can be either 'consensus' or 'maxIt' (Default: consensus)
%  'tol'        - Consensus tolerance, stopping criterion if stoppingCriterion is 'consensus' (Default: 1e-3)
%  'maxIt'      - Maximum number of iterations, stopping criterion if stoppingCriterion is 'maxIt' (Default: 1000)
%  'gendir'     - Directory used to place generated code (Default: gen)
%  'overwriteSolver'   - Whether solver should be overwritten (Default: true)
%  'verbose'    - Determines verbosity of output, 0 is no output, 2 is full output (Default: 2)
%

function [update, xi] = generateSolver(varargin)
    p = inputParser;
    p.addRequired('model', @isstruct);
    p.addRequired('N', @(x)(isnumeric(x) && length(x) == 1 && x > 0 && mod(x,1) == 0));
    p.addRequired('Hx');
    p.addRequired('hxlin');
    p.addRequired('Hu');
    p.addRequired('hulin');
    p.addRequired('xi'); % Proximal scaling
    p.addParameter('solverName', 'hybridMPC', @ischar); % Name of generated solver
    p.addParameter('gamma', 0.5); % Iteration step size
    p.addParameter('maxIt', 1000, @(x)(isnumeric(x) && length(x) == 1 && x > 0 && mod(x,1) == 0)); % Maximum number of iterations
    p.addParameter('ctol', 1e-3, @(x)(isnumeric(x) && length(x) == 1 && x >= eps)); % Consensus tolerance
    p.addParameter('stoppingCriterion', 'consensus', @(s)(any(strcmp(s,{'consensus', 'maxIt'})))); % Stopping criterion
    p.addParameter('overwriteSolver', false, @islogical); % Whether the solver should be overwrittem
    p.addParameter('overwriteProjections', true, @islogical); % Whether the projections should be overwritten
    p.addParameter('lpSolver', 'yalmip', @(s)(any(strcmp(s,{'yalmip', 'cplex', 'gurobi'})))); % LP solver to use in code generation
    p.addParameter('gendir', './gen', @ischar); % Directory used to place generated solver (and auxiliary files)
    p.addParameter('verbose', 0, @(x)(isnumeric(x) && length(x) == 1 && x >= 0 && mod(x,1) == 1));  % Determines verbosity of output
    p.parse(varargin{:});
    options = p.Results;
    
    %% Problem data
    model = options.model;
    N = options.N;
    H = kron(eye(N),blkdiag(options.Hu,0.5*options.Hx,0.5*options.Hx));
    h = repmat([options.hulin; 0.5*options.hxlin; 0.5*options.hxlin],N,1);
    xi = options.xi;
    
    if options.verbose >= 1
        display(['Generating solver "' options.solverName '"...']);
    end
    
    %% Generate iteration data
    [M, ~, W, updateC, xi] = pcg.getIterationData(H, h, N, model.dims, xi);
    update = @(r)(updateC(h-H'*r));
    
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
            if options.overwriteProjections
                pcg.generateProjections(model, N, 'verbose', options.verbose, 'gendir', options.gendir, 'lpSolver', options.lpSolver);
            end
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