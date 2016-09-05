%% Example 1, automatic code generation for hybrid model predictive control with piecewise affine (PWA) dynamics
%
% [x, u, consensus, t, nIt] = EXAMPLE(options)
%
% Generates code for and runs a closed-loop simulation of a hybrid MPC problem
% OUTPUT:
%  x         - Closed-loop state trajectory [x0 ... xN]
%  u         - Closed-loop input trajectory [u0 ... uN-1]
%  consensus - Final consensus violation ||z-y||^2 for each iteration
%  t         - Execution time for each iteration
%  nIt       - Number of iterations for each iteration
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
% This example is taken from Example 4.1, p. 415 of
% [A. Bemporad and M. Morari, “Control of systems integrating logic, dynamics, and constraints”,
%  Automatica, vol. 35, no. 3, pp. 407–427, 1999.]
%
% The model is parametrized by
%  alphas = [pi/3, -pi/3],
%  beta = 0.8,
% and has dynamics
%  [x1(k+1)]        [cos(alpha(1)) -sin(alpha(1))] [x1(k)]   [0]
%  [x2(k+1)] = beta*[sin(alpha(1))  cos(alpha(1))]*[x2(k)] + [1]*u(k)   if x1(k) >= 0 (Region 1),
% and
%  [x1(k+1)]        [cos(alpha(2)) -sin(alpha(2))] [x1(k)]   [0]
%  [x2(k+1)] = beta*[sin(alpha(2))  cos(alpha(2))]*[x2(k)] + [1]*u(k)   if x1(k) <= 0 (Region 2).
%
% It has state and input (box) constraints
%  xmin <= x(k) <= xmax
%  umin <= u(k) <= umax
%
% And the objective is
%  0.5*x(k)'*Hx*x(k) + hx'*x(k) + 0.5*u(k)'*Hu*u(k) + hu'*u(k)
% with
%       [1 0]       [0]
%  Hx = [0 1], hx = [0], Hu = 1, and hu = 0.
%

function [xs, us, consensus, ts, nIts] = example(varargin)
    %% Defining optional inputs as key-value pairs
    p = inputParser;
    p.addParameter('N', 10, @(x)(isnumeric(x) && length(x) == 1 && x > 0 && mod(x,1) == 0)); % Prediction (or control) horizon N
    % Problem parameters
    p.addParameter('alphas', [pi/3, -pi/3], @(x)(isnumeric(x) && length(x) == 2));
    p.addParameter('beta', 0.8, @(x)(isnumeric(x) && length(x) == 1));
    p.addParameter('xmin', [-10; -10], @(x)(isnumeric(x) && size(x,1) == 2 && size(x,2) == 1)); % State lower bounds
    p.addParameter('xmax', [10; 10], @(x)(isnumeric(x) && size(x,1) == 2 && size(x,2) == 1)); % State upper bounds
    p.addParameter('umin', -1, @(x)(isnumeric(x) && length(x) == 1)); % Input lower bound
    p.addParameter('umax', 1, @(x)(isnumeric(x) && length(x) == 1)); % Input upper bound
    p.addParameter('Hx', eye(2), @(x)(isnumeric(x) && size(x,1) == 2 && size(x,2) == 2 && all(eig(x) > eps))); % Quadratic state cost matrix
    p.addParameter('Hu', 1, @(x)(isnumeric(x) && length(x) == 1 && x > eps));  % Quadratic input cost matrix
    % Iteration parameters
    p.addParameter('xi', 10, @(x)(isnumeric(x))); % Proximal scaling
    % Code-generation parameters
    p.addParameter('solverName', 'hybridMPC', @ischar); % Name of generated solver
    p.addParameter('gendir', './gen', @ischar); % Directory used to place generated solver (and auxiliary files)
    p.addParameter('lpSolver', 'yalmip', @(s)(any(strcmp(s,{'yalmip', 'cplex', 'gurobi'})))); % LP solver to use in code generation
    % Simulation parameters
    p.addParameter('x0', [1; 1]); % Initial state
    p.addParameter('maxIt', 1000, @(x)(isnumeric(x) && length(x) == 1 && x > 0 && mod(x,1) == 0)); % Maximum number of iterations
    p.addParameter('ctol', 1e-3, @(x)(isnumeric(x) && length(x) == 1 && x >= eps)); % Consensus tolerance
    p.addParameter('stoppingCriterion', 'consensus', @(s)(any(strcmp(s,{'consensus', 'maxIt'})))); % Stopping criterion
    % Misc. parameters
    p.addParameter('overwriteSolver', false, @islogical); % Set to true to always overwrite the solver
    p.addParameter('verbose', 2, @isnumeric); % Determines verbosity of output
    p.parse(varargin{:});
    options = p.Results;
    
    if options.verbose >= 1
        display('Setting up example...');
    end
    
    %% Problem parameters
    N = options.N;
    beta  = options.beta;
    alphas = options.alphas;
    xmin = options.xmin;
    xmax = options.xmax;
    umin = options.umin;
    umax = options.umax;
    Hx = options.Hx;
    hx = zeros(2,1);
    Hu = options.Hu;
    hu = zeros(1,1);
    %% Parameters of the method
    xi = options.xi;
    %% Simulation parameters
    x0 = options.x0;
    
    %% Sanity checks
    if any(xmin > xmax)
        error(); % TODO
    end
    if any(umin > umax)
        error(); % TODO
    end
    if any(eig(blkdiag(Hx,Hu)) < 1e-4)
        error();
    end

    %% System dimensions stored in model.dims
    nx = 2; % Number of states
    nu = 1; % Number of controls
    nw = nx; % Number of auxiliary variables (always equal to nx)
    nr = 2; % Number of regions of the PWA dynamics
    n = nx+nu+nw; % Total number of optimization variables (per stage)
    nn = N*n; % Total number of optimization variables
    model.dims = struct('nx',nx,'nu',nu,'nw',nw,'nr',nr,'n',n,'nn',nn);
    
    %% System dynamics stored in model.dyn
    % Define affine dynamics, and constraints for each region
    % Region 1:
    % [x1(k+1)]        [cos(alpha(1)) -sin(alpha(1))] [x1(k)]   [0]
    % [x2(k+1)] = beta*[sin(alpha(1))  cos(alpha(1))]*[x2(k)] + [1]*u(k)   if x1(k) >= 0
    % Region 2:
    % [x1(k+1)]        [cos(alpha(2)) -sin(alpha(2))] [x1(k)]   [0]
    % [x2(k+1)] = beta*[sin(alpha(2))  cos(alpha(2))]*[x2(k)] + [1]*u(k)   if x1(k) <= 0
    for i=1:model.dims.nr
        model.dyn.A{i} = beta*[cos(alphas(i)) -sin(alphas(i)); sin(alphas(i)) cos(alphas(i))];
        model.dyn.B{i} = [0; 1];
        model.dyn.c{i} = zeros(2,1); % No affine part
        % Only upper- and lower bound constraints on the regions,
        %  constraints are defined on [x(k); u(k)]
        switch i
            case 1,
                model.dyn.lb{i} = [0; -inf; -inf];
                model.dyn.ub{i} = [inf; inf; inf];
            case 2,
                model.dyn.lb{i} = [-inf; -inf; -inf];
                model.dyn.ub{i} = [0; inf; inf];
        end
        model.dyn.Aineq{i} = {};
        model.dyn.bineq{i} = {};
        model.dyn.Aeq{i} = {};
        model.dyn.beq{i} = {};
    end
    % TODO: Needed?
    for i=1:model.dims.nr
        model.dyn.isin{i} = @(x,u) (all([x;u] >= model.dyn.lb{i} & [x;u] <= model.dyn.ub{i}));
    end
    
    %% Additional state and input constraints for all regions stored in model.const
    model.const = struct('Aineq', zeros(0,model.dims.nx+model.dims.nu), ...
                         'bineq', zeros(0,1), ...
                         'Aeq', zeros(0,model.dims.nx+model.dims.nu), ...
                         'beq', zeros(0,1), ...
                         'lb', [xmin; umin], ...
                         'ub', [xmax; umax]);
	% TODO: Needed?
    model.const.check = @(u, x, e) (all([x;u] >= model.const.lb-e & [x;u] <= model.const.ub+e));
    
    %% Optimization objective (strictly convex quadratic)
    
    %% Generate embedded solver
    [update, xi] = pcg.generateSolver(model, N, Hx, hx, Hu, hu, xi, ...
                                                   'solverName', options.solverName, ...
                                                   'overwriteSolver', options.overwriteSolver, ...
                                                   'gendir', options.gendir, ...
                                                   'stoppingCriterion', options.stoppingCriterion, ...
                                                   'maxIt', options.maxIt, ...
                                                   'ctol', options.ctol, ...
                                                   'lpSolver', options.lpSolver); %#ok
    %% Get function to update  iteration data
    
    %% Run
    % Add path to fixed-point iteration code
    addpath(options.gendir);
    % Warm-up projection
    eval([options.solverName '(zeros(model.dims.nn,1), zeros(model.dims.n,1), zeros(model.dims.nn,1));']);
    
    nSteps = 10;
    % Allocate memory
    nIts = zeros(nSteps,1);
    ts = zeros(nSteps,1);
    consensus = zeros(nSteps,1);
    lambda = zeros(model.dims.nn,nSteps);
    obj = nan(nSteps,N);
    cviol = zeros(nSteps);
    
    % Generate references
    ref.x = zeros(model.dims.nx,200);
    ref.u = zeros(model.dims.nu,200);
    xs = [x0 zeros(model.dims.nx,nSteps)];
    us = zeros(model.dims.nu,nSteps);
    
    for k=1:nSteps
        if options.verbose >= 1
            display(['Time k=' num2str(k) ', running embedded iteration...']);
        end
        % Setting up prox-iteration
        best = inf; %#ok
        % Initializing
        refs = reshape([ref.u(:,k:N+k-1); ref.x(:,k:N+k-1); ref.x(:,k:N+k-1)], [], 1); %#ok
        % Run Fixed-point iteration
        [z, y, s, nIt, t] = eval([options.solverName '(zeros(model.dims.nn, 1), xs(:,k), update(refs));']);
        nIts(k) = nIt;
        ts(k) = t;
        consensus(k) = norm(z-y,2);
        sol = y;
        lambda(:,k) = xi*(sol-s); % Re-constructing multipliers
        switch options.stoppingCriterion
            case 'consenus',
                if nIt >= options.maxIt
                    warning('simulate:Prox:MaximumInterations', ['Embedded Prox: Maximum iterations limit exceeded, consensus is ' num2str(consensus(k)) '.']);
                end
            case 'maxIt',
                if consensus(k) > options.ctol
                    warning('simulate:Prox:InsufficientAccuracy', ['Embedded Prox: Consensus specification violated, consensus is ' num2str(consensus(k)) '.']);
                end
        end
        [x,u,w] = pcg.recoverStatesInputs(sol, model.dims, 'auxiliaries', true);
        
        cviol(k) = norm(x-w,2);
        for l=1:N
            obj(k,l) = 0.5*(x(:,l)-ref.x(:,l))'*Hx*(x(:,l)-ref.x(:,l)) + 0.5*(u(:,l)-ref.u(:,l))'*Hu*(u(:,l)-ref.u(:,l));
        end
        
        % Simulate
        us(:,k) = u(:,1);
        if xs(1,k) >= 0
            xs(:,k+1) = model.dyn.A{1}*xs(:,k) + model.dyn.B{1}*us(:,k) + model.dyn.c{1};
        else
            xs(:,k+1) = model.dyn.A{2}*xs(:,k) + model.dyn.B{2}*us(:,k) + model.dyn.c{2};
        end
    end
    
    if options.verbose >= 2
        figure();
        subplot(2,1,1); plot(xs(1,:),xs(2,:)); axis equal; xlabel('state x_1'); ylabel('state x_2');
        subplot(2,1,2); plot(us); xlabel('time step k'); ylabel('control input u');
        figure();
        semilogy(ts*1000); xlabel('time step k'); ylabel('Execution time in [ms]');
    end

end
