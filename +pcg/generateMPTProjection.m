%% Automatic code generation for parametric projection using MPT3 
%
% info = GENERATEPROJECTIONMPT(name, n, 'Aeq', Aeq, 'beq', beq, 'Aineq', Ainq, 'bineq', bineq, 'lb', lb, 'ub', ub, options)
% 

function info = generateMPTProjection(varargin)
    p = inputParser;
    p.addRequired('name', @ischar);
    p.addRequired('n', @(x) (isnumeric(x) && x>0));
    p.addParameter('Aeq', [], @isnumeric);
    p.addParameter('beq', [], @isnumeric);
    p.addParameter('Aineq', [], @isnumeric);
    p.addParameter('bineq', [], @isnumeric);
    p.addParameter('lb', [], @isnumeric);
    p.addParameter('ub', [], @isnumeric);
    p.addParameter('gendir', './gen', @ischar);
    p.addParameter('verbose', 0, @isnumeric);
    p.addParameter('simplify', false, @islogical);
    p.addParameter('lpSolver', 'yalmip', @(s)(any(strcmp(s,{'yalmip', 'cplex', 'gurobi'})))); % LP solver to use in code generation
    p.addParameter('bintree', true, @islogical);
    p.addParameter('infbound', 10000, @isnumeric);
    p.parse(varargin{:});
    options = p.Results;
    
    %% Extract dimensions
    dims.n = options.n;
    dims.p = 0;
    % get dimensions of parameter (if applicable)
    if ~isempty(options.Aeq)
        dims.p = size(options.Aeq,2)-dims.n;
    elseif ~isempty(options.Aineq)
        dims.p = size(options.Aineq,2)-dims.n;
    elseif ~isempty(options.lb)
        dims.p = length(options.lb)-dims.n;
    elseif ~isempty(options.ub)
        dims.p = length(options.ub)-dims.n;
    else
        throw(MException('MATLAB:pcg:generateMPTProjection:Unbounded', 'No constraints defined.'));
    end
    
    %% Setup YALMIP variables
    x = sdpvar(dims.n,1);
    z = sdpvar(dims.n,1);
    J = (x-z)'*(x-z); % Projection objective
    C = [];
    if dims.n > 0 && dims.p > 0 % With additional parameter
        w = sdpvar(dims.p,1);
        y = [w; x];
        theta = [w; z];
    else
        y = x;
        theta = z;
    end
    % Equality constraint
    if ~isempty(options.Aeq)
        C = [C, options.Aeq*y == options.beq];
    end
    % Inequality constraint
    if ~isempty(options.Aineq)
        C = [C, options.Aineq*y <= options.bineq];
    end
    % Lower bound
    if isempty(options.lb)
        options.lb = -inf(dims.p+dims.n,1);
    end
    % Find non-inf indices
    idx = (~isinf(options.lb)).*(1:dims.p+dims.n)';
    idx = idx(idx > 0);
    if ~isempty(idx)
        C = [C, y(idx) >= options.lb(idx)];
    end
    % Upper bound
    if isempty(options.ub)
        options.ub = -inf(dims.p+dims.n,1);
    end
    % Find non-inf indices
    idx = (~isinf(options.ub)).*(1:dims.p+dims.n)';
    idx = idx(idx > 0);
    if ~isempty(idx)
        C = [C, y(idx) <= options.ub(idx)];
    end
    
    %% Construct optimization problem
    if options.simplify
        % Reset MPT3 options and set infinity box
        clear mptopt; mptopt;
        mptopt('infbound', options.infbound, 'verbose', max(0,options.verbose-1));
    end
    problem = Opt(C, J, theta, x);
    
    %% Generate explicit solution
    explicit = problem.solve();
    if options.simplify
        sol = pcg.simplifyProjection(explicit.xopt.Set, options.infbound, 'lpSolver', options.lpSolver);
    else
        sol = explicit.xopt.Set;
    end
    
    %% Generate c-code
    filename = [options.gendir '/' options.name];
    if options.simplify
        try
            tree = BinTreePolyUnion(PolyUnion(sol));
        catch Me
            display(Me);
            tree = BinTreePolyUnion(PolyUnion(explicit.xopt.Set));
        end
        tree.toC('primal', filename);
    elseif options.bintree
        tree = BinTreePolyUnion(explicit.xopt); % PolyUnion(explicit.xopt.Set) is also ok
        tree.toC('primal', filename);
    else
        explicit.xopt.toC('primal', filename);
    end
    
    %%  Process code
    f = fopen([filename '.c'], 'rt');
    code = fread(f, inf, 'char=>char');
    fclose(f);
    % Modify
    code = strrep(code', '#define MPT_NR', '#include "project.h"\n\n#define MPT_NR');
    code = strrep(code, ['static unsigned long ' options.name], ['unsigned long ' options.name]);
    code = strrep(code, ['static long ' options.name], ['unsigned long ' options.name]);
    code = strrep(code, 'static double MPT_ST', 'static const double MPT_ST');
    code = strrep(code, 'static double MPT_F', 'static const double MPT_F');
    code = strrep(code, 'static double MPT_G', 'static const double MPT_G');
    f = fopen([filename '.c'], 'wt');
    fprintf(f, code);
    fclose(f);
    
    %% Cleanup
    delete([filename '_mex.c']);
        
end