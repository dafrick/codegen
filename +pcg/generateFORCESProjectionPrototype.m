%% Automatic code generation for parametric projection using MPT3 
%
% info = GENERATEFORCESPROJECTION(name, n, 'Aeq', Aeq, 'beq', beq, 'Aineq', Ainq, 'bineq', bineq, 'lb', lb, 'ub', ub, options)
% 

function info = generateFORCESProjectionPrototype(varargin)
    p = inputParser;
    p.addRequired('name', @ischar);
    p.addParameter('B', [], @(x)(isnumeric(x) || iscell(x)))
    p.addParameter('A', [], @(x)(isnumeric(x) || iscell(x)))
    p.addParameter('c', [], @(x)(isnumeric(x) || iscell(x)))
    p.addParameter('Aeq', [], @(x)(isnumeric(x) || iscell(x)));
    p.addParameter('beq', [], @(x)(isnumeric(x) || iscell(x)));
    p.addParameter('Aineq', [], @(x)(isnumeric(x) || iscell(x)));
    p.addParameter('bineq', [], @(x)(isnumeric(x) || iscell(x)));
    p.addParameter('lb', [], @(x)(isnumeric(x) || iscell(x)));
    p.addParameter('ub', [], @(x)(isnumeric(x) || iscell(x)));
    p.addParameter('parametrize', struct(), @isstruct);
    p.addParameter('gendir', './gen', @ischar);
    p.addParameter('verbose', 0, @isnumeric);
    p.parse(varargin{:});
    options = p.Results;
    
    %% Sanity check
    if isempty(options.Aeq) && isempty(options.Aineq) && isempty(options.lb) && isempty(options.ub)
        throw(MException('MATLAB:pcg:generateFORCESProjectionPrototype:Unbounded', 'No constraints defined.'));
    end
    % TODO: Add more sanity checks
    
    %% Setup FORCES problem
    % Distinguish between two cases
    paramIdx = 0;
    if ~isempty(options.A)
        stages = MultistageProblem(2);
        % Find largest necessary "state" dimensions
        if ~iscell(options.A)
            stages(1).dims.n = size(options.A,2);
            stages(2).dims.n = size(options.A,1);
        else
            stages(1).dims.n = 0;
            stages(2).dims.n = 0;
            for i=1:length(options.A)
                stages(1).dims.n = max(stages(1).dims.n, size(options.A{i},2));
                stages(2).dims.n = max(stages(2).dims.n, size(options.A{i},1));
            end
        end
        % Cost
        stages(1).cost.H = eye(stages(1).dims.n);
        paramIdx = paramIdx+1; parameters(paramIdx) = newParam('x', 1, 'cost.f');
        stages(2).cost.H = eye(stages(2).dims.n);
        paramIdx = paramIdx+1; parameters(paramIdx) = newParam('w', 2, 'cost.f');
        % Dynamics
        stages(2).dims.r = stages(2).dims.n;
        if ~iscell(options.A)
            stages(1).C = options.A;
        else
            paramIdx = paramIdx+1; parameters(paramIdx) = newParam('A', 1, 'eq.C');
        end
        if isempty(options.B)
            stages(2).eq.D = -eye(stages(2).dims.n);
        elseif ~iscell(options.B)
            stages(2).eq.D = -options.B;
        else
            paramIdx = paramIdx+1; parameters(paramIdx) = newParam('B', 2, 'eq.D');
        end
        if isempty(options.c)
            stages(2).eq.c = zeros(stages(2).dims.n,1);
        elseif ~iscell(options.c)
            stages(2).eq.c = -options.c;
        else
            paramIdx = paramIdx+1; parameters(paramIdx) = newParam('c', 2, 'eq.c');
        end
        % Equality constraints
        if isempty(options.Aeq)
            stages(1).dims.r = 0;
        elseif ~iscell(options.Aeq)
            stages(1).dims.r = size(options.Aeq,1);
            stages(1).eq.D = options.Aeq;
            stages(1).eq.c = options.beq;
        else
            stages(1).dims.r = 0;
            for i=1:length(options.Aeq)
                stages(1).dims.r = max(stages(1).dims.r, size(options.Aeq{i},1));
            end
            paramIdx = paramIdx+1; parameters(paramIdx) = newParam('Aeq', 1, 'eq.D');
            paramIdx = paramIdx+1; parameters(paramIdx) = newParam('beq', 1, 'eq.c');
        end
        % Inequality constraints
        if isempty(options.Aineq)
            stages(1).dims.p = 0;
        elseif ~iscell(options.Aineq)
            stages(1).dims.p = size(options.Aineq,1);
            stages(1).ineq.p.A = options.Aineq;
            stages(1).ineq.p.b = options.bineq;
        else
            stages(1).dims.p = 0;
            for i=1:length(options.Aineq)
                stages(1).dims.p = max(stages(1).dims.p, size(options.Aineq{i},1));
            end
            paramIdx = paramIdx+1; parameters(paramIdx) = newParam('Aineq', 1, 'ineq.p.A');
            paramIdx = paramIdx+1; parameters(paramIdx) = newParam('bineq', 1, 'ineq.p.b');
        end
        stages(2).dims.p = 0;
        % Lower bounds
        if isempty(options.lb)
            stages(1).dims.l = 0;
        elseif ~iscell(options.lb)
            stages(1).dims.l = sum(~isinf(options.lb));
            stages(1).ineq.b.lbidx = subsref((1:stages(1).dims.n)',struct('type', '()', 'subs', {{~isinf(options.lb)}})); % index vector for lower bounds
            stages(1).ineq.b.lb = options.lb(stages(1).ineq.b.lbidx); % lower bounds
        else
            stages(1).dims.l = 0;
            stages(1).ineq.b.lbidx = [];
            for i=1:length(options.lb)
                stages(1).dims.l = max(stages(1).dims.l, sum(~isinf(options.lb{i})));
                stages(1).ineq.b.lbidx = sort(unique([stages(1).ineq.b.lbidx; subsref((1:length(options.lb{i}))',struct('type', '()', 'subs', {{~isinf(options.lb{i})}}))]));
            end
            paramIdx = paramIdx+1; parameters(paramIdx) = newParam('lb', 1, 'ineq.b.lb');
        end
        stages(2).dims.l = 0;
        % Upper bounds
        if isempty(options.ub)
            stages(1).dims.u = 0;
        elseif ~iscell(options.ub)
            stages(1).dims.u = sum(~isinf(options.ub));
            stages(1).ineq.b.ubidx = subsref((1:stages(1).dims.n)',struct('type', '()', 'subs', {{~isinf(options.ub)}})); % index vector for lower bounds
            stages(1).ineq.b.ub = options.ub(stages(1).ineq.b.ubidx); % lower bounds
        else
            stages(1).dims.u = 0;
            stages(1).ineq.b.ubidx = [];
            for i=1:length(options.ub)
                stages(1).dims.u = max(stages(1).dims.u, sum(~isinf(options.ub{i})));
                stages(1).ineq.b.ubidx = sort(unique([stages(1).ineq.b.ubidx; subsref((1:length(options.ub{i}))',struct('type', '()', 'subs', {{~isinf(options.ub{i})}}))]));
            end
            paramIdx = paramIdx+1; parameters(paramIdx) = newParam('ub', 1, 'ineq.b.ub');
        end
        stages(2).dims.u = 0;
        % Quadtratic constraints
        stages(1).dims.q = 0;
        stages(2).dims.q = 0;
        
        % Set outputs
        outputs(1) = newOutput('z',1,1:stages(1).dims.n);
        outputs(2) = newOutput('u',2,1:stages(2).dims.n);
    else
        stages = MultistageProblem(1);
        % Find largest necessary "state" dimensions
        stages(1).dims.n = 0;
        if ~isempty(options.Aeq)
            if ~iscell(options.Aeq)
                stages(1).dims.n = max(stages(1).dims.n, size(options.Aeq,2));
            else
                for i=1:length(options.Aeq)
                    stages(1).dims.n = max(stages(1).dims.n, size(options.Aeq{i},2));
                end
            end
        end
        if ~isempty(options.Aineq)
            if ~iscell(options.Aineq)
                stages(1).dims.n = max(stages(1).dims.n, size(options.Aineq,2));
            else
                for i=1:length(options.Aineq)
                    stages(1).dims.n = max(stages(1).dims.n, size(options.Aineq{i},2));
                end
            end
        end
        if ~isempty(options.lb)
            if ~iscell(options.lb)
                stages(1).dims.n = max(stages(1).dims.n, length(options.lb));
            else
                for i=1:length(options.lb)
                    stages(1).dims.n = max(stages(1).dims.n, length(options.lb{i}));
                end
            end
        end
        if ~isempty(options.ub)
            if ~iscell(options.ub)
                stages(1).dims.n = max(stages(1).dims.n, length(options.ub));
            else
                for i=1:length(options.ub)
                    stages(1).dims.n = max(stages(1).dims.n, length(options.ub{i}));
                end
            end
        end
        % Cost
        stages(1).cost.H = eye(stages(1).dims.n);
        paramIdx = paramIdx+1; parameters(paramIdx) = newParam('x', 1, 'cost.f');
        % Equality constraints
        if isempty(options.Aeq)
            stages(1).dims.r = 0;
        elseif ~iscell(options.Aeq)
            stages(1).dims.r = size(options.Aeq,1);
            stages(1).eq.D = options.Aeq;
            if ~isfield(options.parametrize, 'beq')
                stages(1).eq.c = options.beq;
            else
                paramIdx = paramIdx+1; parameters(paramIdx) = newParam('beq', 1, 'eq.c');
            end
        else
            stages(1).dims.r = 0;
            for i=1:length(options.Aeq)
                stages(1).dims.r = max(stages(1).dims.r, size(options.Aeq{i},1));
            end
            if stages(1).dims.r > 0
                paramIdx = paramIdx+1; parameters(paramIdx) = newParam('Aeq', 1, 'eq.D');
                paramIdx = paramIdx+1; parameters(paramIdx) = newParam('beq', 1, 'eq.c');
            end
        end
        % Inequality constraints
        if isempty(options.Aineq)
            stages(1).dims.p = 0;
        elseif ~iscell(options.Aineq)
            stages(1).dims.p = size(options.Aineq,1);
            stages(1).ineq.p.A = options.Aineq;
            if ~isfield(options.parametrize, 'beq')
                stages(1).ineq.p.b = options.bineq;
            else
                paramIdx = paramIdx+1; parameters(paramIdx) = newParam('bineq', 1, 'ineq.p.b');
            end
        else
            stages(1).dims.p = 0;
            for i=1:length(options.Aineq)
                stages(1).dims.p = max(stages(1).dims.p, size(options.Aineq{i},1));
            end
            if stages(1).dims.p > 0
                paramIdx = paramIdx+1; parameters(paramIdx) = newParam('Aineq', 1, 'ineq.p.A');
                paramIdx = paramIdx+1; parameters(paramIdx) = newParam('bineq', 1, 'ineq.p.b');
            end
        end
        % Lower bounds
        if isempty(options.lb)
            stages(1).dims.l = 0;
        elseif ~iscell(options.lb)
            stages(1).dims.l = sum(~isinf(options.lb));
            stages(1).ineq.b.lbidx = subsref((1:stages(1).dims.n)',struct('type', '()', 'subs', {{~isinf(options.lb)}})); % index vector for lower bounds
            stages(1).ineq.b.lb = options.lb(stages(1).ineq.b.lbidx); % lower bounds
        else
            stages(1).dims.l = 0;
            stages(1).ineq.b.lbidx = [];
            for i=1:length(options.lb)
                stages(1).dims.l = max(stages(1).dims.l, sum(~isinf(options.lb{i})));
                stages(1).ineq.b.lbidx = sort(unique([stages(1).ineq.b.lbidx; subsref((1:length(options.lb{i}))',struct('type', '()', 'subs', {{~isinf(options.lb{i})}}))]));
            end
            if stages(1).dims.l > 0
                paramIdx = paramIdx+1; parameters(paramIdx) = newParam('lb', 1, 'ineq.b.lb');
            end
        end
        % Upper bounds
        if isempty(options.ub)
            stages(1).dims.u = 0;
        elseif ~iscell(options.ub)
            stages(1).dims.u = sum(~isinf(options.ub));
            stages(1).ineq.b.ubidx = subsref((1:stages(1).dims.n)',struct('type', '()', 'subs', {{~isinf(options.ub)}})); % index vector for lower bounds
            stages(1).ineq.b.ub = options.ub(stages(1).ineq.b.ubidx); % lower bounds
        else
            stages(1).dims.u = 0;
            stages(1).ineq.b.ubidx = [];
            for i=1:length(options.ub)
                stages(1).dims.u = max(stages(1).dims.u, sum(~isinf(options.ub{i})));
                stages(1).ineq.b.ubidx = sort(unique([stages(1).ineq.b.ubidx; subsref((1:length(options.ub{i}))',struct('type', '()', 'subs', {{~isinf(options.ub{i})}}))]));
            end
            if stages(1).dims.u > 0
                paramIdx = paramIdx+1; parameters(paramIdx) = newParam('ub', 1, 'ineq.b.ub');
            end
        end
        % Quadtratic constraints
        stages(1).dims.q = 0;
        
        % Set outputs
        outputs(1) = newOutput('z',1,1:stages(1).dims.n);
    end
    
    %% Options
    codeoptions = getOptions(options.name);
    % TODO: Make customizable
    codeoptions.printlevel = 0;
    codeoptions.maxit = 20;
    codeoptions.timing = 1;
    codeoptions.floattype = 'double';
    codeoptions.overwrite = 1; % Always overwrite
    codeoptions.BuildSimulinkBlock = 0;
    codeoptions.cleanup = 0;
    
    %% Generate projection using FORCES
    generateCode(stages, parameters, codeoptions, outputs);
    
    info.dims = stages(1).dims;
    info.lbidx = stages(1).ineq.b.lbidx;
    info.ubidx = stages(1).ineq.b.ubidx;
        
end