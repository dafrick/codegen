function [s, name] = generateFORCESParameters(varargin)
    p = inputParser;
    p.addRequired('struct', @ischar);
    p.addRequired('set', @isstruct);
    p.addRequired('parametric', @isstruc);
    p.addParameter('name', 'params', @ischar);
    p.addParameter('suffix', '', @ischar);
    p.addParameter('comment', '', @ischar);
    p.addParameter('big', 1e3, @(x)(isnumeric(x) && length(x) == 1 && x > 0));
    p.parse(varargin{:});
    options = p.Results;
    
    set = options.set;
        
    n = options.parametric.n;
    
    vec2strjoin = @(v)(pcg.utils.vec2strjoin(v, ', ', 'format', '%0.17f'));
   
    % Struct name
    name = [options.name, options.suffix];
    
    s = [];
    if isempty(options.comment)
        s = [s, sprintf(['\t' options.struct '_params %s = {\n'], name)];
    else
        s = [s, sprintf(['\t' options.struct '_params %s = { /* %s */\n'], name, options.comment)];
    end
    % Projectant
    s = [s, sprintf(['\t' '\t' '\t{' strjoin(mat2cell(repmat('0.0',1,n),1,3*ones(n,1)), ', ') '}, /* x (projectant) */\n'])];
    % Equality constraint
    if isfield(options.parametric, 'Aeq')
        Aeq = [set.Aeq zeros(size(set.Aeq,1), n-size(set.Aeq,2)); zeros(size(options.parametric.Aeq,1)-size(set.Aeq,1),n)];
        beq = [set.beq; zeros(size(options.parametric.Aeq,1)-length(set.beq),1)];
        s = [s, sprintf(['\t' '\t' '\t{' sprintf('%0.17g, ', Aeq(:,1)) ' /* Aeq (column major) */\n'])];
        for j=2:n-1
            s = [s, sprintf(['\t' '\t' '\t' ' ' sprintf('%0.17g, ', Aeq(:,j)) '\n'])]; %#ok
        end
        s = [s, sprintf(['\t' '\t' '\t' ' ' vec2strjoin(Aeq(:,end)) '},' '\n'])];
        s = [s, sprintf(f, ['\t' '\t' '\t{' vec2strjoin(beq) '}, /* beq */\n'])];
    end
    % Inequality constraint
    if isfield(options.parametric, 'Aineq')
        Aineq = [set.Aineq zeros(size(set.Aineq,1), n-size(set.Aineq,2)); zeros(size(parametric.Aineq,1)-size(set.Aineq,1),n)];
        bineq = [set.bineq; zeros(size(parametric.Aineq,1)-length(set.bineq),1)];
        s = [s, sprintf(['\t' '\t' '\t{' sprintf('%0.17g, ', Aineq(:,1)) ' /* Aineq (column major) */\n'])];
        for j=2:n-1
            s = [s, sprintf(['\t' '\t' '\t' ' ' sprintf('%0.17g, ', Aineq(:,j)) '\n'])]; %#ok
        end
        s = [s, sprintf(['\t' '\t' '\t' ' ' vect2strjoin(Aineq(:,end)) '},' '\n'])];
        s = [s, sprintf(['\t' '\t' '\t{' vect2strjoin(bineq) '}, /* bineq */\n'])];
    end
    % Lower bound
    if isfield(options.parametric, 'lb')
        lb = [set.lb; -inf(n-length(set.lb),1)];
        lb(isinf(lb)) = options.large;
        s = [s, sprintf(['\t' '\t' '\t{' vec2strjoin(lb(options.parametric.lb)) '}, /* lb */\n'])];
    end
    % Upper bound
    if isfield(options.parametric, 'ub')
        ub = [set.ub; -inf(n-length(set.ub),1)];
        ub(isinf(ub)) = options.large;
        s = [s, sprintf(['\t' '\t' '\t{' vec2strjoin(ub(options.parametric.ub)) '}, /* ub */\n'])];
    end
    s = [s, sprintf(['\t' '};'])];

end