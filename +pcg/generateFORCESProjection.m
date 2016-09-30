% GENERATEPROJECTIONS Generate a projection .mex file using MPT3
% 
% generateProjections(model, 'par1', val1, 'par2', val2, ...)
% 
% - model needs to contain structs
%   - dims, with the problem dimensions
%     - nx, number of states
%     - nu, number of inputs
%     - nw, number of auxiliaries - needs to be equal to nx
%     - nr, number of regions of the PWA dyamics
%     - n, number of decision variables per stage (= nx+nu+nw),
%     - nn, number of decision variables (= n*N)
%   - dyn, with the description of the model dynamics
%     - lb, a cell of length nr containing the region lower bounds [xmin; umin]
%     - ub, a cell of length nr containing the region upper bounds [xmax; umax]
%     - Aineq, a cell of length nr with a matrix of dimension (m)x(nx+nu) and
%     - bineq, a cell of length nr with a vector of length m containing the
%        region inequality constraints Aineq*[x; u] <= bineq
%     - Aeq, a cell of length nr with a matrix of dimension (m)x(nx+nu) and
%     - beq, a cell of length nr with a vector of length m containing the
%        region equality constraints Aeq*[x; u] = beq
%     - A, a cell of length nr with a matrix of dimension (nx)x(nx) and
%     - B, a cell of length nr with a matrix of dimension (nx)x(nu) and
%     - c, a cell of length nr with a vector of length nx containing the
%        region dynamics x(k+1) = Ax(k) + Bu(k) +c
%   - const, with the description of the operating constraints
%     - lb, a vector of lower bounds [xmin; umin]
%     - ub, a vector of upper bounds [xmax; umax]
%     - Aineq, matrix of dimension (m)x(nx+nu) and
%     - bineq, a vector of length m containing the 
%        inequality constraints Aineq*[x; u] <= bineq
%     - Aeq, a matrix of dimension (m)x(nx+nu) and
%     - beq, a vector of length m containing the
%        equality constraints Aeq*[x; u] = beq

function info = generateFORCESProjection(varargin)
    p = inputParser;
    p.addRequired('model', @isstruct);
    p.addRequired('N', @(x)(isnumeric(x) && length(x) == 1 && x > 0 && mod(x,1) == 0));
    p.addParameter('gendir', './gen', @ischar);
    p.addParameter('verbose', 0, @isnumeric);
    p.addParameter('generateProjections', true, @islogical);
    p.addParameter('lpSolver', 'yalmip', @(s)(any(strcmp(s,{'yalmip', 'cplex', 'gurobi'})))); % LP solver to use in code generation
    p.addParameter('eps', 1e-4, @isnumeric);
    p.addParameter('large', 1e3, @isnumeric);
    p.addParameter('ignoreConst', false, @islogical);
    p.addParameter('name', 'allinone', @ischar);
    p.parse(varargin{:});
    options = p.Results;
    model = options.model;
    
    % Consistency checks
    if model.dims.nx ~= model.dims.nw
        throw(MException('generateFORCESProjections:InvalidDimensions', ['The dimension of the states nx=' num2str(model.dims.nx) ' needs to match the dimension of the auxiliaries nw=' num2str(model.dims.nw) '.']));
    end

    if options.verbose >= 1
        display('Generating projection...');
    end
    
    j = 1;
    if options.generateProjections
        for i=1:model.dims.nr
            % Combine lower bounds
            if isfield(model.dyn, 'lb') && ~isempty(model.dyn.lb) && length(model.dyn.lb{i}) == model.dims.nx+model.dims.nu
                if (isfield(model.const, 'lb') && ~isempty(model.const.lb) && length(model.const.lb) == model.dims.nx+model.dims.nu)
                    lb{i} = max(model.dyn.lb{i}, model.const.lb);
                end
            elseif (isfield(model.const, 'lb') && ~isempty(model.const.lb) && length(model.const.lb) == model.dims.nx+model.dims.nu)
                lb{i} = model.const.lb;
            end
            % Combine upper bounds
            if isfield(model.dyn, 'ub') && ~isempty(model.dyn.ub) && length(model.dyn.ub{i}) == model.dims.nx+model.dims.nu
                if (isfield(model.const, 'ub') && ~isempty(model.const.ub) && length(model.const.ub) == model.dims.nx+model.dims.nu)
                    ub{i} = min(model.dyn.ub{i}, model.const.ub);
                end
            elseif (isfield(model.const, 'ub') && ~isempty(model.const.ub) && length(model.const.ub) == model.dims.nx+model.dims.nu)
                ub{i} = [];
            end
            % Combine equality constraints
            Aeq{i} = [];
            beq{i} = [];
            if isfield(model.dyn, 'Aeq') && ~isempty(model.dyn.Aeq) && ~isempty(model.dyn.Aeq{i})
                Aeq{i} = model.dyn.Aeq{i};
                beq{i} = model.dyn.beq{i};
            end
            if isfield(model.const, 'Aeq') && ~isempty(model.const.Aeq)
                Aeq{i} = [Aeq{i}; model.const.Aeq];
                beq{i} = [beq{i}; model.const.beq];
            end
            % Combine inequality constraints
            Aineq{i} = [];
            bineq{i} = [];
            if isfield(model.dyn, 'Aineq') && ~isempty(model.dyn.Aineq) && ~isempty(model.dyn.Aineq{i})
                Aineq{i} = model.dyn.Aineq{i};
                bineq{i} = model.dyn.bineq{i};
            end
            if isfield(model.const, 'Aineq') && ~isempty(model.const.Aineq)
                Aineq{i} = [Aineq{i}; model.const.Aineq];
                beq{i} = [bineq{i}; model.const.bineq];
            end
            
            % Initial stage projection
            %          [theta]                   [theta]
            % [-A -B I][u]     = [c]  , [Aineq 0][u]     <= [bineq], lb <= u <= ub
            % [ Aeq  0][w]     = [beq]           [w]
            % Constuct constraint data
            C.Aeq{j} = [eye(model.dims.nx) zeros(model.dims.nx, model.dims.nu+model.dims.nx); -model.dyn.A{i} -model.dyn.B{i} eye(model.dims.nx)];
            C.beq{j} = [zeros(model.dims.nx,1); model.dyn.c{i}];
            if ~isempty(Aeq{i})
                C.Aeq{j} = [C.Aeq; Aeq{i} zeros(size(Aeq{i},1),model.dims.nw)];
                C.beq{j} = [C.beq; beq{i}];
            end
            if ~isempty(Aineq{i})
                C.Aineq{j} = [Aineq{i}, zeros(size(Aineq{i},1),model.dims.nw)];
                C.bineq{j} = bineq{i};
            else
                C.Aineq{j} = [];
                C.bineq{j} = [];
            end
            if ~isempty(lb{i})
                C.lb{j} = [lb{i}; -inf(model.dims.nx,1)];
            else
                C.lb{j} = [];
            end
            if ~isempty(ub{i})
                C.ub{j} = [ub{i}; inf(model.dims.nx,1)];
            else
                C.ub{j} = [];
            end
            j= j+1;
            % Standard stage projection
            %          [x]                 [x]                    [x]
            % [-A -B I][u] = [c], [Aineq 0][u] < = [bineq], lb <= [u] <= ub
            %          [w]                 [w]
            % Constuct constraint data
            C.Aeq{j} = [-model.dyn.A{i} -model.dyn.B{i} eye(model.dims.nx)];
            C.beq{j} = model.dyn.c{i};
            if ~isempty(Aeq{i})
                C.Aeq{j} = [C.Aeq; Aeq{i} zeros(size(Aeq{i},1),model.dims.nw)];
                C.beq{j} = [C.beq; beq{i}];
            end
            if ~isempty(Aineq{i})
                C.Aineq{j} = [Aineq{i}, zeros(size(Aineq{i},1),model.dims.nw)];
                C.bineq{j} = bineq{i};
            else
                C.Aineq{j} = [];
                C.bineq{j} = [];
            end
            if ~isempty(lb{i})
                C.lb{j} = [lb{i}; -inf(model.dims.nx,1)];
            else
                C.lb{j} = [];
            end
            if ~isempty(ub{i})
                C.ub{j} = [ub{i}; inf(model.dims.nx,1)];
            else
                C.ub{j} = [];
            end
            j=j+1;
            % Final stage projection
            %
            %
            % Constuct constraint data
            C.Aeq{j} = [];
            C.beq{j} = [];
            C.Aineq{j} = [];
            C.bineq{j} = [];
            if ~isempty(lb{i})
                C.lb{j} = lb{i}(1:model.dims.nx);
            else
                C.lb{j} = [];
            end
            if ~isempty(ub{i})
                C.ub{j} = ub{i}(1:model.dims.nx);
            else
                C.ub{j} = [];
            end
            j=j+1;
        end
        
        info = pcg.generateFORCESProjectionCode(C, model.dims, options.N, 'gendir', options.gendir, 'verbose', options.verbose);
    else
        info = [];
    end
    
end

