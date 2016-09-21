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

function generateFORCESProjections(varargin)
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
    end
    info = pcg.generateFORCESProjection(options.name, C(:));
    
    %j=6; x=randn(length(options.lb{j}),1); lb = [options.lb{j}; -inf(stages(1).dims.n-length(options.lb{j}),1)]; lb(isinf(lb)) = -1e3; ub = [options.ub{j}; inf(stages(1).dims.n-length(options.lb{j}),1)]; ub(isinf(ub)) = 1e3; [z, f] = allinone(struct('x', [x; zeros(stages(1).dims.n-length(options.lb{j}),1)], 'Aeq', [options.Aeq{j} zeros(size(options.Aeq{j},1), stages(1).dims.n-length(options.lb{j})); zeros(stages(1).dims.r-size(options.Aeq{j},1),stages(1).dims.n)], 'beq', [options.beq{j}; zeros(stages(1).dims.r-size(options.Aeq{j},1),1)], 'lb', lb(stages(1).ineq.b.lbidx), 'ub', ub(stages(1).ineq.b.ubidx))); [y,~,f] = cplexqp(eye(length(options.lb{j})), x, options.Aineq{j}, options.bineq{j}, options.Aeq{j}, options.beq{j}, options.lb{j}, options.ub{j}); y-z.z(1:length(options.lb{j}))
    
    f = fopen([options.gendir '/project.h'], 'w');
    fprintf(f, ['#ifndef _PROJECT_H_\n' ...
                '#define _PROJECT_H_\n\n']);
    fprintf(f, '\n#endif /*_PROJECT_H_*/\n');
    fclose(f);
    
    f = fopen([options.gendir '/project.c'], 'w');
    fprintf(f, ['#include <math.h>\n' ...
                '#include "project.h"\n' ...
                '#include "%s.h"\n' ...
                '#include "mex.h"\n'], options.name);
    fprintf(f, '\n/* Static data */\n');
    fprintf(f, [options.name '_output out; /* Output struct */\n']);
    fprintf(f, [options.name '_info info; /* Solver info struct */\n']);
    for i=1:model.dims.nr
        fprintf(f, ['\t' '/* Region %i */\n'], i);
        projectionSuffix = {'i', 's', 'f'};
        projectionName = {'initial stage', 'standard stage', 'final stage'};
        for p=1:3 
            idx = (i-1)*3+p;
            fprintf(f, ['\t' options.name '_params params%i%s = { /* %s */\n'], i, projectionSuffix{p}, projectionName{p});
            % Projectant
            fprintf(f, ['\t' '\t' '\t{' repmat('0.0, ', 1, info.dims.n-1) '0.0}, /* x (projectant) */\n']);
            % Equality constraint
            if info.dims.r > 0
                Aeq = [C.Aeq{idx} zeros(size(C.Aeq{idx},1), info.dims.n-size(C.Aeq{idx},2)); zeros(info.dims.r-size(C.Aeq{idx},1),info.dims.n)];
                fprintf(f, ['\t' '\t' '\t{' sprintf('%0.17g, ', Aeq(:,1)) ' /* Aeq (column major) */\n']);
                for j=2:info.dims.n-1
                    fprintf(f, ['\t' '\t' '\t' ' ' sprintf('%0.17g, ', Aeq(:,j)) '\n']);
                end
                fprintf(f, ['\t' '\t' '\t' ' ' sprintf('%0.17g, ', Aeq(1:end-1,end)) '%0.17g},' '\n'], Aeq(end,end));
                beq = [C.beq{idx}; zeros(info.dims.r-length(C.beq{idx}),1)];
                fprintf(f, ['\t' '\t' '\t{' sprintf('%0.17g, ', beq(1:end-1)) '%0.17g}, /* beq */\n'], beq(end));
            end
            % Inequality constraint
            if info.dims.p > 0
                Aineq = [C.Aineq{idx} zeros(size(C.Aineq{idx},1), info.dims.n-size(C.Aineq{idx},2)); zeros(info.dims.p-size(C.Aineq{idx},1),info.dims.n)];
                fprintf(f, ['\t' '\t' '\t{' sprintf('%0.17g, ', Aineq(:,1)) ' /* Aineq (column major) */\n']);
                for j=2:info.dims.n-1
                    fprintf(f, ['\t' '\t' '\t' ' ' sprintf('%0.17g, ', Aineq(:,j)) '\n']);
                end
                fprintf(f, ['\t' '\t' '\t' ' ' sprintf('%0.17g, ', Aineq(1:end-1,end)) '%0.17g},' '\n'], Aineq(end,end));
                bineq = [C.bineq{idx}; zeros(info.dims.r-length(C.bineq{idx}),1)];
                fprintf(f, ['\t' '\t' '\t{' sprintf('%0.17g, ', bineq(1:end-1)) '%0.17g}, /* bineq */\n'], bineq(end));
            end
            % Lower bound
            if info.dims.l > 0
                lb = [C.lb{idx}; -inf(info.dims.n-length(C.lb{idx}),1)];
                lb(isinf(lb)) = options.large;
                fprintf(f, ['\t' '\t' '\t{' sprintf('%0.17g, ', lb(info.lbidx(1:end-1))) '%0.17g}, /* lb */\n'], lb(info.lbidx(end)));
            end
            % Upper bound
            if info.dims.l > 0
                ub = [C.ub{idx}; -inf(info.dims.n-length(C.ub{idx}),1)];
                ub(isinf(ub)) = options.large;
                fprintf(f, ['\t' '\t' '\t{' sprintf('%0.17g, ', ub(info.ubidx(1:end-1))) '%0.17g}, /* ub */\n'], ub(info.ubidx(end)));
            end
            fprintf(f, ['\t' '};\n']);
        end
        beq = C.beq{(i-1)*3+1};
        fprintf(f, ['\t' 'double beq%i[%i] = {' sprintf('%0.17g, ', beq(1:end-1)) '%0.17g}; /* unmodified beq */\n'], i, length(beq), beq(end));
    end
    fprintf(f, ['\n/* Helper functions */\n' ...
            'void vcopy(const double* v, unsigned int n, double* u) { /* u=v */\n' ...
            '\t' 'unsigned int i;\n' ...
            '\t' 'for(i=0; i<n; i++) {\n' ...
            '\t' '\t' 'u[i] = v[i];\n' ...
            '\t' '}\n' ...
            '}\n' ...
            'void vcopyminus(const double* v, unsigned int n, double* u) { /* u=-v */\n' ...
            '\t' 'unsigned int i;\n' ...
            '\t' 'for(i=0; i<n; i++) {\n' ...
            '\t' '\t' 'u[i] = -v[i];\n' ...
            '\t' '}\n' ...
            '}\n' ...
            'void vadd(const double* v, const double* w, unsigned int n, double* u) { /* u = v+w */\n' ...
            '\t' 'unsigned int i;\n' ...
            '\t' 'for(i=0; i<n; i++) {\n' ...
            '\t' '\t' 'u[i] = v[i] + w[i];\n' ...
            '\t' '}\n' ...
            '}\n' ...
            'double vdist(const double* v, const double* u, unsigned int n) { /* |v-u|_2^2 */ \n' ...
            '\t' 'unsigned int i;\n' ...
            '\t' 'double dist = 0.0;\n' ...
            '\t' 'for(i=0; i<n; i++) {\n' ...
            '\t' '\t' 'dist += (v[i]-u[i])*(v[i]-u[i]);\n' ...
            '\t' '}\n' ...
            '\t' 'return dist;\n' ...
            '}\n']);
    fprintf(f, ['\n' '/** Projection **/\n' ...
                'void project(const double* x, const double* theta, double* z) {\n' ...
                '\t' 'unsigned int k, i;\n' ...
                ]);
    for i=1:model.dims.nr
        fprintf(f, ['\n\t' '/* Region %i */\n'], i);
        fprintf(f, ['\t' '\t' '/* Initial stage */\n']);
        fprintf(f, ['\t' '\t' 'vadd(beq%i, theta, %i, params%ii.beq); /* Set initial state constraint */\n'], i, model.dims.nx, i);
        fprintf(f, ['\t' '\t' 'vcopyminus(theta, %i, params%ii.x); /* Set projectant */\n'], model.dims.nx, i);
        fprintf(f, ['\t' '\t' 'vcopyminus(&x[0], %i, &params%ii.x[%i]); /* Set projectant */\n'], model.dims.nu+model.dims.nx, i, model.dims.nx);
        if info.dims.n-model.dims.n > 0
            fprintf(f, ['\t' '\t' 'for(i=0; i<%i; i++) { params%ii.x[%i+i] = 0.0; }\n'], info.dims.n-model.dims.n, i, model.dims.n);
        end
        fprintf(f, ['\t' '\t' options.name '_solve(&params%ii, &out, &info, 0); /* Project */\n'], i);
        fprintf(f, ['\t' '\t' 'vcopy(&out.z[%i], %i, &z[%i]); /* Copy solution */\n'], model.dims.nx, model.dims.nu+model.dims.nx, (i-1)*model.dims.nn);
        fprintf(f, ['\t' '\t' '/* Standard stage */\n']);
        fprintf(f, ['\t' '\t' 'for(k=0; k<%i; k++) {\n'], options.N-1);
        fprintf(f, ['\t' '\t' '\t' 'vcopyminus(&x[%i+k*%i], %i, params%is.x); /* Set projectant */\n'], model.dims.nu+model.dims.nx, model.dims.n, model.dims.n, i);
        if info.dims.n-model.dims.n > 0
            fprintf(f, ['\t' '\t' '\t' 'for(i=0; i<%i; i++) { params%is.x[%i+i] = 0.0; }\n'], info.dims.n-model.dims.n, i, model.dims.n);
        end
        fprintf(f, ['\t' '\t' '\t' options.name '_solve(&params%is, &out, &info, 0); /* Project */\n'], i);
        fprintf(f, ['\t' '\t' '\t' 'vcopy(out.z, %i, &z[%i+k*%i]); /* Copy solution */\n'], model.dims.n, (i-1)*model.dims.nn+model.dims.nu+model.dims.nx, model.dims.n);
        fprintf(f, ['\t' '\t' '}\n']);
        fprintf(f, ['\t' '\t' '/* Final stage */\n']);
        fprintf(f, ['\t' '\t' 'vcopyminus(&x[%i], %i, params%if.x); /* Set projectant */\n'], model.dims.nn-model.dims.nx, model.dims.nx, i);
        if info.dims.n-model.dims.nx > 0
            fprintf(f, ['\t' '\t' 'for(i=0; i<%i; i++) { params%if.x[%i+i] = 0.0; }\n'], model.dims.n-model.dims.nx, i, model.dims.nx);
        end
        fprintf(f, ['\t' '\t' options.name '_solve(&params%if, &out, &info, 0); /* Project */\n'], i);
        fprintf(f, ['\t' '\t' 'vcopy(out.z, %i, &z[%i]); /* Copy solution */\n'], model.dims.nx, i*model.dims.nn-model.dims.nx);
    end
    fprintf(f, '}\n');
    fprintf(f, '\nvoid mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {\n');
    fprintf(f, ['\t' '/* Input processing, output preparation */\n' ...
                '\t' 'double* x = mxGetPr(prhs[0]); /* Projectant */\n' ...
                '\t' 'double* theta = mxGetPr(prhs[1]); /* Initial state */\n' ...
                '\t' 'plhs[0] = mxCreateDoubleMatrix(%i,%i,mxREAL); \n' ...
                '\t' 'double* z = mxGetPr(plhs[0]); /* Projection */\n'], model.dims.nn, model.dims.nr);
    fprintf(f, ['\t' 'project(x, theta, z);\n']);
    fclose(f);
    
    if options.verbose >= 1
        display('C-code, done. Compiling...');
    end
    % Compile
    compile = ['mex ' options.gendir '/project.c ' options.name '/obj/' options.name '.o -I./' options.name '/include'];
    compile = [compile ' -outdir ' options.gendir];
    eval(compile);
    if options.verbose >= 1
        display('...done.');
    end
    
end

