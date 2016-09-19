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

function generateProjections(varargin)
    p = inputParser;
    p.addRequired('model', @isstruct);
    p.addRequired('N', @(x)(isnumeric(x) && length(x) == 1 && x > 0 && mod(x,1) == 0));
    p.addParameter('gendir', './gen', @ischar);
    p.addParameter('verbose', 0, @isnumeric);
    p.addParameter('generateProjections', true, @islogical);
    p.addParameter('lpSolver', 'yalmip', @(s)(any(strcmp(s,{'yalmip', 'cplex', 'gurobi'})))); % LP solver to use in code generation
    p.addParameter('eps', 1e-4, @isnumeric);
    p.addParameter('ignoreConst', false, @islogical);
    p.parse(varargin{:});
    options = p.Results;
    model = options.model;
    
    % Consistency checks
    if model.dims.nx ~= model.dims.nw
        throw(MException('generateProjections:InvalidDimensions', ['The dimension of the states nx=' num2str(model.dims.nx) ' needs to match the dimension of the auxiliaries nw=' num2str(model.dims.nw) '.']));
    end

    if options.verbose >= 1
        display('Generating projection...');
    end
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

            % Standard stage projection
            %          [x]                 [x]                    [x]
            % [-A -B I][u] = [c], [Aineq 0][u] < = [bineq], lb <= [u] <= ub
            %          [w]                 [w]
            % Constuct constraint data
            C.Aeq = [-model.dyn.A{i} -model.dyn.B{i} eye(model.dims.nx)];
            C.beq = model.dyn.c{i};
            if ~isempty(Aeq{i})
                C.Aeq = [C.Aeq; Aeq{i} zeros(size(Aeq{i},1),model.dims.nw)];
                C.beq = [C.beq; beq{i}];
            end
            if ~isempty(Aineq{i})
                C.Aineq = [Aineq{i}, zeros(size(Aineq{i},1),model.dims.nw)];
                C.bineq = bineq{i};
            else
                C.Aineq = [];
                C.bineq = [];
            end
            if ~isempty(lb{i})
                C.lb = [lb{i}; -inf(model.dims.nx,1)];
            else
                C.lb = [];
            end
            if ~isempty(ub{i})
                C.ub = [ub{i}; inf(model.dims.nx,1)];
            else
                C.ub = [];
            end
            pcg.generateProjection(['stdproj' num2str(i)], model.dims.n, C(:), 'gendir', options.gendir, 'verbose', options.verbose, 'simplify', true, 'lpSolver', options.lpSolver);

            % Initial stage projection
            %          [theta]                   [theta]
            % [-A -B I][u]     = [c]  , [Aineq 0][u]     <= [bineq], lb <= u <= ub
            % [ Aeq  0][w]     = [beq]           [w]
            % Constuct constraint data
            C.Aeq = [-model.dyn.A{i} -model.dyn.B{i} eye(model.dims.nx)];
            C.beq = model.dyn.c{i};
            if ~isempty(Aeq{i})
                C.Aeq = [C.Aeq; Aeq{i} zeros(size(Aeq{i},1),model.dims.nw)];
                C.beq = [C.beq; beq{i}];
            end
            if ~isempty(Aineq{i})
                C.Aineq = [Aineq{i}, zeros(size(Aineq{i},1),model.dims.nw)];
                C.bineq = bineq{i};
            else
                C.Aineq = [];
                C.bineq = [];
            end
            if ~isempty(lb{i})
                C.lb = [-inf(model.dims.nx,1); lb{i}(model.dims.nx+1:model.dims.nx+model.dims.nu); -inf(model.dims.nx,1)];
            else
                C.lb = [];
            end
            if ~isempty(ub{i})
                C.ub = [inf(model.dims.nx,1); ub{i}(model.dims.nx+1:model.dims.nx+model.dims.nu); inf(model.dims.nx,1)];
            else
                C.ub = [];
            end
            pcg.generateProjection(['iniproj' num2str(i)], model.dims.nu+model.dims.nx, C(:), 'gendir', options.gendir, 'verbose', options.verbose, 'bintree', false, 'lpSolver', options.lpSolver);

            % Final stage projection
            %
            %
            % Constuct constraint data
            C = [];
            if ~isempty(lb{i})
                C.lb = lb{i}(1:model.dims.nx);
            else
                C.lb = [];
            end
            if ~isempty(ub{i})
                C.ub = ub{i}(1:model.dims.nx);
            else
                C.ub = [];
            end
            pcg.generateProjection(['finproj' num2str(i)], model.dims.nx, C(:), 'gendir', options.gendir, 'verbose', options.verbose, 'simplify', false, 'lpSolver', options.lpSolver);
        end
    end
    
    % Generate projection C-file
    f = fopen([options.gendir '/project.h'], 'w');
    fprintf(f, ['#ifndef _PROJECT_H_\n' ...
                '#define _PROJECT_H_\n\n']);
    for i=1:model.dims.nr
        fprintf(f, ['unsigned long stdproj' num2str(i) '(double *X, double *U);\n']);
        fprintf(f, ['unsigned long iniproj' num2str(i) '(double *X, double *U);\n']);
        fprintf(f, ['unsigned long finproj' num2str(i) '(double *X, double *U);\n']);
    end
    fprintf(f, '\n#endif /*_PROJECT_H_*/\n');
    fclose(f);
    f = fopen([options.gendir '/project.c'], 'w');
    fprintf(f, ['#include <math.h>\n' ...
                '#include "project.h"\n' ...
                '#include "mex.h"\n']);
    fprintf(f, ['\n/* Helper functions */\n' ...
                'void vcopy(double* v, unsigned int n, double* u) {\n' ...
                '\t' 'unsigned int i;\n' ...
                '\t' 'for(i=0; i<n; i++) {\n' ...
                '\t' '\t' 'u[i] = v[i];\n' ...
                '\t' '}\n' ...
                '}\n' ...
                'double vdist(double* v, double* u, unsigned int n) {\n' ...
                '\t' 'unsigned int i;\n' ...
                '\t' 'double dist = 0.0;\n' ...
                '\t' 'for(i=0; i<n; i++) {\n' ...
                '\t' '\t' 'dist += (v[i]-u[i])*(v[i]-u[i]);\n' ...
                '\t' '}\n' ...
                '\t' 'return dist;\n' ...
                '}\n']);
    fprintf(f, '\nvoid mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[] ) {\n');
    fprintf(f, ['\t' '/* Input processing, output preparation */\n' ...
                '\t' 'double* x = mxGetPr(prhs[0]);\n' ...
                '\t' 'double* theta = mxGetPr(prhs[1]);\n' ...
                '\t' 'plhs[0] = mxCreateDoubleMatrix(%i,1,mxREAL);\n' ...
                '\t' 'double* z = mxGetPr(plhs[0]);\n'], model.dims.nn);
    fprintf(f, ['\t' 'unsigned int k;\n']);
    fprintf(f, ['\t' 'double temp[%i];\n'], model.dims.n);
    fprintf(f, ['\t' 'double ini[%i];\n'], model.dims.n);
    fprintf(f, ['\t' 'double obj, best;\n']);
    fprintf(f, ['\n\t' '/* Initial state */\n']);
    fprintf(f, ['\t' 'vcopy(theta, %i, ini);\n'], model.dims.nx);
    fprintf(f, ['\t' 'vcopy(x, %i, &ini[%i]);\n'], model.dims.nu+model.dims.nw, model.dims.nx);
    fprintf(f, ['\t' 'best = 1.0e15f;\n']);
    for i=1:model.dims.nr
        fprintf(f, ['\t' '\t' '/* Region %i */\n'], i);
        fprintf(f, ['\t' '\t' '\t' 'if(iniproj%i(ini, temp) >= 1) { /* Project and check feasibility */\n'], i);
        fprintf(f, ['\t' '\t' '\t' '\t' 'obj = vdist(x, temp, %i);\n'], model.dims.nu+model.dims.nw);
        fprintf(f, ['\t' '\t' '\t' '\t' 'if(obj < best) {\n']);
        fprintf(f, ['\t' '\t' '\t' '\t' '\t' 'best = obj;\n']);
        fprintf(f, ['\t' '\t' '\t' '\t' '\t' 'vcopy(temp, %i, z);\n'], model.dims.nu+model.dims.nw);
        fprintf(f, ['\t' '\t' '\t' '\t' '}\n']);
        fprintf(f, ['\t' '\t' '\t' '}\n']);
    end
    fprintf(f, ['\n\t' '/* Iterate over k=1,...,N-1 */\n']);
    fprintf(f, ['\t' 'for(k=0; k<%i; k++) {\n'], options.N-1);
    for i=1:model.dims.nr
        fprintf(f, ['\t' '\t' '/* Region %i */\n'], i);
        if i==1
            fprintf(f, ['\t' '\t' '\t' 'stdproj%i(&x[%i+k*%i], &z[%i+k*%i]); /* Project */\n'], i, model.dims.nu+model.dims.nw, model.dims.n, model.dims.nu+model.dims.nw, model.dims.n);
            fprintf(f, ['\t' '\t' '\t' 'best = vdist(&x[%i+k*%i], &z[%i+k*%i], %i);\n'], model.dims.nu+model.dims.nw, model.dims.n, model.dims.nu+model.dims.nw, model.dims.n, model.dims.n);
        else
            fprintf(f, ['\t' '\t' '\t' 'stdproj%i(&x[%i+k*%i], temp); /* Project */\n'], i, model.dims.nu+model.dims.nw, model.dims.n);
            fprintf(f, ['\t' '\t' '\t' 'obj = vdist(&x[%i+k*%i], temp, %i);\n'], model.dims.nu+model.dims.nw, model.dims.n, model.dims.n);
            fprintf(f, ['\t' '\t' '\t' 'if(obj < best) {\n']);
            fprintf(f, ['\t' '\t' '\t' '\t' 'best = obj;\n']);
            fprintf(f, ['\t' '\t' '\t' '\t' 'vcopy(temp, %i, &z[%i+k*%i]);\n'], model.dims.n, model.dims.nu+model.dims.nw, model.dims.n);
            fprintf(f, ['\t' '\t' '\t' '}\n']);
        end
    end
    fprintf(f, ['\t' '}\n']);
    fprintf(f, ['\n\t' '/* Final state */\n']);
    for i=1:model.dims.nr
        fprintf(f, ['\t' '\t' '/* Region %i */\n'], i);
        if i==1
            fprintf(f, ['\t' '\t' '\t' 'finproj%i(&x[%i], &z[%i]); /* Project */\n'], i, model.dims.nn-model.dims.nx, model.dims.nn-model.dims.nx);
            fprintf(f, ['\t' '\t' '\t' 'best = vdist(&x[%i], &z[%i], %i);\n'], model.dims.nn-model.dims.nx, model.dims.nn-model.dims.nx, model.dims.nx);
        else
            fprintf(f, ['\t' '\t' '\t' 'finproj%i(&x[%i], temp); /* Project */\n'], i, model.dims.nn-model.dims.nx);
            fprintf(f, ['\t' '\t' '\t' 'obj = vdist(&x[%i], temp, %i);\n'], model.dims.nn-model.dims.nx, model.dims.nx);
            fprintf(f, ['\t' '\t' '\t' 'if(obj < best) {\n']);
            fprintf(f, ['\t' '\t' '\t' '\t' 'best = obj;\n']);
            fprintf(f, ['\t' '\t' '\t' '\t' 'vcopy(temp, %i, &z[%i]);\n'], model.dims.nx, model.dims.nn-model.dims.nx);
            fprintf(f, ['\t' '\t' '\t' '}\n']);
        end
    end
    fprintf(f, '}\n');
    fclose(f);
    if options.verbose >= 1
        display('C-code, done. Compiling...');
    end
    % Compile
    compile = ['mex ' options.gendir '/project.c'];
    for i=1:model.dims.nr
        compile = [compile ' ' options.gendir '/stdproj' num2str(i) '.c'];
        compile = [compile ' ' options.gendir '/iniproj' num2str(i) '.c'];
        compile = [compile ' ' options.gendir '/finproj' num2str(i) '.c'];
    end
    compile = [compile ' -outdir ' options.gendir];
    eval(compile);
    if options.verbose >= 1
        display('...done.');
    end
end

