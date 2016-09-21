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

function generateFORCESProjectionCode(varargin)
    p = inputParser;
    p.addRequired('C', @isstruct);
    p.addRequired('dims', @isstruct);
    p.addRequired('N', @(x)(isnumeric(x) && length(x) == 1 && x > 0 && mod(x,1) == 0));
    p.addParameter('gendir', './gen', @ischar);
    p.addParameter('verbose', 0, @isnumeric);
    p.addParameter('large', 1e3, @isnumeric);
    p.addParameter('name', 'project', @ischar);
    p.addParameter('nameFORCES', 'allinone', @ischar);
    p.addParameter('nameHelpers', 'helpers', @ischar);
    p.parse(varargin{:});
    options = p.Results;
    dims = options.dims;
    C = options.C;
    
    options.name = strrep(options.name, ' ', '_');
    options.nameFORCES = strrep(options.nameFORCES, ' ', '_');
    options.nameHelpers = strrep(options.nameHelpers, ' ', '_');
    
    info = pcg.generateFORCESProjectionPrototype(options.nameFORCES, C(:));
    
    %% Generate library of helpers
    f = fopen([options.gendir '/' options.nameHelpers '.h'], 'w');
    fprintf(f, ['#ifndef _%s_H_\n' ...
                '#define _%s_H_\n\n'], upper(options.nameHelpers), upper(options.nameHelpers));
    fprintf(f, ['void vcopy(const double* v, unsigned int n, double* u); /* u=v */\n' ...
                'void vcopyminus(const double* v, unsigned int n, double* u); /* u=-v */\n' ...
                'void vadd(const double* v, const double* w, unsigned int n, double* u); /* u = v+w */\n' ...
                'double vdist(const double* v, const double* u, unsigned int n); /* |v-u|_2^2 */ \n']);
	fprintf(f, '\n#endif /*_%s_H_*/\n', upper(options.nameHelpers));
    fclose(f);
    f = fopen([options.gendir '/' options.nameHelpers '.c'], 'w');
    fprintf(f, ['#include "' strrep(options.nameHelpers, ' ', '_') '.h"\n']);
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
    fclose(f);
    
    f = fopen([options.gendir '/' options.name '.h'], 'w');
    fprintf(f, ['#ifndef _%s_H_\n' ...
                '#define _%s_H_\n\n'], upper(options.name), upper(options.name));
    fprintf(f, 'void %s(const double* x, const double* theta, double* z);\n', options.name);
    fprintf(f, '\n#endif /*_%s_H_*/\n', upper(options.name));
    fclose(f);
    f = fopen([options.gendir '/' options.name '.c'], 'w');
    fprintf(f, ['#include "%s.h"\n' ...
                '#include "%s.h"\n' ...
                '#include "%s.h"\n'], options.name, options.nameHelpers, options.nameFORCES);
    fprintf(f, '\n/* Static data */\n');
    fprintf(f, [options.nameFORCES '_output out; /* Output struct */\n']);
    fprintf(f, [options.nameFORCES '_info info; /* Solver info struct */\n']);
    for i=1:dims.nr
        fprintf(f, ['\t' '/* Region %i */\n'], i);
        projectionSuffix = {'i', 's', 'f'};
        projectionName = {'initial stage', 'standard stage', 'final stage'};
        for p=1:3 
            idx = (i-1)*3+p;
            fprintf(f, ['\t' options.nameFORCES '_params params%i%s = { /* %s */\n'], i, projectionSuffix{p}, projectionName{p});
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
    fprintf(f, '/* DATA MARKER */\n');
    fprintf(f, ['\n' '/** Projection **/\n' ...
                'void project(const double* x, const double* theta, double* z) {\n' ...
                '\t' 'unsigned int k, i;\n' ...
                ]);
    for i=1:dims.nr
        fprintf(f, ['\n\t' '/* Region %i */\n'], i);
        fprintf(f, ['\t' '\t' '/* Initial stage */\n']);
        fprintf(f, ['\t' '\t' 'vadd(beq%i, theta, %i, params%ii.beq); /* Set initial state constraint */\n'], i, dims.nx, i);
        fprintf(f, ['\t' '\t' 'vcopyminus(theta, %i, params%ii.x); /* Set projectant */\n'], dims.nx, i);
        fprintf(f, ['\t' '\t' 'vcopyminus(&x[0], %i, &params%ii.x[%i]); /* Set projectant */\n'], dims.nu+dims.nx, i, dims.nx);
        if info.dims.n-dims.n > 0
            fprintf(f, ['\t' '\t' 'for(i=0; i<%i; i++) { params%ii.x[%i+i] = 0.0; }\n'], info.dims.n-dims.n, i, dims.n);
        end
        fprintf(f, ['\t' '\t' options.nameFORCES '_solve(&params%ii, &out, &info, 0); /* Project */\n'], i);
        fprintf(f, ['\t' '\t' 'vcopy(&out.z[%i], %i, &z[%i]); /* Copy solution */\n'], dims.nx, dims.nu+dims.nx, (i-1)*dims.nn);
        fprintf(f, ['\t' '\t' '/* Standard stage */\n']);
        fprintf(f, ['\t' '\t' 'for(k=0; k<%i; k++) {\n'], options.N-1);
        fprintf(f, ['\t' '\t' '\t' 'vcopyminus(&x[%i+k*%i], %i, params%is.x); /* Set projectant */\n'], dims.nu+dims.nx, dims.n, dims.n, i);
        if info.dims.n-dims.n > 0
            fprintf(f, ['\t' '\t' '\t' 'for(i=0; i<%i; i++) { params%is.x[%i+i] = 0.0; }\n'], info.dims.n-dims.n, i, dims.n);
        end
        fprintf(f, ['\t' '\t' '\t' options.nameFORCES '_solve(&params%is, &out, &info, 0); /* Project */\n'], i);
        fprintf(f, ['\t' '\t' '\t' 'vcopy(out.z, %i, &z[%i+k*%i]); /* Copy solution */\n'], dims.n, (i-1)*dims.nn+dims.nu+dims.nx, dims.n);
        fprintf(f, ['\t' '\t' '}\n']);
        fprintf(f, ['\t' '\t' '/* Final stage */\n']);
        fprintf(f, ['\t' '\t' 'vcopyminus(&x[%i], %i, params%if.x); /* Set projectant */\n'], dims.nn-dims.nx, dims.nx, i);
        if info.dims.n-dims.nx > 0
            fprintf(f, ['\t' '\t' 'for(i=0; i<%i; i++) { params%if.x[%i+i] = 0.0; }\n'], dims.n-dims.nx, i, dims.nx);
        end
        fprintf(f, ['\t' '\t' options.nameFORCES '_solve(&params%if, &out, &info, 0); /* Project */\n'], i);
        fprintf(f, ['\t' '\t' 'vcopy(out.z, %i, &z[%i]); /* Copy solution */\n'], dims.nx, i*dims.nn-dims.nx);
    end
    fprintf(f, '}\n');
    fclose(f);
    
    f = fopen([options.gendir '/' options.name '_mex.c'], 'w');
    fprintf(f, ['#include "%s.h"\n' ...
                '#include "mex.h"\n\n'], options.name);
    fprintf(f, '\nvoid mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {\n');
    fprintf(f, ['\t' '/* Input processing, output preparation */\n' ...
                '\t' 'double* x = mxGetPr(prhs[0]); /* Projectant */\n' ...
                '\t' 'double* theta = mxGetPr(prhs[1]); /* Initial state */\n' ...
                '\t' 'plhs[0] = mxCreateDoubleMatrix(%i,%i,mxREAL); \n' ...
                '\t' 'double* z = mxGetPr(plhs[0]); /* Projection */\n'], dims.nn, dims.nr);
    fprintf(f, ['\t' '%s(x, theta, z);\n'], options.name);
    fprintf(f, ['\t' 'if(nlhs > 1) {\n' ...
                '\t' '\t' 'plhs[1] = mxCreateDoubleMatrix(%i,%i,mxREAL);\n' ...
                '\t' '\t' 'double* C = mxGetPr(plhs[1]); /* Projection costs */\n' ...
                '\t' '\t' 'unsigned int i,k;\n'], options.N+1, dims.nr);
    fprintf(f, ['\t' '\t' 'for(i=0; i<%i; i++) {\n'], dims.nr);
    fprintf(f, ['\t' '\t' '\t' 'C[i*%i] = vdist(&x[0], &z[i*%i], %i); /* Initial stage */\n'], options.N+1, dims.nn, dims.nx+dims.nu);
    fprintf(f, ['\t' '\t' '\t' 'for(k=0; k<%i; k++) { C[i*%i+k+1] = vdist(&x[%i+k*%i], &z[i*%i+%i+k*%i], %i); } /* Standard stage */\n'], options.N-1, options.N+1, dims.nx+dims.nu, dims.n, dims.nn, dims.nx+dims.nu, dims.n, dims.n);
    fprintf(f, ['\t' '\t' '\t', 'C[i*%i+%i] = vdist(&x[%i], &z[%i*i+%i], %i); /* Final stage */\n'], options.N+1, options.N, dims.nn-dims.nx, dims.nn, dims.nn-dims.nx, dims.nx);
    fprintf(f, ['\t' '\t' '}\n' ...
                '\t' '}\n']);
    fprintf(f, '}\n');
    fclose(f);
    
    if options.verbose >= 1
        display('C-code, done. Compiling...');
    end
    % Compile
    % TODO: Properly handel relative and absolute paths in gendir
    compile = ['mex ' options.gendir '/' options.name '_mex.c'];
    compile = [compile ' -I./' options.gendir ...
                       ' -I./' options.nameFORCES '/include'];
    compile = [compile ' ' options.gendir '/' options.name '.c' ... 
                       ' ' options.gendir '/' options.nameHelpers '.c' ...
                       ' ' options.nameFORCES '/obj/' options.nameFORCES '.o'];
    compile = [compile ' -outdir ' options.gendir ' -output ' options.name];
    eval(compile);
    if options.verbose >= 1
        display('...done.');
    end
    
end

