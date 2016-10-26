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

function info = generateProjectionCode(varargin)
    checkProjector = @(s)(any(strcmp(s, {'explicit', 'forces'})));
    checkProjectors = @(s)(isstruct(s) && ...
                           isfield(s, 'solver') && ( (~iscell(s.solver) && checkProjector(s.solver)) || (iscell(s.solver) && all(cellfun(checkProjector, s.solver))) ) && ...
                           isfield(s, 'name') && ( (~iscell(s.name) && ischar(s.name)) || (iscell(s.name) && all(cellfun(@ischar, s.name)))) && ...
                           (~isfield(s, 'indices') || ( (~iscell(s.indices) && isnumeric(s.indices)) || (iscell(s.indices) && all(cellfun(@isnumeric, s.indices))) )) && ...
                           (~isfield(s, 'separateInitial') || islogical(s.separateInitial)) && ...
                           xor(iscell(s.solver), iscell(s.name)) && (~iscell(s.solver) || length(s.solver) == length(s.name)) && ...
                           ( ~isfield(s, 'indices') || ( xor(iscell(s.solver),iscell(s.indices)) && (~iscell(s.solver) || length(s.solver) == length(s.indices))) ) && ...
                           ( ~isfield(s, 'separateInitial') || (~iscell(s.solver) || length(s.separateInitial) == 1 || length(s.solver) == length(s.separateInitial)) ) );

    p = inputParser;
    p.addRequired('sets', @(s)(isstruct(s) || (iscell(s) && all(cellfun(@isstruct, s)))));
    p.addRequired('dims', @isstruct);
    p.addRequired('N', @(x)(isnumeric(x) && length(x) == 1 && x > 0 && mod(x,1) == 0));
    p.addParameter('dir', './gen', @ischar);
    p.addParameter('verbose', 0, @isnumeric);
    p.addParameter('big', 1e3, @isnumeric);
    p.addParameter('type', 'union', @(s)(any(strcmp(s, {'union', 'all'}))));
    p.addParameter('name', 'project', @ischar);
    p.addParameter('nameUtils', 'utils', @ischar);
    p.addParameter('enableFlag', false, @islogical);
    p.addParameter('enableTimings', false, @islogical);
    p.addParameter('enableWarmstart', false, @islogical);
    p.addParameter('initialState', true, @islogical);
    p.addParameter('projectors', struct('solver', 'explicit', 'name', 'eproj'), checkProjectors);
    p.addParameter('')
    p.parse(varargin{:});
    options = p.Results;
    
    [st,~] = dbstack();
    display(st.name);
    
    %% Input processing
    dims = options.dims;
    
    % Simplify sets
    sets = options.sets;
    if ~iscell(sets)
        sets = {sets};
    end
    dims.ns = length(sets);
    
    % Simplify projectors
    projectors = options.projectors;
    if ~iscell(projectors.solver)
        projectors.solver = {projectors.solver};
        projectors.name = {projectors.name};
        if isfield(projectors, 'indices')
            projectors.indices = {projectors.indices};
        end
        if isfield(projectors, 'separateInitial') && length(projectors.solver) > 1 && length(projectors.separateInitial) == 1
            projectors.separateInitial = repmat(projectors.separateInitial,1,length(projectors.solver));
        end
    end
    dims.np = length(projectors.name);
    % TODO, ensure uniqueness of names
    
    % Default values for projector indices
    if ~isfield(projectors, 'indices')
        if dims.np == 1
            projectors.indices = 1:dims.ns;
        elseif dims.np == dims.ns;
            for i=1:dims.np
                projectors.indices{i} = i;
            end
        else
            throw(MException('pcg:generateProjectionCode:UnspecifiedProtectorIndices', 'The indices for the projectors is not specified and the number of sets does not match the number of projectors.'));
        end
    end
    % Default values for projector separateInitial
    if ~isfield(projectors, 'separateInitial')
        for i=1:dims.np
            switch projectors.solver{i}
                case 'explicit',
                    projectors.separateInitial(i) = true;
                case 'forces',
                    projectors.separateInitial(i) = false;
            end
        end
    end

    % Make sure name strings are consistent
    sanitize = @(s)(strrep(s, ' ',  '_'));
    options.name = sanitize(options.name);
    capsname = upper(options.name);
    options.nameHelpers = sanitize(options.nameHelpers);
    for i=1:dims.np
        projectors.name{i} = sanitize(projectors.name{i});
    end
    
    % Define used strings for consistency
    naming.initial = '_initial';
    naming.final = '_final';
    
    % Check missing features
    if strcmp(options.type, 'union')
        throw(MException('pcg:generateProjectionCode:MissingImplementation', 'Unions of sets are not yet implemented.'));
    end
    if options.enableTimings
        throw(MException('pcg:generateProjectionCode:MissingImplementation', 'Timings are not yet implemented.'));
    end
    if options.enableWarmstart
        throw(MException('pcg:generateProjectionCode:MissingImplementation', 'Warmstarting is not yet implemented.'));
    end
    for i=1:dims.np
        if strcmp(projectors.solver{i}, 'explicit')
            throw(MException('pcg:generateProjectionCode:MissingImplementation', 'Explicit projectors are not yet implemented.'));
        elseif strcmp(projectors.solver{i}, 'forces')
            if any(cellfun(@(s)(isfield(s,'end')), sets(projectors.indices{i})))
                throw(MException('pcg:generateProjectionCode:MissingImplementation', 'Final sets are not yet implemented for FORCES projectors.'));
            end
            if projectors.separateInitial(i)
                throw(MException('pcg:generateProjectionCode:MissingImplementation', 'Separating initial sets is not yet implemented for FORCES projectors.'));
            end
        end
    end
    
    %% Generate individual projectors/projection prototypes
    for i=1:dims.np
        switch projectors.solver{i}
            case 'explicit',
                % TODO
            case 'forces',
                if projectors.separateInitial(i)
                    % TODO
                else
                    projectors.info{i} = pcg.generateFORCESProjectionPrototype(options.nameFORCES, sets(projectors.indices{i}));
                end
        end
    end
    
    %% Generate library of helpers
    utilsInfo = pcg.generateUtils('name', options.nameUtils, 'overwrite', true);
    
    %% Generate projection H file
     f = fopen([options.dir '/' options.name '.h'], 'w');
    fprintf(f, ['#ifndef _%s_H_\n' ...
                '#define _%s_H_\n' ...
                '\n'], capsname, capsname);
    % TODO insert code documentation
    if options.enableFlag
        fprintf(f, 'void %s(const double* x, const double* theta, double* z, int* flag);\n', options.name);
    else
        fprintf(f, 'void %s(const double* x, const double* theta, double* z);\n', options.name);
    end
    fprintf(f, '\n#endif /*_%s_H_*/\n', capsname);
    fclose(f);
    
    %% Generate projection C file
    f = fopen([options.dir '/' options.name '.c'], 'w');
    fprintf(f, ['#include "%s.h"\n' ...
                '#include "%s.h"\n' ...
                '#include "%s.h"\n'], options.name, utilsInfo.name);
    for i=1:dims.np
        fprintf(f, '#include "%s.h"\n', projectors.info{i}.name);
    end
    % Static data
    fprintf(f, ['\n' ...
                '/*** Projectors static data ***/\n']);
    for i=1:dims.np
         switch projectors.solver
            case 'explicit',
                % TODO
            case 'forces',
                if length(projectors.indices{i}) == 1
                    fprintf(f, ['/* Projector %s, #%i interface (Set ' num2str(projectors.indices{i}) ') */\n'], projectors.name{i}, i);
                else
                    fprintf(f, ['/* Projector %s, #%i interface (Sets ' mat2str(projectors.indices{i}) ') */\n'], projectors.name{i}, i);
                end
                % Interface
                if projectors.separateInitial(i)
                    % TODO
                end
                fprintf(f, [projectors.info{i}.name '_output out; /* Output struct for FORCES projector %s */\n'], projectors.name{i});
                fprintf(f, [projectors.info{i}.name '_info info; /* Solver info struct for FORCES projector %s */\n'], projectors.name{i});
                if any(cellfun(@(s)(isfield(s,'end')), sets(projectors.indices{i})))
                    % TODO
                end
                % Sets
                for j=1:length(projectors.indices{i})
                    idx = projectors.indices{i}(j);
                    fprintf(f, '/* Set %i parameters (projector %s, #%i) */\n', idx, projectors.name{i}, i);
                    if projectors.separateInitial(i)
                        [s, name] = generateFORCESParameters(projectors.info{i}.initial.name, sets{idx}.init, projectors.info{i}.initial.parametric{idx}, ...
                                                        'suffix', [num2str(idx), naming.initial], ...
                                                        'comment', sprintf('Set %i (init) parameters (projector %s, #%i)', idx, projectors.name{i}, i), ...
                                                        'big', options.big);
                        projectors.info{i}.initial.parametersName{j} = name;
                        fprintf(f, '%s\n', s);
                    else
                        [s, name] = generateFORCESParameters(projectors.info{i}.name, sets{idx}.init, projectors.info{i}.parametric{idx}, ...
                                                            'suffix', [num2str(idx), naming.initial], ...
                                                            'comment', sprintf('Set %i (init) parameters (projector %s, #%i)', idx, projectors.name{i}, i), ...
                                                            'big', options.big);
                        projectors.info{i}.initial.parametersName{j} = name;
                        fprintf(f, '%s\n', s);
                    end
                    [s, name] = generateFORCESParameters(projectors.name{i}, sets{idx}, projectors.info{i}.parametric{idx}, ...
                                                        'suffix', num2str(idx), ...
                                                        'comment', sprintf('Set %i parameters (projector %s, #%i)', idx, projectors.name{i}, i), ...
                                                        'big', options.big);
                    projectors.info{i}.parametersName{j} = name;
                    fprintf(f, '%s\n', s);
                    if any(cellfun(@(s)(isfield(s,'end')), sets(projectors.indices{i})))
                        [s, name] = generateFORCESParameters(projectors.info{i}.name, sets{idx}.end, projectors.info{i}.parametric{idx}, ...
                                                            'suffix', [num2str(idx), naming.final], ...
                                                            'comment', sprintf('Set %i (final) parameters (projector %s, #%i)', idx, projectors.name{i}, i), ...
                                                            'big', options.big);
                        projectors.info{i}.initial.parametersName{j} = name;
                        fprintf(f, '%s\n', s);
                    end
                end
         end
    end
    fprintf(f, '/* DATA MARKER */\n');
    % Projection function
    fprintf(f, ['\n' ...
                '/*** Projection ***/\n']);
    if options.enableFlag
        fprintf(f, 'void %s(const double* x, const double* theta, double* z, int* flag) {\n', options.name);
    else
        fprintf(f, 'void %s(const double* x, const double* theta, double* z) {\n', options.name);
    end
    
    switch options.type
        case 'union', 
            % TODO
        case 'all',
            fprintf(f, ['\t' 'unsigned int k, i;\n']);
            for i=1:dims.nr
                fprintf(f, ['\n\t' '/* Region %i */\n'], i);
                fprintf(f, ['\t' '\t' '/* Initial stage */\n']);
                if options.initialState
                    fprintf(f, ['\t' '\t' 'vadd(beq%i, theta, %i, params%ii.beq); /* Set initial state constraint */\n'], i, dims.nx, i);
                    fprintf(f, ['\t' '\t' 'vcopyminus(theta, %i, params%ii.x); /* Set projectant */\n'], dims.nx, i);
                    fprintf(f, ['\t' '\t' 'vcopyminus(&x[0], %i, &params%ii.x[%i]); /* Set projectant */\n'], dims.nu+dims.nx, i, dims.nx);
                else
                    fprintf(f, ['\t' '\t' 'vcopyminus(&x[0], %i, &params%ii.x[0]); /* Set projectant */\n'], dims.nu+dims.nx, i);
                end
                if info.dims.n-dims.n > 0
                    fprintf(f, ['\t' '\t' 'for(i=0; i<%i; i++) { params%ii.x[%i+i] = 0.0; }\n'], info.dims.n-dims.n, i, dims.n);
                end
                if options.returnFlag
                    fprintf(f, ['\t' '\t' 'flag[(%i-1)*(%i+1)] = ' options.nameFORCES '_solve(&params%ii, &out, &info, 0); /* Project */\n'], i, options.N, i);
                else
                    fprintf(f, ['\t' '\t' options.nameFORCES '_solve(&params%ii, &out, &info, 0); /* Project */\n'], i);
                end
                if options.initialState
                    fprintf(f, ['\t' '\t' 'vcopy(&out.z[%i], %i, &z[%i]); /* Copy solution */\n'], dims.nx, dims.nu+dims.nx, (i-1)*dims.nn);
                else
                    fprintf(f, ['\t' '\t' 'vcopy(&out.z[0], %i, &z[%i]); /* Copy solution */\n'], dims.nu+dims.nx, (i-1)*dims.nn);
                end
                fprintf(f, ['\t' '\t' '/* Standard stage */\n']);
                fprintf(f, ['\t' '\t' 'for(k=0; k<%i; k++) {\n'], options.N-1);
                fprintf(f, ['\t' '\t' '\t' 'vcopyminus(&x[%i+k*%i], %i, params%is.x); /* Set projectant */\n'], dims.nu+dims.nx, dims.n, dims.n, i);
                if info.dims.n-dims.n > 0
                    fprintf(f, ['\t' '\t' '\t' 'for(i=0; i<%i; i++) { params%is.x[%i+i] = 0.0; }\n'], info.dims.n-dims.n, i, dims.n);
                end
                if options.returnFlag
                    fprintf(f, ['\t' '\t' '\t' 'flag[(%i-1)*(%i+1)+k+1] = ' options.nameFORCES '_solve(&params%is, &out, &info, 0); /* Project */\n'], i, options.N, i);
                else
                    fprintf(f, ['\t' '\t' '\t' options.nameFORCES '_solve(&params%is, &out, &info, 0); /* Project */\n'], i);
                end
                fprintf(f, ['\t' '\t' '\t' 'vcopy(out.z, %i, &z[%i+k*%i]); /* Copy solution */\n'], dims.n, (i-1)*dims.nn+dims.nu+dims.nx, dims.n);
                fprintf(f, ['\t' '\t' '}\n']);
                fprintf(f, ['\t' '\t' '/* Final stage */\n']);
                fprintf(f, ['\t' '\t' 'vcopyminus(&x[%i], %i, params%if.x); /* Set projectant */\n'], dims.nn-dims.nx, dims.nx, i);
                if info.dims.n-dims.nx > 0
                    fprintf(f, ['\t' '\t' 'for(i=0; i<%i; i++) { params%if.x[%i+i] = 0.0; }\n'], dims.n-dims.nx, i, dims.nx);
                end
                if options.returnFlag
                    fprintf(f, ['\t' '\t' 'flag[(%i-1)*(%i+1)+%i] = ' options.nameFORCES '_solve(&params%if, &out, &info, 0); /* Project */\n'], i, options.N, options.N, i);
                else
                    fprintf(f, ['\t' '\t' options.nameFORCES '_solve(&params%if, &out, &info, 0); /* Project */\n'], i);
                end
                fprintf(f, ['\t' '\t' 'vcopy(out.z, %i, &z[%i]); /* Copy solution */\n'], dims.nx, i*dims.nn-dims.nx);
            end
    end
    fprintf(f, '}\n');
    fclose(f);
    
    %% Generate projection MEX file
    f = fopen([options.dir '/' options.name '_mex.c'], 'w');
    fprintf(f, ['#include "%s.h"\n' ...
                '#include "mex.h"\n\n'], options.name);
    fprintf(f, '\nvoid mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {\n');
    fprintf(f, ['\t' '/* Input processing, output preparation */\n' ...
                '\t' 'double* x = mxGetPr(prhs[0]); /* Projectant */\n' ...
                '\t' 'double* theta = mxGetPr(prhs[1]); /* Initial state */\n' ...
                '\t' 'plhs[0] = mxCreateDoubleMatrix(%i,%i,mxREAL); \n' ...
                '\t' 'double* z = mxGetPr(plhs[0]); /* Projection */\n' ...
                '\t' 'plhs[1] = mxCreateNumericMatrix(%i,%i,mxINT32_CLASS,mxREAL);\n' ...
                '\t' 'int* feas = mxGetData(plhs[1]); /* Projection costs */\n'], dims.nn, dims.nr, options.N+1, dims.nr);
    fprintf(f, ['\t' '%s(x, theta, z, feas);\n'], options.name);
    fprintf(f, ['\t' 'if(nlhs > 2) {\n' ...
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
    % TODO: Properly handle relative and absolute paths in dir
    includes = {['./' options.dir], ['./' options.nameFORCES '/include']};
    objects = {[options.dir '/' options.name '.c'], [options.dir '/' options.nameHelpers '.c'], [options.nameFORCES '/obj/' options.nameFORCES '.o']};
    compile = ['mex ' options.dir '/' options.name '_mex.c'];
    compile = [compile ' -I' strjoin(includes, ' -I')];
    compile = [compile ' ' strjoin(objects, ' ')];
    compile = [compile ' -outdir ' options.dir ' -output ' options.name];
    eval(compile);
    if options.verbose >= 1
        display('...done.');
    end
    
    info.name = options.name;
    info.nameHelpers = options.nameHelpers;
    info.nameFORCES = options.nameFORCES;
    info.compile.includes = includes;
    info.compile.objects = objects;
    
end

