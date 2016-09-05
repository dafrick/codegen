function [simple, info] = simplifyProjection(varargin)
    p = inputParser;
    p.addRequired('sol');
    p.addRequired('infbound', @isnumeric);
    p.addParameter('eps', 10,@(x)(isnumeric(x) && x > 0)); % Epsilon to shrink infinity bound
    p.addParameter('verbose', 2, @isnumeric);
    p.parse(varargin{:});
    options = p.Results;
    
    sol = options.sol;
    
    % get dimensions
    n = sol.Dim;
    
    % Setup the box constraints for the MPT3 infinity bound as Ax<=b
    box.A = [eye(n); -eye(n)];
    box.b = repmat(options.infbound,2*n,1) - options.eps;
    
    %% Iterate over all regions of solution
    if options.verbose >= 2
        display('Simplifying...');
    end
    simple  = [];
    nConstraints = 0;
    nPurged = 0;
    for i=1:length(sol)
        % Get inequality constraints of region
        A = sol(i).A; b = sol(i).b;
        nConstraints = nConstraints + size(A,1);
        isArtifact = false(size(A,1),1);
        % Iterate over all constraints
        for j=1:size(A,1)
            %[~,~,exitflag] = cplexlp(zeros(n,1), box.A, box.b, A(j,:), b(j,:));
            x = sdpvar(n,1);
            C = [box.A*x<=box.b, A(j,:)*x == b(j,:)];
            J = 0;
            result = optimize(C,J,sdpsettings('verbose',0));
            exitflag = result.problem;
            if exitflag == 1
                isArtifact(j) = true;
            end
        end
        % If some constraints of one region aren't artifacts
        if any(isArtifact)
            % Purge the offending constraints
            idx = find(~isArtifact);
            nPurged = nPurged+size(A,1)-length(idx);
            As = A(idx,:); bs = b(idx);
            P = Polyhedron('A', As, 'b', bs, 'Ae', sol(i).Ae, 'be', sol(i).be);
            P.addFunction(sol(i).getFunction('primal'), 'primal'); % Attach solution of projection
            P.addFunction(sol(i).getFunction('obj'), 'obj'); % Attach objective value of projection
            % Add to simplified solution
            simple = [simple P];
        % Otherwise complete region is purged
        end
    end
    if options.verbose >= 2
        display(['...done. Removed ' num2str(length(sol)-length(simple)) '/'  num2str(length(sol)) ' regions and ' num2str(nPurged) '/' num2str(nConstraints) ' constraints.']);
    end
    
    info.constraints = [nConstraints-nPurged, nConstraints];
    info.regions = [length(simple), length(sol)];
    
end
