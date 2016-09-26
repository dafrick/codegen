function [M, c, W, update, xi] = getIterationData(H, h, N, dims, xi, varargin)
    p = inputParser;
    p.addParameter('Aeq', [], @isnumeric);
    p.addParameter('beq', [], @isnumeric);
    p.addParameter('eps', 100*eps, @isnumeric);
    p.parse(varargin{:});
    options = p.Results;
    
    if isempty(options.Aeq)
        % Build equality constraints
        Aeq = kron(eye(N), [zeros(dims.nx,dims.nu) eye(dims.nx) -eye(dims.nw)]);
        beq = zeros(N*dims.nx,1);
    else
        Aeq = options.Aeq;
        beq = options.beq;
    end
    m = sum(svd(Aeq) > options.eps);
    v = Aeq\beq;

    K = [H Aeq'; Aeq zeros(size(Aeq,1),size(Aeq,1))]^(-1);
    R = K(1:dims.nn,1:dims.nn);
    % Checking whether xi satisfies xi > 1/[smallest positive eigenvalue of R]
    ximin = 1/min(subsref(sort(eig(R),1,'descend'), struct('type', '()', 'subs', {{1:m}})));
    if xi <= ximin
        warning('getIterationData:InvalidProximalScaling', ['The proximal scaling xi is not large enough, needs to be larger than ' num2str(ximin) '. Setting to xi ' num2str(ximin*2) ' and proceeding']);
        xi = 2*ximin;
    end
    M = xi*((xi*R-eye(size(R)))\R);
    c = (xi*R-eye(size(R)))\(R*(h+H*v)-v);
    % Compute weigthing matrix W
    [V,D] = eig((M+M')/2); % Ensure M is symmetric (small asymmetry can lead to large errors!)
    d = diag(D);
    [~,i] = sort(real(d));
    d = 2*d;
    d(i(1:m)) = -1;
    d = 1./d;
    W = V*diag(d)*V';
    W = real(W);
    
    % Sparsify
    R(abs(R)< options.eps) = 0;
    R = sparse(R);
    H(abs(H)< options.eps) = 0;
    H = sparse(H);
    M(abs(M)< options.eps) = 0;
    M = sparse(M);
    W(abs(W)< options.eps) = 0;
    W = sparse(W);
    
    % Update function for changing h
    update = @(hn) (xi*R-eye(size(R)))\(R*(hn+H*v)-v);

end

