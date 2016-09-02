function [M, c, W, update] = getIterationData(H, h, N, dims, xi)

    % Build equality constraints
    Aeq = kron(eye(N), [zeros(dims.nx,dims.nu) eye(dims.nx) -eye(dims.nw)]);
    beq = zeros(N*dims.nx,1);
    v = Aeq\beq;

    K = [H Aeq'; Aeq zeros(size(Aeq,1),size(Aeq,1))]^(-1);
    R = K(1:dims.nn,1:dims.nn);
    % Checking whether xi satisfies xi > 1/[smallest positive eigenvalue of R]
    ximin = 1/min(subsref(eig(R), struct('type', '()', 'subs', {{eig(R)>0}})));
    if xi <= ximin
        warning('getIterationData:InvalidProximalScaling', ['The proximal scaling xi is not large enough, needs to be larger than ' num2str(ximin) '. Setting to xi ' num2str(ximin*2) ' and proceeding']);
        xi = 2*ximin;
    end
    M = xi*((xi*R-eye(size(R)))\R);
    c = (xi*R-eye(size(R)))\(R*(h+H*v)-v);
    % Compute weigthing matrix W
    [V,D] = eig(M);
    d = diag(D);
    [~,i] = sort(d);
    d = 2*d;
    d(i(1:size(Aeq,1))) = -1;
    d = 1./d;
    W = V*diag(d)*V';
    
    % Update function for changing h
    update = @(hn) (xi*R-eye(size(R)))\(R*(hn+H*v)-v);

end

