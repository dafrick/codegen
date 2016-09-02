function [x, u, w] = recoverStatesInputs(z, dims, varargin)
    p = inputParser;
    p.addParameter('auxiliaries', false, @islogical);
    p.parse(varargin{:});
    options = p.Results;
    
    if ~options.auxiliaries
        z = reshape(z, dims.nu+dims.nx, []);
        u = z(1:dims.nu,:);
        x = z(dims.nu+1:end,:);
        w = [];
    else
        z = reshape(z, dims.nu+dims.nw+dims.nx, []);
        u = z(1:dims.nu,:);
        w = z(dims.nu+1:dims.nu+dims.nw,:);
        x = z(dims.nu+dims.nw+1:end,:);
    end

end

