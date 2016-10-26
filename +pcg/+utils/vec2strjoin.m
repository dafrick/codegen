function s = vec2strjoin(varargin)
    p = inputParser;
    p.addRequired('v', @isnumeric);
    p.addRequired('sep', @ischar);
    p.addParameter('format', '%i', @ischar);
    p.parse(varargin{:});
    options = p.Results;

    s = strjoin(cellfun(@(s)(sprintf(options.format, s)), num2cell(options.v), 'UniformOutput', false), options.sep);

end
