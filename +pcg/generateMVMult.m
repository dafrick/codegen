% GENERATEMVMULT Generate a fixed-point iteration .mex file
% 
% generateProxCode(model, 'par1', val1, 'par2', val2, ...)
% 

function [data, code] = generateMVMult(varargin)
    p = inputParser;
    p.addRequired('M', @isnumeric);
    p.addRequired('Mname', @ischar)
    p.addRequired('vname', @ischar);
    p.addRequired('rname', @ischar);
    p.addParameter('prefix', @ischar);
    p.addParameter('sym', false, @islogical);
    p.addParameter('nDataTabs', 0);
    p.addParameter('nCodeTabs', 0);
    p.parse(varargin{:});
    options = p.Results;
    M = options.M;
    
    % Ensure symmetry
    if options.sym
        M = (M+M')/2;
    end
    
    if options.nDataTabs > 0
        dataTabs = repmat('\t',1,options.nDataTabs);
    else
        dataTabs = '';
    end
    if options.nCodeTabs > 0
        codeTabs = repmat('\t',1,options.nCodeTabs);
    else
        codeTabs = '';
    end
    data = [dataTabs 'static const double ' options.Mname '[] = {\n'];
    dataIdx = 0;
    code = [];
    firstRow = true;
    for i=1:size(M,1)
        firstCol = true;
        codeline = [codeTabs options.rname '[' num2str(i-1) '] = ' options.rname '[' num2str(i-1) ']'];
        for j=1:size(M,2)
            if M(i,j) ~= 0
                % Write matrix entry into data
                if firstCol && firstRow
                    data = [data dataTabs '\t' sprintf('%0.17g', M(i,j))];
                elseif firstCol
                    data = [data ',\n' dataTabs '\t' sprintf('%0.17g', M(i,j))];
                else
                    data = [data dataTabs ', ' sprintf('%0.17g', M(i,j))];
                end
                dataIdx = dataIdx+1;
                % Generate code
                codeline = [codeline ' + ' options.Mname '[' num2str(dataIdx-1) ']*' options.vname '[' num2str(j-1) ']'];
                % Update firstCol and firstRow
                if firstCol && firstRow
                    firstCol = false;
                    firstRow = false;
                elseif firstCol
                    firstCol = false;
                end
            end
        end
        if ~firstCol % Add line if there was something to add
            code = [code codeline ';\n'];
        end
    end
    data = [data ' };'];

end

