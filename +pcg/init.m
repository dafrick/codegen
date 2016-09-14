function init()

    %% Check whether preferences file is present
    

    %% Check if MPT is installed
    p = which('mpt_init');
    if isempty(p)
        error(sprintf(['The multi-parametric toolbox (MPT) could not be found. Please make sure it is set-up correcly and try again.\n' ...
               'Instructions on how to install it can be found on http://control.ee.ethz.ch/~mpt/3/Main/Installation']));
    end
    %% Check for lp solvers
    py = which('yalmip');
    pc = which('cplexlp');
    pg = which('gurobi');
    if isempty(py) && isempty(pc) && isempty(pg)
       error('No linear programming (LP) solvers found. Please make sure that you have either YALMIP (with appropriate solvers), CPLEX or Gurobi installed.');
    end
    if ~isempty(py) + ~isempty(pc) + ~isempty(pg) >= 2
        display('The following solvers were found:');
        choice = [];
        if ~isempty(py)
            display('YALMIP');
            choice = 'yalmip';
        end
        if ~isempty(pc)
            display('CPLEX');
            if isempty(choice)
                choice = 'CPLEX';
            else
                choice = [choice '|cplex'];
            end
        end
        if ~isempty(pg)
            display('Gurobi');
            choice = [choice '|gurobi'];
        end
        yn = input(['Choose one of the solvers to use [' choice '|default|abort]'], 's');
        while ~any(strcmp(yn,{'yalmip', 'cplex', 'gurobi', 'default', 'abort'}))
            yn = input(['Please type [' choice '|default|abort] in order to proceed.'], 's');
        end
        switch yn
            case 'yalmip',
            case 'cplex',
            case 'gurobi',
            case 'default',
            case 'abort',
                error('Aborting...');
        end
    else
        if ~isempty(py)
            display('Setting YALMIP as linear programming (LP) solver. If you want to use a different solver make sure it is installed and re-run pcg.init(''overwrite'', true).');
        end
        if ~isempty(pc)
            display('Setting CPLEX as linear programming (LP) solver. If you want to use a different solver make sure it is installed and re-run pcg.init(''overwrite'', true).');
        end
        if ~isempty(pg)
            display('Setting Gurobi as linear programming (LP) solver. If you want to use a different solver make sure it is installed and re-run pcg.init(''overwrite'', true).');
        end
    end


end
