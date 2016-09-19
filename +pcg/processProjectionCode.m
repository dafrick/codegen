function processProjectionCode(dir, filename)

    % Open file and get code in large string
    f = fopen([dir '/' filename '.c'], 'rt');
    code = fread(f, inf, 'char=>char');
    fclose(f);
    % Modify
    code = strrep(code', '#define MPT_NR', '#include "project.h"\n\n#define MPT_NR');
    code = strrep(code, ['static unsigned long ' filename], ['unsigned long ' filename]);
    code = strrep(code, ['static long ' filename], ['unsigned long ' filename]);
    code = strrep(code, 'static double MPT_ST', 'static const double MPT_ST');
    code = strrep(code, 'static double MPT_F', 'static const double MPT_F');
    code = strrep(code, 'static double MPT_G', 'static const double MPT_G');
    f = fopen([dir '/' filename '.c'], 'wt');
    fprintf(f, code);
    fclose(f);

end