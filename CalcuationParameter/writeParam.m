function writeParam(Parameter)

fname = strcat(Parameter.proj_path,'\', 'Param.dat');

fid = fopen( fname,'w+');

fprintf(fid, '%s\n','DISROC-5-8 200118');
fprintf(fid, '%s\n','ProblemType');
%% Problem type information
fprintf(fid, '%d\n', Parameter.problem.physics);
fprintf(fid, '%d\n', Parameter.problem.time);
fprintf(fid, '%d\n', Parameter.problem.type);
fprintf(fid, '%f\n', Parameter.problem.axesymmetry);
fprintf(fid, '%d\n', Parameter.problem.planetype);
fprintf(fid, '%f\n', Parameter.problem.generalized);
fprintf(fid, '%d\n', Parameter.problem.hyro.matrix);
fprintf(fid, '%d\n', Parameter.problem.hydro.gravity.active);
fprintf(fid, '%d\n', Parameter.problem.hydro.gravity.value);
fprintf(fid, '%d\n', Parameter.problem.user);

fprintf(fid, '%s\n','EndSection');
%% Special parameter
fprintf(fid, '%s\n','SpecialParameters');   % ????
fprintf(fid, '%d\n', Parameter.specpara.staging);
fprintf(fid, '%d\n', Parameter.specpara.stepnum);
fprintf(fid, '%d\n', Parameter.specpara.boundaryforce);
fprintf(fid, '%s\n','EndSection');

%% Load parameters
fprintf(fid, '%s\n','LoadParameters');
fprintf(fid, '%d\n',Parameter.load.resumption.active);
fprintf(fid, '%f\n',Parameter.load.maxratio);
fprintf(fid, '%d\n',Parameter.load.volumeforce.active);
fprintf(fid, '%f\n',Parameter.load.volumeforce.gx);
fprintf(fid, '%f\n',Parameter.load.volumeforce.gy);
fprintf(fid, '%d\n', 0);  % 39
fprintf(fid, '%d\n', 1);  % 40
fprintf(fid, '%d\n', 1);  % 15
fprintf(fid, '%d\n', 0);  % 51
fprintf(fid, '%d\n', 0);  % 60
fprintf(fid, '%d\n', 0);  % 52 fluid pressure forces ??
fprintf(fid, '%d\n', 0);  % 70 thermal stress forces ??
fprintf(fid, '%d\n', 0);  % 86
fprintf(fid, '%d\n', 0);  % 89
fprintf(fid, '%f\n', 1.5);% 23
fprintf(fid, '%f\n', 3.5);% 24
fprintf(fid, '%d\n', 4);  % 90
fprintf(fid, '%d\n', Parameter.load.resumption.stepnum);    % 91
fprintf(fid, '%d\n', Parameter.load.resumption.stepactive); % 99
fprintf(fid, '%d\n', 0);  % 24
fprintf(fid, '%d\n', 0);  % 67
fprintf(fid, '%d\n', 0);  % 101
fprintf(fid, '%d\n', 1);  % 103
fprintf(fid, '%s\n','EndSection');

%% Calculation parameters
fprintf(fid, '%s\n','CalculationParameters');
fprintf(fid, '%d\n', Parameter.calpara.loadincrement);
fprintf(fid, '%d\n', 1);
fprintf(fid, '%d\n', Parameter.calpara.itermax);
fprintf(fid, '%f\n', Parameter.calpara.tolerance.criteria);
fprintf(fid, '%d\n', 1);  % conjugate gradient
fprintf(fid, '%d\n', 10000);
fprintf(fid, '%f\n', Parameter.calpara.tolerance.convergence);
fprintf(fid, '%f\n', Parameter.calpara.tolerance.displacement);
fprintf(fid, '%f\n', 1.e-200 );
fprintf(fid, '%f\n', Parameter.calpara.time.start);
fprintf(fid, '%f\n', Parameter.calpara.time.end);
fprintf(fid, '%f\n', Parameter.calpara.time.increment);
fprintf(fid, '%d\n', 0);
fprintf(fid, '%d\n', 0);
fprintf(fid, '%s\n','EndSection');
%% Output parameters
fprintf(fid, '%s\n','OutPutParameters');

fprintf(fid, '%d\n', 0);   % 22
fprintf(fid, '%d\n', 1);   % 35
fprintf(fid, '%d\n', 1);   % 36
fprintf(fid, '%d\n', 0);   % 20
fprintf(fid, '%d\n', 1);   % 71
fprintf(fid, '%d\n', 1);   % 43
fprintf(fid, '%d\n', 0);   % 201
fprintf(fid, '%d\n', 0);   % 202 
fprintf(fid, '%d\n', 0);   % 203
fprintf(fid, '%d\n', 1);   % 204
fprintf(fid, '%d\n', 0);   % 205
fprintf(fid, '%d\n', 0);   % 206
fprintf(fid, '%d\n', 1);   % 207
fprintf(fid, '%d\n', 1);   % 208
fprintf(fid, '%d\n', 0);   % 209
fprintf(fid, '%d\n', 0);   % 92
fprintf(fid, '%d\n', 0);   % 210
fprintf(fid, '%d\n', 0);   % 211
fprintf(fid, '%d\n', 0);   % 212
fprintf(fid, '%d\n', 0);   % 213
fprintf(fid, '%d\n', 0);   % 26
fprintf(fid, '%d\n', 0);   % 27
fprintf(fid, '%d\n', 0);   % 214
fprintf(fid, '%s\n','EndSection');
fprintf(fid, '%s\n','EndFile');

fclose(fid);






