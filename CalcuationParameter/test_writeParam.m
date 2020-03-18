% % File operation
% fid = fopen('prj_name.txt','r');
% l1  = fscanf(fid,'%s');
% fclose(fid);


proj_name = 'SingleFracture0.gid';

% Parameter file
Parameter.proj_path = strcat(proj_path,'\',proj_name);
Parameter.proj_name = proj_name ;
%% Problem type information
Parameter.problem.physics      = 1; %1:Mechanics, 2:Hydraulic, 3:Thermal, 4:HM, 5:TM, 6:THM
Parameter.problem.time         = 1; %1:Time Independent, 2:Transient
Parameter.problem.type         = 1;
Parameter.problem.axesymmetry  = 0;
Parameter.problem.planetype    = 1;
Parameter.problem.generalized  = 0;  %?
Parameter.problem.hyro.matrix  = 1;
Parameter.problem.hydro.gravity.active = 0;
Parameter.problem.hydro.gravity.value  = 0.0098;
Parameter.problem.user         = 1;
%% Special parameters

 Parameter.specpara.staging     = 0;
 Parameter.specpara.stepnum     = 2; % savee all the results in this step
 Parameter.specpara.boundaryforce = 0;

%% Load parameters

Parameter.load.resumption.active = 1;
Parameter.load.maxratio          = 1.0;
Parameter.load.volumeforce.active= 0;
Parameter.load.volumeforce.gx    = 0;
Parameter.load.volumeforce.gy    = -0.0981;
Parameter.load.resumption.stepnum = 103;    % 91 Step num for resumption
Parameter.load.resumption.stepactive = 1;       % 99

%% Calculation parameters
Parameter.calpara.loadincrement         = 1000;
Parameter.calpara.itermax               = 1000;
Parameter.calpara.tolerance.criteria    = 1e-6;
Parameter.calpara.tolerance.convergence = 1e-6;
Parameter.calpara.tolerance.displacement= 1e-4;

Parameter.calpara.time.start            = 0;
Parameter.calpara.time.end              = 1;
Parameter.calpara.time.increment        = 1e-1;

writeParam(Parameter)

