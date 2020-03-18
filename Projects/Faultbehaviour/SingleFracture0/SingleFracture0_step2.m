SingleFracture0_step1

%% Resumption
resumption.active      = 1;
resumption.step.active = 1;
resumption.step.number = 1;
saveas_step.number     = 2;
%% write boudary conditions

clear boundary;
boundary.Ux.count = 4;
boundary.Ux.text{1} = '3 , 2e2, 0.0, 0.0';
boundary.Ux.text{2} = '4 , 2e2, 0.0, 0.0';
boundary.Ux.text{3} = '7 , 0.0, 0.0, 0.0';
boundary.Ux.text{4} = '8 , 0.0, 0.0, 0.0';

boundary.Uy.count = 2;
boundary.Uy.text{1} = '7 , 0.0, 0.0, 0.0';
boundary.Uy.text{2} = '8 , 0.0, 0.0, 0.0';


writeBoundaryConditions( prob_info,boundary )


%% write calculation parameters


% Parameter file
Parameter.proj_path = prob_info.proj_path;
Parameter.proj_name = prob_info.proj_name;
% Problem type information
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
% Special parameters

Parameter.specpara.staging     = 0;
Parameter.specpara.stepnum     =  saveas_step.number; % savee all the results in this step
Parameter.specpara.boundaryforce = 0;

% Load parameters

Parameter.load.resumption.active = resumption.active; 
Parameter.load.maxratio          = 1.0;
Parameter.load.volumeforce.active= 0;
Parameter.load.volumeforce.gx    = 0;
Parameter.load.volumeforce.gy    = -0.0981;
Parameter.load.resumption.stepnum = resumption.step.number;    % 91 Step num for resumption
Parameter.load.resumption.stepactive =  resumption.step.active;       % 99

% Calculation parameters
Parameter.calpara.loadincrement         = 1000;
Parameter.calpara.itermax               = 1000;
Parameter.calpara.tolerance.criteria    = 1e-6;
Parameter.calpara.tolerance.convergence = 1e-6;
Parameter.calpara.tolerance.displacement= 1e-4;

Parameter.calpara.time.start            = 0;
Parameter.calpara.time.end              = 1;
Parameter.calpara.time.increment        = 1e-1;



%% run Disroc
runDisroc(Parameter,Material,Disroc_path)
%% save to step 1
foldername = [Parameter.proj_path,'\STEP-',num2str(saveas_step.number)];
if not(exist(foldername,'dir'))
        mkdir(foldername)
end
cmd_txt = ['copy', ' ', Parameter.proj_path,'\RepM.dat', ' ',...
     foldername,'\RepM.dat'];
system(cmd_txt)
cmd_txt = ['copy', ' ', Parameter.proj_path,'\',...
     erase(Parameter.proj_name,'gid'), '*', ' ',...
     foldername,'\',erase(Parameter.proj_name,'gid'), '*'];
 system(cmd_txt)
 
 %% Plot
 singleFracturePlot

