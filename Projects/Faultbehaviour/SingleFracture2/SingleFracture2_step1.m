clear;
Path_Control
%%
Kt   = 2e8;      % Kt 
Kn   = 5e8;      % Kn
e   = 10;       % e
sigma_R   = 1.2e6;    % sigma_R
tau_c     = 5.8e6;

phi   = 31;       % phi
Cohesion   = sqrt(2*sigma_R*tau_c*tan(phi/180*pi) -(sigma_R*tan(phi/180*pi))^2);    % C
hr   = 0.33;     % h
beta   = 0.1;      % beta
beta_angle   = 1;        % beta'
k0t  = 1e8;      % k0t
k0n  = 1e10;      % k0n

Sn_load  = -30e6;

resumption.active = 0;
resumption.step.active=0;
resumption.step.number = 0;
saveas_step.number = 1;
proj_name = 'SingleFracture2.gid';

%% problem information



prob_info.proj_name    = proj_name;
prob_info.proj_path = strcat(proj_path,'\',proj_name);
prob_info = readDatFile(prob_info);

%% write boudary conditions

clear boundary;
boundary.Ux.count = 1;
boundary.Ux.text{1} = '8 , 0.0, 0.0, 0.0';

boundary.Uy.count = 2;
boundary.Uy.text{1} = '7 , 0.0, 0.0, 0.0';
boundary.Uy.text{2} = '8 , 0.0, 0.0, 0.0';

boundary.Sn.count = 1;
boundary.Sn.text{1} = strcat('4,3,',num2str(Sn_load),',0,0');

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
Parameter.specpara.stepnum     = saveas_step.number; % savee all the results in this step
Parameter.specpara.boundaryforce = 0;

% Load parameters

Parameter.load.resumption.active = resumption.active; 
Parameter.load.maxratio          = 1.0;
Parameter.load.volumeforce.active= 0;
Parameter.load.volumeforce.gx    = 0;
Parameter.load.volumeforce.gy    = -0.0981;
Parameter.load.resumption.stepnum = resumption.step.number;    % 91 Step num for resumption
Parameter.load.resumption.stepactive = resumption.step.active;       % 99

% Calculation parameters
Parameter.calpara.loadincrement         = 1000;
Parameter.calpara.itermax               = 1000;
Parameter.calpara.tolerance.criteria    = 1e-6;
Parameter.calpara.tolerance.convergence = 1e-6;
Parameter.calpara.tolerance.displacement= 1e-4;

Parameter.calpara.time.start            = 0;
Parameter.calpara.time.end              = 1;
Parameter.calpara.time.increment        = 1e-1;


%% write material parameters


% Parameter file
Material.proj_path = prob_info.proj_path;
Material.proj_name = prob_info.proj_name;

Material.total_number = 2;

% LinearElastic

Material.type{1}.name = 'LinearElastic';
Material.type{1}.nature = 30000;
Material.type{1}.mecha.modelnum = 31100;
Material.type{1}.mecha.numPara  = 2;
Material.type{1}.mecha.Para(1)  = 1e10;  % Young's modulus
Material.type{1}.mecha.Para(2)  = 0.25;  % Poisson
Material.type{1}.hydro.modelnum = 0;
Material.type{1}.hydro.numPara  = 0;
Material.type{1}.couplingPar(1) = 0;    % selfweight
Material.type{1}.couplingPar(2) = 0;  % biot's coefficient
Material.type{1}.couplingPar(3) = 0;    % thermoespansion


% Fracture1

% Fracture1

Material.type{2}.name = 'Fracture1';
Material.type{2}.nature = 20000;
Material.type{2}.mecha.modelnum  = 21510;
Material.type{2}.mecha.numPara   = 12;
Material.type{2}.mecha.Para(1)   = Kt;      % Kt 
Material.type{2}.mecha.Para(2)   = Kn;      % Kn
Material.type{2}.mecha.Para(3)   = e;       % e
Material.type{2}.mecha.Para(4)   = sigma_R;    % sigma_R
Material.type{2}.mecha.Para(5)   = Cohesion;    % C
Material.type{2}.mecha.Para(6)   = phi;       % phi
Material.type{2}.mecha.Para(7)   = hr;     % hr
Material.type{2}.mecha.Para(8)   = beta;      % beta
Material.type{2}.mecha.Para(9)   = beta_angle;        % beta'
Material.type{2}.mecha.Para(10)  = k0t;      % k0t
Material.type{2}.mecha.Para(11)  = k0n;      % knt
Material.type{2}.mecha.Para(12)  = 1;        % plasticity

Material.type{2}.hydro.modelnum = 0;
Material.type{2}.hydro.numPara  = 0;
Material.type{2}.couplingPar(1) = 0;    % selfweight
Material.type{2}.couplingPar(2) = 0;  % biot's coefficient
Material.type{2}.couplingPar(3) = 0;    % thermoespansion

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
