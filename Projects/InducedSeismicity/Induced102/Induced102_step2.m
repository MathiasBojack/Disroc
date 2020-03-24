clear;clc;
Path_Control
% Induced101_step1
% Step 1, applied normal stress
%% Elastic fracture for stress initiation 



P_load         = 19e6;  

% Matrix material parameter
E              = 1e10;  % young's modulus
nu             = 0.25;
k_permeability = 100 * 1e-15; % mDarcy
fluids         = 'water';
phi_porosity   = 0.1;
Ks             = 4e10; 

switch fluids
    case 'water'
        mu_viscosity   = 5e-4;  % Pa.s 
        Kf             = 2e9;   % Pa
    case 'CO2'
        mu_viscosity   = 5e-5;  % Pa.s 
        Kf             = 2e8;   % Pa
end
K_bulk         = E/3/(1-2*nu);
biot           = 1 - K_bulk/Ks;
k_intrinsic    = k_permeability/mu_viscosity;
CM_bulk        = phi_porosity/Kf +(biot - phi_porosity)/Ks;

% Fracture material parameter
Kt             = 2e8;      % Kt 
Kn             = 5e8;      % Kn
e              = 10;       % e
sigma_R        = 1.2e6;    % sigma_R
tau_c          = 5.8e6;

phi            = 31;       % phi
Cohesion       = sqrt(2*sigma_R*tau_c*tan(phi/180*pi) -(sigma_R*tan(phi/180*pi))^2);    % C
hr             = 0.33;     % h
beta           = 0.1;      % beta
beta_angle     = 1;        % beta'
k0t            = 1e8;      % k0t
k0n            = 1e10;      % k0n

Ct_fracture    = k_intrinsic*e;
Cn_fracture    = k_intrinsic/e;
CM_fracture    = CM_bulk*e;

%% Resumption
resumption.active      = 1;
resumption.step.active = 1;
resumption.step.number = 1;
saveas_step.number     = 2;
proj_name = 'Induced102.gid';
%% problem information

prob_info.proj_name    = proj_name;
prob_info.proj_path    = strcat(proj_path,'\',proj_name);
prob_info              = readDatFile(prob_info);



%% boundary conditions

clear boudnary
eps = 1e-1;
joint_nodes_I = prob_info.connectivity.joint(:,1);
joint_nodes_J = prob_info.connectivity.joint(:,2);
joint_lower   = union(joint_nodes_I,joint_nodes_J);
joint_nodes_K = prob_info.connectivity.joint(:,3);
joint_nodes_L = prob_info.connectivity.joint(:,4);
joint_upper   = union(joint_nodes_K,joint_nodes_L);

nodeList          = 1:prob_info.numNode;

% bottom
nodes_bottomBC      = nodeList(abs(prob_info.coordinate(:,2))<eps);
boundary.Uy.count = length(nodes_bottomBC);

for i = 1:length(nodes_bottomBC)
     boundary.Uy.text{i} = [ num2str(nodes_bottomBC(i)) ', 0.0, 0.0, 0.0'];
end

% Left
nodes_leftBC        = nodeList(abs(prob_info.coordinate(:,1))<eps);
boundary.Ux.count = length(nodes_leftBC);
for i=1:length(nodes_leftBC)
     boundary.Ux.text{i} = [ num2str(nodes_leftBC(i)) ', 0.0, 0.0, 0.0'];
end

% Pressure prescribed on all four sides

nodes_topBC        = nodeList(abs(prob_info.coordinate(:,2)-5000)<eps);
nodes_rightBC      = nodeList(abs(prob_info.coordinate(:,1)-5000)<eps);

nodes_pressure     = union_several(nodes_bottomBC, nodes_leftBC,...
                                    nodes_topBC, nodes_rightBC);

boundary.P.count = length( nodes_pressure);

for i=1:length(nodes_pressure)
     boundary.P.text{i} = [ num2str(nodes_pressure(i)) ',' ...
          ',' num2str(P_load) ', 0.0, 0.0'];
end



writeBoundaryConditions( prob_info,boundary )

%% write calculation parameters

% Parameter file
Parameter.proj_path = prob_info.proj_path;
Parameter.proj_name = prob_info.proj_name;
% Problem type information
Parameter.problem.physics      = 2; %1:Mechanics, 2:Hydraulic, 3:Thermal, 4:HM, 5:TM, 6:THM
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
Parameter.calpara.loadincrement         = 10;
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
Material.type{1}.mecha.Para(1)  = E;  % Young's modulus
Material.type{1}.mecha.Para(2)  = nu;  % Poisson
% 
Material.type{1}.hydro.modelnum = 32110;
Material.type{1}.hydro.numPara  = 2;
Material.type{1}.hydro.Para(1)  = k_intrinsic;
Material.type{1}.hydro.Para(2)  = CM_bulk; 

Material.type{1}.couplingPar(1) = 0;    % selfweight
Material.type{1}.couplingPar(2) = biot;  % biot's coefficient
Material.type{1}.couplingPar(3) = 0;    % thermoexpansion

% % Fracture 21510
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

% Material.type{2}.hydro.modelnum = 22210;
% Material.type{2}.hydro.numPara  = 3;
% Material.type{2}.hydro.Para(1)  = Ct_fracture;
% Material.type{2}.hydro.Para(2)  = Cn_fracture;
% Material.type{2}.hydro.Para(3)  = CM_fracture;

Material.type{2}.hydro.modelnum = 22110;
Material.type{2}.hydro.numPara  = 2;
Material.type{2}.hydro.Para(1)  = Ct_fracture;
Material.type{2}.hydro.Para(2)  = CM_fracture;


Material.type{2}.couplingPar(1) = 0;    % selfweight
Material.type{2}.couplingPar(2) = biot;  % biot's coefficient
Material.type{2}.couplingPar(3) = 0;    % thermoespansion




%% run Disroc

runDisroc(Parameter,Material,Disroc_path)

%% save to step 1

foldername = [Parameter.proj_path,'\STEP-',num2str(saveas_step.number)];
if not(exist(foldername,'dir'))
        mkdir(foldername)
end


cmd_txt = ['copy', ' ', Parameter.proj_path,'\',...
     erase(Parameter.proj_name,'gid'), '*', ' ',...
     foldername,'\',erase(Parameter.proj_name,'gid'), '*'];
system(cmd_txt)


cmd_txt = ['copy', ' ', Parameter.proj_path,'\',...
     '*.dat', ' ', foldername,'\', '*.dat'];
system(cmd_txt)

%% plot joint results along the fault

% fnameJointHydro = strcat(proj_path,'\',proj_name,...
%     '\STEP-',num2str(saveas_step.number),'\','203.jointHydro.dat');
% fidJointHydro   = fopen(fnameJointHydro,'r');
% cellJoinHydro   = textscan(fidJointHydro, '%f %d %f %f','Headerlines',1);
% 
% faultplot.X = prob_info.coordinate(joint_nodes_I,1);
% faultplot.Y = prob_info.coordinate(joint_nodes_I,2);
% 
% 
% figure(1);clf;hold on;
% xlabel('$\sigma _n$','interpreter','latex')
% ylabel('$Y$','interpreter','latex')
% 
% 
%     faultplot.step = step_list(i);
%     faultplot.index = ((faultplot.step-1)*prob_info.numJointElem + (1:prob_info.numJointElem))';
%     faultplot.Ut   = cellJoinHydro{3}(faultplot.index);
%     faultplot.Un   = cellJoinHydro{4}(faultplot.index);
%     faultplot.Tau  = cellJoinHydro{5}(faultplot.index);
%     faultplot.Sn   = cellJoinHydro{6}(faultplot.index);
%     figure(2);plot(faultplot.Sn/1e6,faultplot.Y)
%     figure(3);plot(faultplot.Tau/1e6 ,faultplot.Y)
%     figure(4);plot(faultplot.Tau/faultplot.Sn ,faultplot.Y)
%     figure(5);plot(faultplot.Ut ,faultplot.Y)
