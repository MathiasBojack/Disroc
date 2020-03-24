clear;
Path_Control

% Step 1, applied normal stress
%% 
sigma_R      = 1.2e6;
tau_c        = 5.8e6;
phi          = 31;
Cohesion   = sqrt(2*sigma_R*tau_c*tan(phi/180*pi) -(sigma_R*tan(phi/180*pi))^2);    % C

hr           = 0.33;
beta         = 0.1;
beta_angle   = 1;


Sn_load  = 30e6;   % The applied normal stress is considered positive when the direction is aligned with the outer normal
% Tau_load = -10e6; % The applied shear stress is considered positve when Ut is positive.
resumption.active = 0;
resumption.step.active=0;
resumption.step.number = 0;
saveas_step.number = 1;
proj_name = 'Induced101.gid';
%% problem information

prob_info.proj_name    = proj_name;
prob_info.proj_path = strcat(proj_path,'\',proj_name);
prob_info = readDatFile(prob_info);


%% boundary conditions

clear boudnary
eps = 1e-5;
joint_nodes_I = prob_info.connectivity.joint(:,1);
joint_nodes_J = prob_info.connectivity.joint(:,2);
joint_lower   = union(joint_nodes_I,joint_nodes_J);
joint_nodes_K = prob_info.connectivity.joint(:,3);
joint_nodes_L = prob_info.connectivity.joint(:,4);
joint_upper   = union(joint_nodes_K,joint_nodes_L);
% figure(1);hold on;
% plot(prob_info.coordinate(joint_nodes_I,1),prob_info.coordinate(joint_nodes_I,2),'*');
% plot(prob_info.coordinate(joint_nodes_J,1),prob_info.coordinate(joint_nodes_J,2),'o');
% plot(prob_info.coordinate(joint_nodes_K,1),prob_info.coordinate(joint_nodes_K,2),'s');
% plot(prob_info.coordinate(joint_nodes_L,1),prob_info.coordinate(joint_nodes_L,2),'d');

joint_nodes_I_left = joint_nodes_I(abs(prob_info.coordinate(joint_nodes_I,1))<eps);

boundary.Ux.count = 1;
boundary.Ux.text{1} = [ num2str(joint_nodes_I_left) ', 0.0, 0.0, 0.0'];
boundary.Uy.count = length(joint_lower);
for i=1:length(joint_lower)
    boundary.Uy.text{i} = [ num2str(joint_lower(i)) ', 0.0, 0.0, 0.0'];
end

% boundary.Uy.count = length(joint_lower)*2;
% boundary.Uy.text = cell(length(joint_lower)*2,1);
% for i=1:length(joint_lower)
%     boundary.Uy.text{i} = [ num2str(joint_lower(i)) ', 0.0, 0.0, 0.0'];
%     boundary.Uy.text{length(joint_lower) + i} = [ num2str(joint_upper(i)) ', ' num2str(5e-1) ', 0.0, 0.0'];
% end

boundary.Sn.count = prob_info.numJointElem;
for i = 1 : prob_info.numJointElem
    boundary.Sn.text{i} = [num2str(joint_nodes_L(i)) ', ' num2str(joint_nodes_K(i)) ...
        ', ' num2str(Sn_load),',0,0'];
end

% nodes_list = 1:prob_info.numNode;
% nodes_top = nodes_list(abs(prob_info.coordinate(:,2)- 5000)<eps);
% [~,sortedIndex] = sort(prob_info.coordinate(nodes_top,1));
% nodes_top_sorted = nodes_top(sortedIndex);
% boundary.Tau.count = length(nodes_top_sorted)-1;
% for i = 1 : (length(nodes_top_sorted)-1)
%     boundary.Tau.text{i} = [num2str(nodes_top_sorted(i)) ', ' num2str(nodes_top_sorted(i+1)) ...
%         ', ' num2str(Tau_load),',0,0'];
% end




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
Material.type{1}.mecha.Para(1)  = 1e10;  % Young's modulus
Material.type{1}.mecha.Para(2)  = 0.25;  % Poisson
Material.type{1}.hydro.modelnum = 0;
Material.type{1}.hydro.numPara  = 0;
Material.type{1}.couplingPar(1) = 2;    % selfweight
Material.type{1}.couplingPar(2) = 0.5;  % biot's coefficient
Material.type{1}.couplingPar(3) = 0;    % thermoespansion


% Fracture1

Material.type{2}.name = 'Fracture1';
Material.type{2}.nature = 20000;
Material.type{2}.mecha.modelnum  = 21510;
Material.type{2}.mecha.numPara   = 12;
Material.type{2}.mecha.Para(1)   = 2e8;      % Kt 
Material.type{2}.mecha.Para(2)   = 5e8;      % Kn
Material.type{2}.mecha.Para(3)   = 10;       % e
Material.type{2}.mecha.Para(4)   = sigma_R;    % sigma_R 1.2e6
Material.type{2}.mecha.Para(5)   = Cohesion;    % C   2.8e6 => tau_c = 5.8
Material.type{2}.mecha.Para(6)   = phi;       % phi
Material.type{2}.mecha.Para(7)   = hr;     % h  0.33
Material.type{2}.mecha.Para(8)   = beta;      % beta  0.1
Material.type{2}.mecha.Para(9)   = beta_angle;        % beta' 1
Material.type{2}.mecha.Para(10)  = 1e8;      % k0t
Material.type{2}.mecha.Para(11)  = 1e10;      % knt
Material.type{2}.mecha.Para(12)  = 1;        % 

Material.type{2}.hydro.modelnum = 0;
Material.type{2}.hydro.numPara  = 0;
Material.type{2}.couplingPar(1) = 0;    % selfweight
Material.type{2}.couplingPar(2) = 1;  % biot's coefficient
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
%% Plot

TimeIncrement = 1;
TimeEndRatio = 1;
plotJointElemNo = 1282;
%
fnameJointMecha = strcat(Material.proj_path,'\','201.jointMecha.dat');
fidJointMecha   = fopen(fnameJointMecha,'r');
cellJoinMecha   = textscan(fidJointMecha, '%f %d %f %f %f %f %f %f %f','Headerlines',1);
TimeEnd = floor(length(cellJoinMecha{1})/TimeEndRatio);
Time    = cellJoinMecha{1}(1:TimeIncrement:TimeEnd);
NoElem  = cellJoinMecha{2}(1:TimeIncrement:TimeEnd);
Ut      = cellJoinMecha{3}(1:TimeIncrement:TimeEnd);
Un      = cellJoinMecha{4}(1:TimeIncrement:TimeEnd);
Tau     = cellJoinMecha{5}(1:TimeIncrement:TimeEnd);
Sn      = cellJoinMecha{6}(1:TimeIncrement:TimeEnd);
Utp     = cellJoinMecha{7}(1:TimeIncrement:TimeEnd);
Unp     = cellJoinMecha{8}(1:TimeIncrement:TimeEnd);
Damage  = cellJoinMecha{9}(1:TimeIncrement:TimeEnd);
fclose(fidJointMecha);

crit_Parameter.beta       = beta;
crit_Parameter.beta_angle = beta_angle;
crit_Parameter.hr         = hr;
crit_Parameter.Damage     = Damage;
crit_Parameter.c          = Cohesion;
crit_Parameter.sR         = sigma_R;
crit_Parameter.phi        = phi;

crit_Parameter.tau        = Tau;
crit_Parameter.sn         = Sn;
Fun_crit                  = damageCriterion(crit_Parameter);
% Time evolution plot
isPlot = NoElem==plotJointElemNo;
f1 = figure(1); 
clf;
hold on;
plot(Time(isPlot), Tau(isPlot)/1e6,'r')
% plot asymptotic line
% plot([Time(1),Time(end)], Tau(end)*[1,1]/1e6)
plot(Time(isPlot), Sn(isPlot)/1e6,'b')
% % plot(Time(isPlot), Tau(isPlot)./Sn(isPlot),'-k')
% plot(Time(isPlot), Utp(isPlot)/Utp(end))
% plot(Time(isPlot), Damage(isPlot))
plot(Time(isPlot), Fun_crit(isPlot))



xlabel('Shear displacement [/m]','interpreter','latex')
ylabel('Stress [/MPa]','interpreter','latex')
% title('Fault with Mohr-Coulomb Failure criteria','interpreter','latex')
legtex{1} = '$\tau$';
legtex{2} = '$\sigma_n$';
legend(legtex,'interpreter','latex','Location','best');
grid on;
box on;
% saveas(f1,'.\Projects\Faultbehaviour\SingleFracture2\CohesiveFracture.pdf')
