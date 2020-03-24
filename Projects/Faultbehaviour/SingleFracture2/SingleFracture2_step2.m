SingleFracture2_step1

Imposed_disp = 1e-1;

%% Resumption
resumption.active      = 1;
resumption.step.active = 1;
resumption.step.number = 1;
saveas_step.number     = 2;
%% write boudary conditions

clear boundary;
boundary.Ux.count = 4;
joint_nodes_I = prob_info.connectivity.joint(:,1);
joint_nodes_J = prob_info.connectivity.joint(:,2);
joint_lower   = union(joint_nodes_I,joint_nodes_J);
joint_nodes_K = prob_info.connectivity.joint(:,3);
joint_nodes_L = prob_info.connectivity.joint(:,4);
joint_upper   = union(joint_nodes_K,joint_nodes_L);

boundary.Ux.text{1} = ['3 ,' num2str(Imposed_disp) ', 0.0, 0.0'];
boundary.Ux.text{2} = ['4 ,' num2str(Imposed_disp) ', 0.0, 0.0'];
boundary.Ux.text{3} = '7 , 0.0, 0.0, 0.0';
boundary.Ux.text{4} = '8 , 0.0, 0.0, 0.0';

% boundary.Ux.count = 2;
% boundary.Ux.text{1} = '7 , 0.0, 0.0, 0.0';
% boundary.Ux.text{2} = '8 , 0.0, 0.0, 0.0';
% boundary.Tau.count = 1;
% boundary.Tau.text{1} = '3,4, 23.3e6,0,0';


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

TimeIncrement = 1;
TimeEndRatio = 1;
plotJointElemNo = 5;
%
fnameJointMecha = strcat(Material.proj_path,'\','201.jointMecha.dat');
fidJointMecha   = fopen(fnameJointMecha,'r');
cellJoinMecha   = textscan(fidJointMecha, '%f %d %f %f %f %f %f %f %f','Headerlines',1);
TimeEnd = floor(length(cellJoinMecha{1})/TimeEndRatio);
Time    = cellJoinMecha{1}(1:TimeIncrement:TimeEnd)*Imposed_disp;
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
plot(Ut(isPlot), Tau(isPlot)/1e6,'r')
% plot asymptotic line
% plot([Time(1),Time(end)], Tau(end)*[1,1]/1e6)
plot(Time(isPlot), Sn(isPlot)/1e6,'b')
% % plot(Time(isPlot), Tau(isPlot)./Sn(isPlot),'-k')
% plot(Time(isPlot), Utp(isPlot)/Utp(end))
% plot(Time(isPlot), Damage(isPlot))
% plot(Time(isPlot), Fun_crit(isPlot))



xlabel('Shear displacement [/m]','interpreter','latex')
ylabel('Stress [/MPa]','interpreter','latex')
title('Fault with Cohesive damage criteria','interpreter','latex')
legtex{1} = '$\tau$';
legtex{2} = '$\sigma_n$';
legend(legtex,'interpreter','latex','Location','best');
grid on;
box on;
% saveas(f1,'.\Projects\Faultbehaviour\SingleFracture2\CohesiveFracture.pdf')

%% plot joint results along the fault


faultplot.X = prob_info.coordinate(joint_nodes_I,1);
numStepCalculated = length(cellJoinMecha{1})/prob_info.numJointElem;
faultplot.stepInitial = 1;
faultplot.indexInitial = ((faultplot.stepInitial-1)*prob_info.numJointElem + (1:prob_info.numJointElem))';
faultplot.UtInitial   = cellJoinMecha{3}(faultplot.indexInitial);
faultplot.TauInitial   = cellJoinMecha{5}(faultplot.indexInitial);
faultplot.SnInitial   = cellJoinMecha{6}(faultplot.indexInitial);

figure(2);clf;hold on;
ylabel('$Ut$ [/m]','interpreter','latex')

figure(3);clf;hold on;
ylabel('$\tau$ [/MPa]','interpreter','latex')

figure(4);clf;hold on;
ylabel('$\sigma_n$ [/MPa]','interpreter','latex')

step_list = [1 10 20 34];
for i = 1:length(step_list)
    faultplot.step = step_list(i);
    faultplot.index = ((faultplot.step-1)*prob_info.numJointElem + (1:prob_info.numJointElem))';
    faultplot.Ut   = cellJoinMecha{3}(faultplot.index);
    faultplot.Un   = cellJoinMecha{4}(faultplot.index);
    faultplot.Tau  = cellJoinMecha{5}(faultplot.index);
    faultplot.Sn   = cellJoinMecha{6}(faultplot.index);
    figure(2);plot(faultplot.Ut - faultplot.UtInitial,faultplot.X)
    figure(3);plot((faultplot.Tau -faultplot.TauInitial)/1e6 ,faultplot.X);
    figure(4);plot((faultplot.Sn-faultplot.SnInitial)/1e6 ,faultplot.X);
end