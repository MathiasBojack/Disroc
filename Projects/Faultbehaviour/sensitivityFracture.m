%%

% sR          = 1.2e6;
% c           = 5.8e6;
% phi         = 31;
% hr          = 0.33;
% beta        = 0.1;
% beta_angle  = 1;

beta_list                 = kron(ones(1,4),[1,0.1,1.8])';
beta_angle_list           = ones(1,12)';
plasticity_indicator_list = [zeros(1,3), ones(1,3), zeros(1,3), ones(1,3)]';
sn_indicator_list         = [103*ones(1,6), 104*ones(1,6)]';

for i = 1:length(beta_list)
    
sR           = 1.2e6;
c            = 5.8e6;
phi          = 31;
hr           = 0.33;
beta         = beta_list(i);
beta_angle   = beta_angle_list(i);
sn_indicator = sn_indicator_list(i);

%% write calculation parameters
fid = fopen('prj_name.txt','r');
l1  = fscanf(fid,'%s');
fclose(fid);
% Parameter file
Parameter.project_path = strcat(project_path,'\',l1);
Parameter.project_name = l1;
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
 Parameter.specpara.stepnum     = 2; % savee all the results in this step
 Parameter.specpara.boundaryforce = 0;

% Load parameters

Parameter.load.resumption.active = 1;
Parameter.load.maxratio          = 1.0;
Parameter.load.volumeforce.active= 0;
Parameter.load.volumeforce.gx    = 0;
Parameter.load.volumeforce.gy    = -0.0981;
Parameter.load.resumption.stepnum = sn_indicator;    % 91 Step num for resumption
Parameter.load.resumption.stepactive = 1;       % 99

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
Material.project_path = strcat(project_path,'\',l1);
Material.project_name = l1;


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
Material.type{2}.mecha.Para(4)   = sR;    % sigma_R 1.2e6
Material.type{2}.mecha.Para(5)   = c;    % C   5.8e6
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
Material.type{2}.couplingPar(2) = 0;  % biot's coefficient
Material.type{2}.couplingPar(3) = 0;    % thermoespansion


%% run Disroc

runDisroc(Parameter,Material,Disroc_path)

%% plot for single Fracture

TimeIncrement = 1;
TimeEndRatio = 10;
plotJointElemNo = 5;
%
fnameJointMecha = strcat(path,'\','201.jointMecha.dat');
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
% calculate the damage criterion

crit_Parameter.beta       = beta;
crit_Parameter.beta_angle = beta_angle;
crit_Parameter.hr         = hr;
crit_Parameter.Damage     = Damage;
crit_Parameter.c          = c;
crit_Parameter.sR         = sR;
crit_Parameter.phi        = phi;
crit_Parameter.tau        = Tau;
crit_Parameter.sn         = Sn;
Fun_crit                  = damageCriterion(crit_Parameter);
% Time evolution plot
isPlot = NoElem==plotJointElemNo;
f1 = figure(1); 
clf;hold on;
plot(Time(isPlot), Tau(isPlot)/1e6,'r')
% plot asymptotic line
% plot([Time(1),Time(end)], Tau(end)*[1,1]/1e6)
% plot(Time(isPlot), Sn(isPlot)/1e6,'b')
% % plot(Time(isPlot), Tau(isPlot)./Sn(isPlot),'-k')
% plot(Time(isPlot), Utp(isPlot)/Utp(end))
plot(Time(isPlot), Damage(isPlot))
% plot(Time(isPlot), Fun_crit(isPlot))

ax = gca;
leg_txt{1} = '$Tau$';
leg_txt{2} = '$Damage$';
xlabel('Load','interpreter','latex')
ylabel('Shear stress','interpreter','latex')
legend(leg_txt,'interpreter','latex')

figname = strcat('F:\13.Code\Git\Disroc\Projects\Faultbehaviour\SingleFracutre',num2str(i),'.png') ;
saveas(f1,figname)

end