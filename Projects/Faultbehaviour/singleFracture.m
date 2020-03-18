clc

project_name = 'SingleFracture0';
path = strcat(project_path,'\',project_name,'.gid');

% command = strcat(Disroc_path,'\uDisroc.exe',{' '},project_name);
% command_exec = command{1};
% system(command_exec)
%%
TimeIncrement = 1;
TimeEndRatio = 1;
plotJointElemNo = 5;
%%
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

% Time evolution plot
isPlot = NoElem==plotJointElemNo;
f1 = figure(1); hold on;
plot(Time(isPlot), Tau(isPlot)/1e6,'r')
% plot asymptotic line
% plot([Time(1),Time(end)], Tau(end)*[1,1]/1e6)
plot(Time(isPlot), Sn(isPlot)/1e6,'b')
% % plot(Time(isPlot), Tau(isPlot)./Sn(isPlot),'-k')
plot(Time(isPlot), Utp(isPlot)/Utp(end))
plot(Time(isPlot), Damage(isPlot))


%% Optimized parameters:


% Material.type{2}.mecha.Para(1)   = 2e8;      % Kt 
% Material.type{2}.mecha.Para(2)   = 5e8;      % Kn
% Material.type{2}.mecha.Para(3)   = 10;       % e
% Material.type{2}.mecha.Para(4)   = 1.2e6;    % sigma_R
% Material.type{2}.mecha.Para(5)   = 5.8e6;    % C
% Material.type{2}.mecha.Para(6)   = 31;       % phi
% Material.type{2}.mecha.Para(7)   = 0.33;     % h
% Material.type{2}.mecha.Para(8)   = 0.2;      % beta
% Material.type{2}.mecha.Para(9)   = 1;        % beta'
% Material.type{2}.mecha.Para(10)  = 1e7;      % k0t
% Material.type{2}.mecha.Para(11)  = 1e9;      % knt
% Material.type{2}.mecha.Para(12)  = 1;

