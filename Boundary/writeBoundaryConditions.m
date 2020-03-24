function writeBoundaryConditions( prob_info,boundary )
%This function writes the boundary condition for the Disroc dat file
%   Detailed explanation goes here

fname = strcat(prob_info.proj_path,'\',erase(prob_info.proj_name,'.gid'),'.dat');
fid = fopen(fname,'r+');

fname_temp = strcat(prob_info.proj_path,'\',erase(prob_info.proj_name,'.gid'),'_temp.dat');
fid_temp = fopen(fname_temp,'w+');


I = prob_info.numNode + prob_info.numElem + prob_info.numJointElem + prob_info.numMat + 6;
for i = 1: I
    tline = fgets(fid);
    fprintf(fid_temp,'%s',tline);
end

%% Ux


if isfield(boundary,'Ux')==0 || (boundary.Ux.count ==0)
    tline = '1, 0  -----------Ux';
    fprintf(fid_temp,'%s',tline);fprintf(fid_temp,'\n');
elseif boundary.Ux.count~=0
    tline = ['1, ', num2str(boundary.Ux.count), ' -----------Ux'];
    fprintf(fid_temp,'%s',tline);fprintf(fid_temp,'\n');
    for k = 1:boundary.Ux.count
        tline = boundary.Ux.text{k};
        fprintf(fid_temp,'%s',tline);fprintf(fid_temp,'\n');
    end
end

%% Uy
if isfield(boundary,'Uy')==0 || (boundary.Uy.count ==0)
    tline = '2, 0  -----------Uy';
    fprintf(fid_temp,'%s',tline);fprintf(fid_temp,'\n');
elseif boundary.Uy.count~=0
    tline = ['2, ', num2str(boundary.Uy.count), ' -----------Uy'];
    fprintf(fid_temp,'%s',tline);fprintf(fid_temp,'\n');
    for k = 1:boundary.Uy.count
        tline = boundary.Uy.text{k};
        fprintf(fid_temp,'%s',tline);fprintf(fid_temp,'\n');
    end
end
%% P


if isfield(boundary,'P')==0 || (boundary.P.count ==0)
    tline = '3, 0  -----------P';
    fprintf(fid_temp,'%s',tline);fprintf(fid_temp,'\n');
elseif boundary.P.count~=0
    tline = ['3, ', num2str(boundary.P.count), ' -----------P'];
    fprintf(fid_temp,'%s',tline);fprintf(fid_temp,'\n');
    for k = 1:boundary.P.count
        tline = boundary.P.text{k};
        fprintf(fid_temp,'%s',tline);fprintf(fid_temp,'\n');
    end
end

%% Sn


if isfield(boundary,'Sn')==0 || (boundary.Sn.count ==0)
    tline = '4, 0  -----------Sn';
    fprintf(fid_temp,'%s',tline);fprintf(fid_temp,'\n');
elseif boundary.Sn.count~=0
    tline = ['4, ', num2str(boundary.Sn.count), ' -----------Sn'];
    fprintf(fid_temp,'%s',tline);fprintf(fid_temp,'\n');
    for k = 1:boundary.Sn.count
        tline = boundary.Sn.text{k};
        fprintf(fid_temp,'%s',tline);fprintf(fid_temp,'\n');
    end
end

%% Tau
if isfield(boundary,'Tau')==0 || (boundary.Tau.count ==0)
    tline = '5, 0  -----------Tau';
    fprintf(fid_temp,'%s',tline);fprintf(fid_temp,'\n');
elseif boundary.Tau.count~=0
    tline = ['5, ', num2str(boundary.Tau.count), ' -----------Tau'];
    fprintf(fid_temp,'%s',tline);fprintf(fid_temp,'\n');
    for k = 1:boundary.Tau.count
        tline = boundary.Tau.text{k};
        fprintf(fid_temp,'%s',tline);fprintf(fid_temp,'\n');
    end
end

%% Vn
tline = '6, 0  -----------Vn';
fprintf(fid_temp,'%s',tline);fprintf(fid_temp,'\n');
%% Fx
tline = '7, 0  -----------Fx';
fprintf(fid_temp,'%s',tline);fprintf(fid_temp,'\n');

%% Fy
tline = '8, 0  -----------Fy';
fprintf(fid_temp,'%s',tline);fprintf(fid_temp,'\n');

%% Qn
tline = '9, 0  -----------Qn';
fprintf(fid_temp,'%s',tline);fprintf(fid_temp,'\n');

%% Theta
tline = '10, 0  -----------Theta';
fprintf(fid_temp,'%s',tline);fprintf(fid_temp,'\n');
%% Bending Moment
tline = '11, 0  -----------Bending Moment';
fprintf(fid_temp,'%s',tline);fprintf(fid_temp,'\n');
%% Temperature
tline = '12, 0  -----------Temperature';
fprintf(fid_temp,'%s',tline);fprintf(fid_temp,'\n');
%% HeatFlux
tline = '13, 0  -----------HeatFlux';
fprintf(fid_temp,'%s',tline);fprintf(fid_temp,'\n');
%% Concentration
tline = '14, 0  -----------Concentration';
fprintf(fid_temp,'%s',tline);fprintf(fid_temp,'\n');
%% ChemicalFlux
tline = '15, 0  -----------ChemicalFlux';
fprintf(fid_temp,'%s',tline);fprintf(fid_temp,'\n');

fprintf(fid_temp,'%s','End Data');

fclose(fid);
fclose(fid_temp);
fclose('all');
movefile(fname_temp,fname)
end

