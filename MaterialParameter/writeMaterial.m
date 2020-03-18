function writeMaterial(Material)

%%
% Material.
%          project_path
%          project_name
%          total_number
%          type{i}.
%                  name
%                  nature
%                  mecha.
%                        modelnum
%                        numPara
%                        Para(1:numPara)
%                  hydro.
%                        modelnum
%                        numPara
%                        Para(1:numPara)
% %                  couplingPar(1:3)
fname = strcat(Material.project_path,'\', 'Mats.dat');
fid = fopen( fname,'w+');

fprintf(fid, '%s\n','Disroc5');
fprintf(fid, '%d\n',Material.total_number);

for i = 1: Material.total_number
    clear l;
    l = cell(Material.type{i}.mecha.numPara + Material.type{i}.hydro.numPara +10,1);
    l{1} = num2str(i);
    fprintf(fid, '%s\n',l{1});
    temp = strcat(Material.type{i}.name,{' '}, num2str(Material.type{i}.nature));
    l{2} = temp{1};
    fprintf(fid, '%s\n',l{2});
    l{3} = 'Mechanics:';
    fprintf(fid, '%s\n',l{3});
    temp = strcat(num2str(Material.type{i}.mecha.modelnum), {' '},num2str(Material.type{i}.mecha.numPara));
    l{4} = temp{1};
    fprintf(fid, '%s\n',l{4});
    for j=1:Material.type{i}.mecha.numPara
        l{j+4} = num2str(Material.type{i}.mecha.Para(j));
        fprintf(fid, '%s\n',l{j+4});
    end
    l{Material.type{i}.mecha.numPara+5} = 'Hydraulic:' ;
    fprintf(fid, '%s\n',l{Material.type{i}.mecha.numPara+5});
    temp = strcat(num2str(Material.type{i}.hydro.modelnum), {' '},num2str(Material.type{i}.hydro.numPara));
    l{Material.type{i}.mecha.numPara+6} = temp{1};
    fprintf(fid, '%s\n',l{Material.type{i}.mecha.numPara+6});
    for j=1:Material.type{i}.hydro.numPara
        l{Material.type{i}.mecha.numPara+6+j} = num2str(Material.type{i}.hydro.Para(j));
        fprintf(fid, '%s\n',l{Material.type{i}.mecha.numPara+6+j});
    end
    l{Material.type{i}.mecha.numPara + Material.type{i}.hydro.numPara + 7}...
        = 'Thermal:';
    fprintf(fid, '%s\n',...
        l{Material.type{i}.mecha.numPara + Material.type{i}.hydro.numPara + 7});
    l{Material.type{i}.mecha.numPara + Material.type{i}.hydro.numPara + 8}...
        = '0 0 ! Model, Nb of parameters';
    fprintf(fid, '%s\n',...
        l{Material.type{i}.mecha.numPara + Material.type{i}.hydro.numPara + 8});
    l{Material.type{i}.mecha.numPara + Material.type{i}.hydro.numPara + 9}...
        = 'Couplings-rho-biot-alpha:';    
    fprintf(fid, '%s\n',...
        l{Material.type{i}.mecha.numPara + Material.type{i}.hydro.numPara + 9});
    temp = strcat(num2str(Material.type{i}.nature),{' '},'3');
    l{Material.type{i}.mecha.numPara + Material.type{i}.hydro.numPara + 10}...
        = temp{1};
    fprintf(fid, '%s\n',...
        l{Material.type{i}.mecha.numPara + Material.type{i}.hydro.numPara + 10});
    
    for j =1 :3
        fprintf(fid, '%s\n', num2str(Material.type{i}.couplingPar(j)) ) ;
    end
    
end

fclose(fid);


end