% File operation
fid = fopen('prj_name.txt','r');
l1  = fscanf(fid,'%s');
fclose(fid);
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
Material.type{1}.mecha.Para(2)  = 0.25;
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
Material.type{2}.mecha.Para(4)   = 1.2e6;    % sigma_R
Material.type{2}.mecha.Para(5)   = 5.8e6;    % C
Material.type{2}.mecha.Para(6)   = 31;       % phi
Material.type{2}.mecha.Para(7)   = 0.33;     % h
Material.type{2}.mecha.Para(8)   = 0.1;      % beta
Material.type{2}.mecha.Para(9)   = 1;        % beta'
Material.type{2}.mecha.Para(10)  = 1e7;      % k0t
Material.type{2}.mecha.Para(11)  = 1e9;      % knt
Material.type{2}.mecha.Para(12)  = 1;        % 

Material.type{2}.hydro.modelnum = 0;
Material.type{2}.hydro.numPara  = 0;
Material.type{2}.couplingPar(1) = 0;    % selfweight
Material.type{2}.couplingPar(2) = 0;  % biot's coefficient
Material.type{2}.couplingPar(3) = 0;    % thermoespansion

writeMaterial(Material)

