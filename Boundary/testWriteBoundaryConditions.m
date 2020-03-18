% fid = fopen('prj_name.txt','r');
% proj_name  = fscanf(fid,'%s');
% fclose(fid);
proj_name = 'SingleFracture0.gid';

prob_info.proj_name    = proj_name;
prob_info.proj_path = strcat(proj_path,'\',proj_name);

prob_info = readDatFile(prob_info);

%%

boundary.Ux.count = 4;
boundary.Ux.text{1} = '3 , 0.5e1, 0.0, 0.0';
boundary.Ux.text{2} = '4 , 0.5e1, 0.0, 0.0';
boundary.Ux.text{3} = '7 , 0.0, 0.0, 0.0';
boundary.Ux.text{4} = '8 , 0.0, 0.0, 0.0';

boundary.Uy.count = 2;
boundary.Uy.text{1} = '7 , 0.0, 0.0, 0.0';
boundary.Uy.text{2} = '8 , 0.0, 0.0, 0.0';

writeBoundaryConditions( prob_info,boundary )
