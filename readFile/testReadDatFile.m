fid = fopen('prj_name.txt','r');
proj_name  = fscanf(fid,'%s');
fclose(fid);

prob_info.proj_name    = proj_name;
prob_info.proj_path = strcat(proj_path,'\',proj_name);

prob_info = readDatFile(prob_info);

