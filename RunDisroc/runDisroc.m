function runDisroc(Parameter,Material,Disroc_path)

writeParam(Parameter)
writeMaterial(Material)
temp = strcat('.\uDisroc.exe', {' '}, erase(Parameter.proj_name,'.gid')); 

path = pwd;
cd(Disroc_path);

 
cmdtext = temp{1};
system(cmdtext)

cd(path)
end