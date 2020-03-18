function prob_info = readDatFile(prob_info)


proj_name       = erase(prob_info.proj_name,'.gid');
project_path    = prob_info.project_path;

fname = strcat(project_path,'\',proj_name,'.gid','\',proj_name,'.dat');

fid = fopen(fname,'r+');

l1              = textscan(fid,'%s',1);
l2_1              = textscan(fid,'%d',5, 'Delimiter',{','});
l2_2              = textscan(fid,'%s',5, 'Delimiter',{','});
numNode         = l2_1{1}(1);
numElem         = l2_1{1}(3);
numJointElem    = l2_1{1}(4); 
numMat          = l2_1{1}(5); 
l3   = textscan(fid,'%s',1);  % coordinates
l4   = textscan(fid, '%d %f %f', numNode,'Delimiter',{','});

prob_info.coordinate = zeros(numNode,2);
prob_info.coordinate(:,1) = l4{2}(:);
prob_info.coordinate(:,2) = l4{3}(:);

l5   = textscan(fid,'%s',1);  % connectivity
l6   = textscan(fid, '%d %d %d %d %d %d ', numElem,'Delimiter',{','});
prob_info.connectivity.matrix = zeros(numElem, 3);
prob_info.connectivity.matrix(:,1) = l6{3}(:);
prob_info.connectivity.matrix(:,2) = l6{4}(:);
prob_info.connectivity.matrix(:,3) = l6{5}(:);

if numJointElem >0
    l7 = textscan(fid, '%d %d %d %d %d %d %d', numJointElem,'Delimiter',{','});
    prob_info.connectivity.joint = zeros(numJointElem, 4);
    prob_info.connectivity.joint(:,1) = l7{3}(:);
    prob_info.connectivity.joint(:,2) = l7{4}(:);
    prob_info.connectivity.joint(:,3) = l7{5}(:);
    prob_info.connectivity.joint(:,3) = l7{6}(:);
end
l8 = textscan(fid,'%s',1);  % materials
l9 = textscan(fid,'%d %d %s',numMat,'Delimiter',',');
end