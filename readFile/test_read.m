clc
y = [(1:12)', rand(12,1)]';
fid = fopen('.\readFile\test_read.txt','w+');
fprintf(fid,'%d %4.2f\n', y);
fclose(fid);
%%
fid = fopen('.\readFile\test_read.txt','r+');
% ftell(fid)
% fl  = fgetl(fid)
% ftell(fid)
% fl  = fgetl(fid)
% ftell(fid)

% fl  = fgetl(fid)
% fl  = fgetl(fid)
% fl  = fgetl(fid)

A = fscanf(fid,'%g', [4  3]);
ftell(fid)
% fl  = fgetl(fid)
% fl  = fgetl(fid)
% fl  = fgetl(fid)
B = fscanf(fid,'%g', [2  3]);

fclose(fid);