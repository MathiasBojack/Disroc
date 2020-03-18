%% set calculation parameter

test_writeParam

%% set Material parameter

test_writeMaterial

%% run the simulation
 runDisroc(Parameter,Material,Disroc_path)

%% post-treatment of the results

% for singlefracture0
singleFracture

