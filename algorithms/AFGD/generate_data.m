% Matrix Completion Simulation data generation

clc
clear 
close all

% ### Data Generation ###

% Set the parameter first
d1 = 1000;
d2 = 1000;
r = 5;

% Generate the low rank matrix 
X_star = data_generation_ms(d1,d2,r);


save(['./Data/data_d1_' num2str(d1) '_d2_' num2str(d2) '_r_' num2str(r) '.mat'],'X_star');
