clear all
close all
clc

%% Install noiselets
cd('./utility');
mex('realnoiselet.c');
cd('..');

cd('./PROPACK');
install_mex;
cd('..');