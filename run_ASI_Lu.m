clc;clear;
close all;

%% set parameters
addpath Funcs
tau_p = 0.095;
tau_h = 0.95;

%% read input image
imagepath='imgs\';
imagename='074_0030.jpg';
filename = [imagepath,imagename];
I = imread(filename);

%% run ASI-Lu
[ellipses]  = ASI_Lu(tau_p, tau_h, I);

%% draw detected ellipses
drawEllipses(ellipses', I);

%%
rmpath Funcs