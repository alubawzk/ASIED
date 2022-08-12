clc;clear;
close all;

%%
imgPath = 'img\test.jpg';
I = imread(imgPath);
elps = ASIED(I);
drawEllipses(elps',I);
