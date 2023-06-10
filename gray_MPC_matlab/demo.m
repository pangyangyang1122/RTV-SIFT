
clear all;
close all;

%% 选择图片
img_path = 'D:\Disk_G\My_DataHub\testdata/SAR_pool16_PPB.png';
im = imread(img_path);

%% 常数参数设定
nscale = 4;
minWaveLength = 3;
mult = 2;
sigmaOnf = 0.55;
k = 2;
cutOff = 0.3;
deviationGain = 1.5;
noiseMethod = -1;

%% MPC提取边缘
[M or ft T] = phasecongmono(im,nscale,minWaveLength,mult,sigmaOnf,k,cutOff,10,deviationGain,noiseMethod);

%% 画图
figure;
imshow(M);
title('MPC edge');