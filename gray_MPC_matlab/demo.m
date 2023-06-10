
clear all;
close all;

%% ѡ��ͼƬ
img_path = 'D:\Disk_G\My_DataHub\testdata/SAR_pool16_PPB.png';
im = imread(img_path);

%% ���������趨
nscale = 4;
minWaveLength = 3;
mult = 2;
sigmaOnf = 0.55;
k = 2;
cutOff = 0.3;
deviationGain = 1.5;
noiseMethod = -1;

%% MPC��ȡ��Ե
[M or ft T] = phasecongmono(im,nscale,minWaveLength,mult,sigmaOnf,k,cutOff,10,deviationGain,noiseMethod);

%% ��ͼ
figure;
imshow(M);
title('MPC edge');