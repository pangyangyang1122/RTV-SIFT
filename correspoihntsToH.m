clear;
close all;
% addpath 'RTV';
% addpath 'gray_MPC_matlab';
file_image=  '../finaldataset/'; %  'E:\image_registration\pooling16'
[filename_1,pathname_1]=uigetfile({'*.*','All Files(*.*)'},'Reference image',...
                          file_image);
image_1=imread(strcat(pathname_1,filename_1));
[filename_2,pathname_2]=uigetfile({'*.*','All Files(*.*)'},'Image to be registered',...
                          file_image);
image_2=imread(strcat(pathname_2,filename_2));

% cor1=load(strcat(file_image,'8_opt.txt'));
% cor2=load(strcat(file_image,'8_sar.txt'));
% [solution,rmse,cor11,cor22]=FSC(cor2,cor1,'affine',1); 
% fhand=appendimages(image_1,image_2,cor1,cor2);
currentFile = sprintf('H_%s.mat',filename_1(1));
H = load(currentFile);   %加载标准变换矩阵
H = struct2cell(H);
solution = cell2mat(H);

f_3 = image_fusion(image_1,image_2,solution);
f = figure;
imshow(f_3,'border','tight','initialmagnification','fit');
str = 'man_9';
saveas(f,str,'png');
% % 保存变换矩阵
% currentFile = sprintf('H_%s2to1.mat',filename_1(1));
% save (currentFile,'solution');
