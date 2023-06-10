close all;
clear;
addpath 'RTV';
addpath 'gray_MPC_matlab';

for n=1:4
    disp(['--------------------------',num2str(n),'----------------------']);
    im1_name = ['F:\Registration_baseline\testdata_2\',num2str(n),'_opt.png'];
    im2_name = ['F:\Registration_baseline\testdata_2\',num2str(n),'_sar.png'];
    image_1 = imread(im1_name);
    image_2 = imread(im2_name);
    
    % range_v = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9];
    fid = fopen('F:\Registration_baseline\noise_intensity\rtv_noise_a.txt','a');

    for i=0.0:0.01:0.1
        disp(['--------------------------',num2str(i),'----------------------']);
        [N,S] = RTV_sift(image_1,image_2,i);
        fprintf(fid,'%4f\t',[N,S]);
        fprintf(fid,'\r');
%          figure, imshow(I_1);
%          str1 = ['D:\Disk_G\PROJECTS\Registration_baseline\testdata\',num2str(i),'_noise.png'];
%          imwrite(I_1,str1);
    end
    fclose(fid);
end


function [N,S] = RTV_sift(image_1,image_2,level)

    t1=clock;
    %% ת������ͼ���ʽ
    [~,~,num1]=size(image_1);
    [~,~,num2]=size(image_2);
    if(num1==3)
        image_11=rgb2gray(image_1);
    else
        image_11=image_1;
    end
    if(num2==3)
        image_22=rgb2gray(image_2);
    else
        image_22=image_2;
    end

    %ת��Ϊ��������
    image_11=im2double(image_11);
    image_22=im2double(image_22); 
    [M,N,c] = size(image_11);
    noise = 0+level*randn(M,N,c);
    image_11 = image_11 + noise;
    
    
    %% ��ʼ�����趨
    sigma_opt=2;
    sigma_sar=4;
    sigma_2=6;
    POED = 1;
    ratio=2^(-1/3);%�߶ȱ�
    Mmax=6;%�߶ȿռ�Ĳ���
    first_layer=1;%��ֵ���⿪ʼ����
    d=0.04;%HARRIS�������ⳣ��Ĭ����0.04
    d_SH_1=0.1;%�ο�ͼ����ֵ�����scharr�˲�ʱ��ȡֵ�ϴ�500�������sobel�˲�ȡֵ��С
    d_SH_2=0.1;%����׼ͼ����ֵ
    change_form='similarity';%���������Ʊ任������任��
    sift_or_log_polar='����������������';%�����ǡ����������������ӡ��͡�SIFT�����ӡ�

     

    %% �������������RTV�߶ȿռ䣬
    tic;
    [nonelinear_space_1] = RTV_space_opt(image_11,Mmax,sigma_opt);
    [nonelinear_space_2] = RTV_space_sar(image_22,Mmax,sigma_sar);
    % disp(['����RTV�߶ȿռ仨��ʱ�䣺',num2str(toc),'��']);

    %% �߶ȿռ��ÿ��ͼ���ʾharris����
    tic;
    [harris_function_1,gradient_1,angle_1,gradient_sobel_1,gradient_11]=...
        harris_scale(nonelinear_space_1,d,sigma_2,ratio);  
    [harris_function_2,gradient_2,angle_2,gradient_sobel_2,gradient_22]=...
        harris_scale(nonelinear_space_2,d,sigma_2,ratio);                                                                                          
    % disp(['����HARRIS�����߶ȿռ仨��ʱ�䣺',num2str(toc),'��']);                    

    %% ��RTV-HARRIS�����в��Ҽ�ֵ��
    tic;
    [position_1]=find_scale_extreme(harris_function_1,d_SH_1,sigma_2,ratio,...
                 gradient_1,angle_1,first_layer);
    [position_2]=find_scale_extreme(harris_function_2,d_SH_2,sigma_2,ratio,...
                 gradient_2,angle_2,first_layer);

    %% ����ο�ͼ��ʹ���׼ͼ���������
    tic;
    [descriptors_1,locs_1]=calc_descriptors(gradient_1,angle_1,...
                                            position_1,sift_or_log_polar);                                     
    [descriptors_2,locs_2]=calc_descriptors(gradient_2,angle_2,...
                                            position_2,sift_or_log_polar);   


    %% ��ʼƥ��
    tic;
    [solution,rmse,cor1,cor2]=match(descriptors_2,locs_2,...
                                    descriptors_1,locs_1,change_form,POED);
    t2=clock;
    N = size(cor1,1);
    disp(['RTV-SIFT RMSE is��',num2str(rmse)]);
    disp(['RTV-SIFT total time is��',num2str(etime(t2,t1)),'s']);

    %% ͼ���ں�
    S = show_results(image_1,image_2,cor1,cor2,0);
    disp(['RTV-SIFT S is��',num2str(S)]);
end
                                              


