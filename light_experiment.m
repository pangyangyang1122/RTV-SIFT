addpath 'RTV';
addpath 'gray_MPC_matlab';

for n=1:4
    disp(['--------------------------',num2str(n),'----------------------']);
    im1_name = ['F:\Registration_baseline\testdata_2\',num2str(n),'_opt.png'];
    im2_name = ['F:\Registration_baseline\testdata_2\',num2str(n),'_sar.png'];
    image_1 = imread(im1_name);
    image_2 = imread(im2_name);

    range_v = [0.7,0.76,0.82,0.89,0.94,1.0,1.12,1.24,1.36,1.48,1.6];
    fid = fopen('F:\Registration_baseline\Light_intensity_variation\rtv_ray.txt','a');

    for i=1:11
        disp(['--------------------------',num2str(i),'----------------------']);
%         I_hsv = rgb2hsv(image_1);
%         I_hsv(:,:,3) = range_v(i)*I_hsv(:,:,3);
%         I_1 = hsv2rgb(I_hsv);
        I_1 = range_v(i)*image_1;
        [rmse,N,S] = RTV_sift(I_1,image_2);
        fprintf(fid,'%4f\t',[rmse,N,S]);
        fprintf(fid,'\r');
         % figure, imshow(I_1);
    end
    fclose(fid);
end


function [rmse,N,S] = RTV_sift(image_1,image_2)

    t1=clock;
    %% ��ʼ�����趨
    sigma_opt=2;
    sigma_sar=4;
    POED = 1;
    sigma_2=6;
    ratio=2^(-1/3);%�߶ȱ�
    Mmax=6;%�߶ȿռ�Ĳ���
    first_layer=1;%��ֵ���⿪ʼ����
    d=0.04;%HARRIS�������ⳣ��Ĭ����0.04
    d_SH_1=0.1;%�ο�ͼ����ֵ�����scharr�˲�ʱ��ȡֵ�ϴ�500�������sobel�˲�ȡֵ��С
    d_SH_2=0.1;%����׼ͼ����ֵ
    change_form='similarity';%���������Ʊ任������任��
    sift_or_log_polar='����������������';%�����ǡ����������������ӡ��͡�SIFT�����ӡ�

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
    [solution,rmse,cor2,cor1]=match(descriptors_2,locs_2,...
                                    descriptors_1,locs_1,change_form,POED);
    t2=clock;
    N = size(cor1,1);
    disp(['RTV-SIFT RMSE is��',num2str(rmse)]);
    disp(['RTV-SIFT total time is��',num2str(etime(t2,t1)),'s']);

    %% ͼ���ں�
    S = show_results(image_1,image_2,cor1,cor2,0);
    disp(['RTV-SIFT S is��',num2str(S)]);
end
                                              


