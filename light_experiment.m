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
    %% 初始参数设定
    sigma_opt=2;
    sigma_sar=4;
    POED = 1;
    sigma_2=6;
    ratio=2^(-1/3);%尺度比
    Mmax=6;%尺度空间的层数
    first_layer=1;%极值点检测开始层数
    d=0.04;%HARRIS函数任意常数默认是0.04
    d_SH_1=0.1;%参考图像阈值如果是scharr滤波时候取值较大500，如果是sobel滤波取值较小
    d_SH_2=0.1;%待配准图像阈值
    change_form='similarity';%可以是相似变换，仿射变换，
    sift_or_log_polar='对数极坐标描述子';%可以是‘对数极坐标描述子’和‘SIFT描述子’

    %% 转换输入图像格式
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

    %转换为浮点数据
    image_11=im2double(image_11);
    image_22=im2double(image_22);  

    %% 这里仅仅创建了RTV尺度空间，
    tic;
    [nonelinear_space_1] = RTV_space_opt(image_11,Mmax,sigma_opt);
    [nonelinear_space_2] = RTV_space_sar(image_22,Mmax,sigma_sar);
    % disp(['构造RTV尺度空间花费时间：',num2str(toc),'秒']);

    %% 尺度空间的每层图像表示harris函数
    tic;
    [harris_function_1,gradient_1,angle_1,gradient_sobel_1,gradient_11]=...
        harris_scale(nonelinear_space_1,d,sigma_2,ratio);  
    [harris_function_2,gradient_2,angle_2,gradient_sobel_2,gradient_22]=...
        harris_scale(nonelinear_space_2,d,sigma_2,ratio);                                                                                          
    % disp(['构造HARRIS函数尺度空间花费时间：',num2str(toc),'秒']);                    

    %% 在RTV-HARRIS函数中查找极值点
    tic;
    [position_1]=find_scale_extreme(harris_function_1,d_SH_1,sigma_2,ratio,...
                 gradient_1,angle_1,first_layer);
    [position_2]=find_scale_extreme(harris_function_2,d_SH_2,sigma_2,ratio,...
                 gradient_2,angle_2,first_layer);

    %% 计算参考图像和待配准图像的描述符
    tic;
    [descriptors_1,locs_1]=calc_descriptors(gradient_1,angle_1,...
                                            position_1,sift_or_log_polar);                                     
    [descriptors_2,locs_2]=calc_descriptors(gradient_2,angle_2,...
                                            position_2,sift_or_log_polar);   


    %% 开始匹配
    tic;
    [solution,rmse,cor2,cor1]=match(descriptors_2,locs_2,...
                                    descriptors_1,locs_1,change_form,POED);
    t2=clock;
    N = size(cor1,1);
    disp(['RTV-SIFT RMSE is：',num2str(rmse)]);
    disp(['RTV-SIFT total time is：',num2str(etime(t2,t1)),'s']);

    %% 图像融合
    S = show_results(image_1,image_2,cor1,cor2,0);
    disp(['RTV-SIFT S is：',num2str(S)]);
end
                                              


