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
    
    fid = fopen('F:\Registration_baseline\fog_intensity\rtv_fog.txt','a');

    for i=0:0.004:0.04
        disp(['--------------------------',num2str(i),'----------------------']);
        [N,S] = RTV_sift(image_1,image_2,i);
        fprintf(fid,'%4f\t',[N,S]);
        fprintf(fid,'\r');
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
    image_11 = fog(image_11,level);
    image_11=im2double(image_11);
    image_22=im2double(image_22); 
    
    
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
    [harris_function_1,gradient_1,angle_1,~,~]=...
        harris_scale(nonelinear_space_1,d,sigma_2,ratio);  
    [harris_function_2,gradient_2,angle_2,~,~]=...
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
                                              
function Iw = fog(I,beta)

    I=double(I)/255;

    [row,col,z] = size(I);
    landline = 0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Iw = I;
    A = 0.8;
    m = floor(row/2);
    n = floor(col/2);

    for i=1:z
        for j=landline+1:row
            for l=1:col
                d(j,l) = 1/((j-landline)^.05 + 0.0001);
                d2(j,l) = d(j,l)*8;
                d(j,l) = -0.04*sqrt((j-m).^2+(l-n).^2) + 17;
                td(j,l) = exp(-beta*d(j,l));
                Iw(j,l,i) = I(j,l,i)*td(j,l) + A*(1-td(j,l));
            end
        end
    end

end

