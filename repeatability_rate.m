
addpath 'F:\registration\Registration_baseline\RTV-SIFT\RTV';
addpath 'F:\registration\Registration_baseline\RTV-SIFT\gray_MPC_matlab';

img_indx = textread('F:\registration\Registration_baseline\img_of_osdataset.txt');

fid = fopen('F:\registration\Registration_baseline\RTV_Harris_repeatability_246.txt','a');
indexs = size(img_indx,1);
for n=1:indexs
    for error=1:10
        p = img_indx(n);
        disp(['----------------------------',num2str(p),'----------------------------']);
        im1_name = ['F:\registration\SAR-OPT-data\OSdataset\512\train\','opt',num2str(p),'.png'];
        im2_name = ['F:\registration\SAR-OPT-data\OSdataset\512\train\','sar',num2str(p),'.png'];
        r_rtv = RTV_Harris_repeatability(im1_name,im2_name,error);
        fprintf(fid,'%4f\t',[round(p),r_rtv]);
        fprintf(fid,'\r');
    end
end
fclose(fid);

function r = RTV_Harris_repeatability(im1_name,im2_name,error)
    %% image path
    image_1 = imread(im1_name);
    image_2 = imread(im2_name);
%     imwrite(image_1,['D:\Disk_G\PROJECTS\Registration_baseline\img_in_osdataset\',num2str(p),'opt.png']);
%     imwrite(image_2,['D:\Disk_G\PROJECTS\Registration_baseline\img_in_osdataset\',num2str(p),'sar.png']);
    t1=clock;
    
    %% ��ʼ�����趨
    sigma_opt=1;
    sigma_sar=2;
    sigma_2=1.6;
    ratio=2^(1/3);%�߶ȱ�
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

    uni1=position_1(:,[1,2,3,4]);
    [~,i,~]=unique(uni1,'rows','first');
    cor1=position_1(sort(i)',:);
    points_1 = cor1(:,[1,2]);
    uni2=position_2(:,[1,2,3,4]);
    [~,i,~]=unique(uni2,'rows','first');
    cor2=position_2(sort(i)',:);
    points_2 = cor2(:,[1,2]);
    Nc = 0;
    n1 = size(points_1,1);
    n2 = size(points_2,1);
    for i=1:n1
        for j=1:n2
            diff = points_2(j,:)-points_1(i,:);
            distance = sqrt(sum(diff.^2));
            if distance<error
                Nc = Nc+1;
                points_2(j,:) = [0 0];
                points_1(i,:) = [0 0];
                break;
            end
        end
    end
    r = Nc/min(n1,n2);
    disp(['repeatability rate is��',num2str(r)]);
    % showpoint_detected(image_1,image_2,position_1,position_2,1);        

%     t2=clock;
%     disp(['time is��',num2str(etime(t2,t1)),'s']);
end


                                              
                         
