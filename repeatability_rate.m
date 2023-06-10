
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
    
    %% 初始参数设定
    sigma_opt=1;
    sigma_sar=2;
    sigma_2=1.6;
    ratio=2^(1/3);%尺度比
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
    disp(['repeatability rate is：',num2str(r)]);
    % showpoint_detected(image_1,image_2,position_1,position_2,1);        

%     t2=clock;
%     disp(['time is：',num2str(etime(t2,t1)),'s']);
end


                                              
                         
