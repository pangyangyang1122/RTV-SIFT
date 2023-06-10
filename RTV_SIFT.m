 clear;
close all;
addpath 'RTV';
addpath 'gray_MPC_matlab';

%�ú�������RTV���۹���RTV�߶ȿռ䣬����harris���нǵ���

%% image path
file_image=  '../finaldataset'; %  'E:\image_registration\pooling16'
[filename_1,pathname_1]=uigetfile({'*.*','All Files(*.*)'},'Reference image',...
                          file_image);
image_1=imread(strcat(pathname_1,filename_1));
[filename_2,pathname_2]=uigetfile({'*.*','All Files(*.*)'},'Image to be registered',...
                          file_image);
image_2=imread(strcat(pathname_2,filename_2));

% % ��������ͼƬ
% [M,N,c] = size(image_1);
% noise = 0+0.1*randn(M,N,c);
% image_noise = im2double(image_1) + noise;


% % ���������ڵ�
% image_fog = fog(image_1,0.09);
% image_f = uint8(image_fog * 255);
% % ������SARͼ��
% dark_sar = image_2*0.3;

% figure;
% imshow(dark_sar);
% imwrite(dark_sar,'dark_sar.png');

%% ��ʼ�����趨
sigma_opt=1;
sigma_sar=2;
POED =0;
sigma_2=5;
ratio=2^(-1/3);%�߶ȱ�
Mmax=6;%�߶ȿռ�Ĳ���
first_layer=1;%��ֵ���⿪ʼ����
d=0.04;%HARRIS�������ⳣ��Ĭ����0.04����ǰ��ĳ���k
d_SH_1=0.1;%�ο�ͼ����ֵ�����scharr�˲�ʱ��ȡֵ�ϴ�500�������sobel�˲�ȡֵ��С
d_SH_2=0.1;%����׼ͼ����ֵ
change_form='affine';%������similarity��affine��͸�ӱ任��
sift_or_log_polar='����������������';%�����ǡ����������������ӡ��͡�SIFT�����ӡ�
threshold=3; %�ж���ȷƥ�����ֵ

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

%% �������������RTV�߶ȿռ䣬
tic;
[nonelinear_space_1] = RTV_space_opt(image_11,Mmax,sigma_opt);
[nonelinear_space_2] = RTV_space_sar(image_22,Mmax,sigma_sar);
% disp(['����RTV�߶ȿռ仨��ʱ�䣺',num2str(toc),'��']);

%% �߶ȿռ��ÿ��ͼ���ʾharris����
tic;
[harris_function_1,gradient_1,angle_1]=...
    harris_scale(nonelinear_space_1,d,sigma_2,ratio);  
[harris_function_2,gradient_2,angle_2]=...
    harris_scale(nonelinear_space_2,d,sigma_2,ratio);                                                                                          
% disp(['����HARRIS�����߶ȿռ仨��ʱ�䣺',num2str(toc),'��']);                    

%% ��RTV-HARRIS�����в��Ҽ�ֵ��
tic;
[position_1]=find_scale_extreme(harris_function_1,d_SH_1,sigma_2,ratio,...
             gradient_1,angle_1,first_layer);
[position_2]=find_scale_extreme(harris_function_2,d_SH_2,sigma_2,ratio,...
             gradient_2,angle_2,first_layer);
% disp(['�߶ȿռ���Ҽ�ֵ�㻨��ʱ�䣺',num2str(toc),'��']);

%% ��ʾ��⵽�Ľǵ��λ���ڲο�ͼ��ʹ���׼ͼ����
% showpoint_detected(image_1,image_2,position_1,position_2,1);
% r12 = points_num_1/points_num_2;r21 = points_num_2/points_num_1; % ������һ������Ӧ�ķ���

%% ����ο�ͼ��ʹ���׼ͼ���������
tic;
[descriptors_1,locs_1]=calc_descriptors(gradient_1,angle_1,...
                                        position_1,sift_or_log_polar);                                     
[descriptors_2,locs_2]=calc_descriptors(gradient_2,angle_2,...
                                        position_2,sift_or_log_polar);   
% disp(['��������������ʱ�䣺',num2str(toc),'��']);
                                              
%% ��ʼƥ��  image2��Ϊ����ͼ��
tic;
[solution,rmse,cor2,cor1]=match(descriptors_2,locs_2,...
                                descriptors_1,locs_1,change_form,POED);
t2=clock; 
disp(['RTV-SIFT Total time is��',num2str(etime(t2,t1)),'s']); 
% %% �����ֹ���ѡ��ı任���� 
% fprintf('RMSE%f\n:',rmse);
% [M,N]=size(cor2);
% match2_xy=cor2(:,1:2)';
% match2_xy=[match2_xy;ones(1,M)];
% t_match2_xy=solution*match2_xy;
% match1_xy=cor1(:,1:2)';
% match1_xy=[match1_xy;ones(1,M)];
% diff_match_xy=match1_xy-t_match2_xy;
% diff_match_xy=sqrt(sum(diff_match_xy.^2));
% index_cor=find(diff_match_xy<1);
% fprintf('%d/%d\n',size(index_cor,2),size(cor1,1));
% man_cor1=cor1(index_cor,1:2);
% man_cor2=cor2(index_cor,1:2);
% 
% [solution_man,~,~,~]=FSC(man_cor2,man_cor1,change_form,1); 
% fhand=appendimages(image_1,image_2,man_cor1,man_cor2);
% f_3 = image_fusion(image_1,image_2,solution_man);
% figure(3);
% imshow(f_3,'border','tight','initialmagnification','fit');
% 
% ����任����
% currentFile = sprintf('H_%s.mat',filename_1(1));
% save (currentFile,'solution');

%% ������ȷƥ����
% currentFile = sprintf('H_%s.mat',filename_1(1));
% H = load(currentFile);   %���ر�׼�任����
% H = struct2cell(H);
% solution_1 = cell2mat(H);
% [M,N]=size(cor2);
% match2_xy=cor2(:,1:2)';
% match2_xy=[match2_xy;ones(1,M)];
% t_match2_xy=solution_1*match2_xy;
% match1_xy=cor1(:,1:2)';
% match1_xy=[match1_xy;ones(1,M)];
% diff_match_xy=match1_xy-t_match2_xy;
% diff_match_xy=sqrt(sum(diff_match_xy.^2));
% index_cor=find(diff_match_xy<threshold);
% cor_opt = cor1(index_cor, 1:2);
% cor_sar = cor2(index_cor, 1:2);
% fprintf('RTV-SIFT CMN is: %d\n',size(index_cor,2));
% fprintf('RTV-SIFT CMR is: %d\n', size(index_cor,2)/size(cor2,1));
% 
% % �����˹�RMSE
% Man_RMSE=sqrt(sum(sum(diff_match_xy.^2))/M);
% disp(['RTV-SIFT manual RMSE is��',num2str(Man_RMSE)]);

% disp([filename_1,filename_2,num2str(sigma_opt),num2str(sigma_sar)]);
disp(['RTV-SIFT RMSE is��',num2str(rmse)]);
disp(['RTV-SIFT total time is��',num2str(etime(t2,t1)),'s']);
% % showpoint_detected(image_1,image_2,locs_1,locs_2,filename_1(:,1:2),1,1);
fhand=appendimages(image_1,image_2,cor1,cor2);
str1=['.\save_image_final\',filename_1(1:2),'����ƥ����.png'];
% saveas(fhand,str1,'png');
% %% ͼ���ں�
S = show_results(image_1,image_2,cor1,cor2,0);        %����ɢ��
disp(['RTV-SIFT S is��',num2str(S)]);
f_3 = image_fusion(image_1,image_2,solution);
figure;
imshow(f_3,'border','tight','initialmagnification','fit');
str=['.\save_image_final\',filename_1(1:2),'testdata��׼�ںϺ������ͼ��.png'];
% imwrite(f_3,str);
                                              
% ��ʾÿ��������ֲ�
% button=disp_points_distribute_1(locs_1,locs_2,cor1,cor2,Mmax);
% 
% str1=['.\save_image\',filename_1(1:2),'����ֲ�','.png'];
% saveas(button,str1,'png');
% showpoint_detected(image_1,image_2,cor1,cor2,filename_1(:,1:2),0,2);     

%% ������Ӧ����Ƕ�Ӧ��ķֲ�
disp(['OPT_corr:',num2str(size(cor1,1))]);
disp(['SAR_corr:',num2str(size(cor2,1))]);
fals1=setdiff(locs_1(:,[1,2]),cor1(:,[1,2]),'row');
disp(['OPT_false:',num2str(size(fals1,1))]);
fals2=setdiff(locs_2(:,[1,2]),cor2(:,[1,2]),'row');
disp(['SAR_false:',num2str(size(fals2,1))]);
% cor1_x=cor1(:,1);cor1_y=cor1(:,2);  
% cor2_x=cor2(:,1);cor2_y=cor2(:,2);
% fals1_x=fals1(:,1);fals1_y=fals1(:,2); 
% fals2_x=fals2(:,1);fals2_y=fals2(:,2); 
% f=figure(3);colormap('gray');imshow(image_1);
% hold on;
% scatter(fals1_x,fals1_y,20,'r','filled');hold on;
% scatter(cor1_x,cor1_y,30,'g','filled');hold on;
% set(gca,'ticklength',[0 0])
% set(gca,'xtick',[],'xticklabel',[]);set(gca,'ytick',[],'yticklabel',[])
% set(gca,'Position',[0 0 1 1]);
% str1=['.\save_image_final\',filename_1,'TaF.png'];
% saveas(f,str1,'png');
% 
% f=figure(4);colormap('gray');imshow(image_2);
% hold on;
% scatter(fals2_x,fals2_y,20,'r','filled');hold on;
% scatter(cor2_x,cor2_y,30,'g','filled');hold on;
% set(gca,'ticklength',[0 0])
% set(gca,'xtick',[],'xticklabel',[]);set(gca,'ytick',[],'yticklabel',[])
% set(gca,'Position',[0 0 1 1]);
% str1=['.\save_image_final\',filename_2,'TaF.png'];
% saveas(f,str1,'png');
% 
% function Iw = fog(I,beta)
% 
%     I=double(I)/255;
% 
%     [row,col,z] = size(I);
%     landline = 0;
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Iw = I;
%     A = 0.8;
%     m = floor(row/2);
%     n = floor(col/2);
% 
%     for i=1:z
%         for j=landline+1:row
%             for l=1:col
%                 d(j,l) = 1/((j-landline)^.05 + 0.0001);
%                 d2(j,l) = d(j,l)*8;
%                 d(j,l) = -0.04*sqrt((j-m).^2+(l-n).^2) + 17;
%                 td(j,l) = exp(-beta*d(j,l));
%                 Iw(j,l,i) = I(j,l,i)*td(j,l) + A*(1-td(j,l));
%             end
%         end
%     end
% 
% end

                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
