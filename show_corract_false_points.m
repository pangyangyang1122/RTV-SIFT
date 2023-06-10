function [cor1_corr,cor2_corr]=show_corract_false_points(image_1,image_2,loc1,loc2,cor1,cor2,threshold)
%求出正确匹配的点
% 求两点间的空间欧氏距离来判断是否匹配
%cor1,cor2是经过NNDR匹配的
%loc1,loc2是原始检测出的特征点
% H = load('H_trans.mat');   %加载标准变换矩阵
% H = struct2cell(H);
% solution = cell2mat(H);
% 
% [M,N]=size(cor1);
% match1_xy=cor1(:,1:2)';
% match1_xy=[match1_xy;ones(1,M)];
% t_match1_xy=solution*match1_xy;
% match2_xy=cor2(:,1:2)';
% match2_xy=[match2_xy;ones(1,M)];
% diff_match2_xy=t_match1_xy-match2_xy;
% diff_match2_xy=sqrt(sum(diff_match2_xy.^2));
% index_in=find(diff_match2_xy<threshold);
% cor1_corr=cor1(index_in,:);
% cor1_x=cor1_corr(:,1);cor1_y=cor1_corr(:,2);    %正确匹配点的坐标点
% cor2_corr=cor2(index_in,:);
% % cor2_x=cor2_corr(:,1);cor2_y=cor2_corr(:,2);
% 
% [~,i,~]=unique(cor1_corr(:,[1,2]),'rows','first');   %i是排序后在原位置的元素索引值，重复值中取第一次出现该值的索引值
% cor1=cor1_corr(sort(i)',:);
% [~,i,~]=unique(cor2_corr(:,[1,2]),'rows','first');
% cor2=cor2_corr(sort(i)',:);
% disp(['SAR_corr_pre:',num2str(size(cor1,1))]);
% disp(['OPT_corr_pre:',num2str(size(cor2,1))]);

% 去除原始检测到的重复点,(setdiff能直接去重）
uni1=loc1(:,[1,2,3,4]);
[~,i,~]=unique(uni1,'rows','first');
loc1=loc1(sort(i)',:);
uni2=loc2(:,[1,2,3,4]);
[~,i,~]=unique(uni2,'rows','first');
loc2=loc2(sort(i)',:);


% 画出NNDR的匹配点
cor1_x=cor1(:,1);cor1_y=cor1(:,2);  
cor2_x=cor2(:,1);cor2_y=cor2(:,2);
% %找出错误匹配点
%  %C=setdiff(A,B)函数返回在向量A中却不在向量B中的元素，并且C中不包含重复元素，并且从小到大排序
fals1=setdiff(loc1(:,[1,2]),cor1(:,[1,2]),'row');  
disp(['SAR_falese_pre:',num2str(size(fals1,1))]);
fals2=setdiff(loc2(:,[1,2]),cor2(:,[1,2]),'row');
disp(['OPT_falese_pre:',num2str(size(fals2,1))]);
fals1_x=fals1(:,1);fals1_y=fals1(:,2); 
fals2_x=fals2(:,1);fals2_y=fals2(:,2); 
f=figure(1);colormap('gray');imshow(image_1);
hold on;
scatter(fals1_x,fals1_y,20,'r','filled');hold on;
scatter(cor1_x,cor1_y,30,'g','filled');hold on;
set(gca,'ticklength',[0 0])
set(gca,'xtick',[],'xticklabel',[]);set(gca,'ytick',[],'yticklabel',[])
set(gca,'Position',[0 0 1 1]);
% str1=['.\save_image\','sar_c_f_rmse10_halfEMPC','.png'];
% saveas(f,str1,'png');

f=figure(2);colormap('gray');imshow(image_2);
hold on;
scatter(fals2_x,fals2_y,20,'r','filled');hold on;
scatter(cor2_x,cor2_y,30,'g','filled');hold on;
set(gca,'ticklength',[0 0])
set(gca,'xtick',[],'xticklabel',[]);set(gca,'ytick',[],'yticklabel',[])
set(gca,'Position',[0 0 1 1]);
% str1=['.\save_image\','opt_c_f_rmse10_halfEMPC','.png'];
% saveas(f,str1,'png');


