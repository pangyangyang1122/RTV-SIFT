function [solution,rmse,cor1,cor2]= match(des1,loc1,...
       des2,loc2,change_form,POED)
distRatio=0.9;
des2t = des2';
M_des1 = size(des1,1);
M_des2=size(des2,1);
Td = 0.9;
cos_dis = zeros(M_des1,M_des2);
match = zeros(1,M_des1);

%对于参考图像中的每一个点寻找和待匹配图像中的相似点
tic;
for i = 1 : M_des1
  dotprods = des1(i,:) * des2t;    
  [vals,indx] = sort(acos(dotprods));  
  cos_dis(i,:) = acos(dotprods);
%     temp_des1=des1(i,:);
%     temp_des1=repmat(temp_des1,M_des2,1);
%     diff_des1=temp_des1-des2;
%     ED_distance=sqrt(sum(diff_des1.^2,2));  
%     [vals,indx] = sort(ED_distance);
%     cos_dis(i,:)=ED_distance';
    
  if (vals(1) < distRatio * vals(2))
     match(i) = indx(1);%%match保存着des2t中的对应的特征点的索引
  else
      match(i) = 0;
  end
end
% disp(['寻找和待匹配图像中的相似点所用时间：',num2str(toc),'秒']);

%输出参考图像和待配准图像中特征点个数
tic;
% fprintf('参考图像特征描述子数目%d.\n待配准图像特征描述子数目是%d.\n', size(des2,1),size(des1,1));
num1 = sum(match > 0);%%匹配的个数
% fprintf('初始距离比Found %d matches.\n', num);
[~,point1,point2]=find(match);
%保存【x,y,尺度，layer，角度】
loc1(find(loc1(:,5)>180)',5) = loc1(find(loc1(:,5)>180)',5)-360;
loc2(find(loc2(:,5)>180)',5) = loc2(find(loc2(:,5)>180)',5)-360;
cor1=loc1(point1,[1 2 3 4 5]);
cor2=loc2(point2,[1 2 3 4 5]);
cor1=[cor1 point1'];cor2=[cor2 point2'];%point2保存着最开始的特征点编号
%% 移除重复点对
uni1=[cor1(:,[1,2]),cor2(:,[1,2])];
[~,i,~]=unique(uni1,'rows','first');   %i是排序后在原位置的元素索引值，重复值中取第一次出现该值的索引值
cor1=cor1(sort(i)',:);
cor2=cor2(sort(i)',:);
num = size(cor1,1);

% fprintf('删除重复点对后Found %d matches.\n', num);
% disp(['移除重复点对所用时间：',num2str(toc),'秒']);

%%  显示初始匹配点对正确错误匹配的分布
% [cor1_corr,cor2_corr]=show_corract_false_points(image_1,image_2,loc1,loc2,cor1,cor2,10);

%% 初始移除错误点对后使用FSC算法
tic;
[solution,rmse,cor1,cor2]=FSC(cor1,cor2,change_form,1);                                     

% disp(['FSC所用时间：',num2str(toc),'秒']);

%% scale_orien_position_joint_restriction
tic;
if POED
    [solution,cor1,cor2,rmse]=scale_orien_joint_restriction(solution,loc1,loc2,cor1,cor2,des1,des2,Td,change_form,cos_dis);
end

fprintf('RTV-SIFT MN is: %d\n', size(cor1,1));
fprintf('RTV-SIFT MR2 is: %d\n', size(cor1,1)/num);
% disp(['二次匹配所用时间：',num2str(toc),'秒']);

% %画出正确点错误点的分布
% disp(['SAR_corr:',num2str(size(cor1,1))]);
% disp(['OPT_corr:',num2str(size(cor2,1))]);
% fals1=setdiff(loc1(:,[1,2]),cor1(:,[1,2]),'row');
% disp(['SAR_false:',num2str(size(fals1,1))]);
% fals2=setdiff(loc2(:,[1,2]),cor2(:,[1,2]),'row');
% disp(['OPT_false:',num2str(size(fals2,1))]);
% cor1_x=cor1(:,1);cor1_y=cor1(:,2);  
% cor2_x=cor2(:,1);cor2_y=cor2(:,2);
% fals1_x=fals1(:,1);fals1_y=fals1(:,2); 
% fals2_x=fals2(:,1);fals2_y=fals2(:,2); 
% f=figure(1);colormap('gray');imshow(im1);
% hold on;
% scatter(fals1_x,fals1_y,20,'r','filled');hold on;
% scatter(cor1_x,cor1_y,30,'g','filled');hold on;
% set(gca,'ticklength',[0 0])
% set(gca,'xtick',[],'xticklabel',[]);set(gca,'ytick',[],'yticklabel',[])
% set(gca,'Position',[0 0 1 1]);
% str1=['.\save_image_final\',file1,'TaF.png'];
% % saveas(f,str1,'png');
% 
% f=figure(2);colormap('gray');imshow(im2);
% hold on;
% scatter(fals2_x,fals2_y,20,'r','filled');hold on;
% scatter(cor2_x,cor2_y,30,'g','filled');hold on;
% set(gca,'ticklength',[0 0])
% set(gca,'xtick',[],'xticklabel',[]);set(gca,'ytick',[],'yticklabel',[])
% set(gca,'Position',[0 0 1 1]);
% str1=['.\save_image_final\',file2,'TaF.png'];
% % saveas(f,str1,'png');


