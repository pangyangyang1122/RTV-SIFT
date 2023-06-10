function showpoint_detected(im1,im2,cor1,cor2,filename,show_print,j)
uni1=cor1(:,[1,2,3,4]);
[~,i,~]=unique(uni1,'rows','first');
cor1=cor1(sort(i)',:);
cor1_x=cor1(:,1);cor1_y=cor1(:,2);
name = char('detect', 'correct');
f=figure;colormap('gray');imshow(im1);title('1');
hold on;
scatter(cor1_x,cor1_y,60,'r','filled');hold on;
set(gca,'ticklength',[0 0])
set(gca,'xtick',[],'xticklabel',[]);set(gca,'ytick',[],'yticklabel',[])
set(gca,'Position',[0 0 1 1]);
str1=['.\save_image\',filename,name(j,:),'.png'];
% saveas(f,str1,'png');
% points_num_1 = size(cor1,1);
if show_print ==1
    fprintf('参考图像检测到特征点个数是%d\n', size(cor1,1));
end

uni1=cor2(:,[1,2,3,4]);
[~,i,~]=unique(uni1,'rows','first');
cor2=cor2(sort(i)',:);
cor2_x=cor2(:,1);cor2_y=cor2(:,2);
f=figure;colormap('gray');imshow(im2);title('2');
hold on;
scatter(cor2_x,cor2_y,70,'r','filled');hold on;
set(gca,'ticklength',[0 0])
set(gca,'xtick',[],'xticklabel',[]);set(gca,'ytick',[],'yticklabel',[])
set(gca,'Position',[0 0 1 1]);
str1=['.\save_image\',filename,name(j,:),'.png'];
%saveas(f,str1,'png');
% points_num_2 = size(cor2,1);
if show_print ==1
    fprintf('待配准图像检测到特征点个数是%d\n', size(cor2,1));
end


end



