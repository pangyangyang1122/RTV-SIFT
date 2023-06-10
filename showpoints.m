% This function draws the matching points in the scene image and the model image
% Input
%    im1 im2 the scene image and the model image
%    cor1 cor2 the matching points e.g cor1(x,y,label)

%��ͼ���ϸ����������������һ��ԲȦ
function [fhand_1,fhand_2]=showpoints(im1,im2,cor1,cor2)
%% ����im1�ǲο�ͼ��im2�Ǵ���׼ͼ��
%fhand_1�ǲο�ͼ������fhand_2�Ǵ���׼ͼ����

cor1_x=cor1(:,1);cor1_y=cor1(:,2);point2=cor1(:,6);
%cor_x1=loc1(point1,1);cor_y1=loc1(point1,2);
fhand_1=figure;colormap('gray');imagesc(im1);
title(['�ο�ͼ��',num2str(size(point2,1)),'������']);hold on;
scatter(cor1_x,cor1_y,'r');hold on;%scatter���������ɢ��ͼ
for i=1:size(point2,1)
text(cor1_x(i),cor1_y(i),num2str(point2(i)),'color','y');
end

cor2_x=cor2(:,1);cor2_y=cor2(:,2);
%cor_x2=loc2(point2,1);cor_y2=loc2(point2,2);
fhand_2=figure;colormap('gray');imagesc(im2);
title(['����׼ͼ��',num2str(size(point2,1)),'������']);hold on;
scatter(cor2_x,cor2_y,'r');hold on;
for i=1:size(point2,1)
text(cor2_x(i),cor2_y(i),num2str(point2(i)),'color','y');
end
