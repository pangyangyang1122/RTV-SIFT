function f_3=image_fusion(image_1,image_2,solution)
image_1 = im2double(image_1);
image_2 = im2double(image_2);
%该函数完成简单的图像融合算法，image_1和image_2是待融合的两幅图像
%fusion_image是融合之后的图像
[M1,N1,num1]=size(image_1);
[M2,N2,num2]=size(image_2);
if(num1==3 && num2==3)
    fusion_image=zeros(3*M1,3*N1,num1);
elseif(num1==1 && num2==3)
    fusion_image=zeros(3*M1,3*N1);
    image_2=rgb2gray(image_2);
elseif(num1==3 && num2==1)
    fusion_image=zeros(3*M1,3*N1);
    image_1=rgb2gray(image_1);
elseif(num1==1 && num2==1)
    fusion_image=zeros(3*M1,3*N1);
end

%% 对参考图像和待配准图像进行融合

solution_1=[1,0,N1;0,1,M1;0,0,1];
tform=maketform('projective',solution_1');
f_1=imtransform(image_1,tform, 'XYScale',1,'XData',[1 3*N1], 'YData',[1 3*M1]);
tform=maketform('projective',(solution_1*solution)');
f_2=imtransform(image_2,tform, 'XYScale',1,'XData',[1 3*N1], 'YData',[1 3*M1]);

%% 方法二
same_index=find(f_1 & f_2);%相同的区域
index_1=find(f_1 & ~f_2);
index_2=find(~f_1 & f_2);
fusion_image(same_index)=f_1(same_index)./2+f_2(same_index)./2;
% fusion_image(same_index)=f_1(same_index);
fusion_image(index_1)=f_1(index_1);
fusion_image(index_2)=f_2(index_2);
fusion_image=uint8(fusion_image);

%删除多余的区域
left_up=(solution_1*solution)*[1,1,1]';%左上角坐标
left_down=(solution_1*solution)*[1,M2,1]';%左下角坐标
right_up=(solution_1*solution)*[N2,1,1]';%右上角坐标
right_down=(solution_1*solution)*[N2,M2,1]';%右下角坐标
X=[left_up(1),left_down(1),right_up(1),right_down(1)];
Y=[left_up(2),left_down(2),right_up(2),right_down(2)];
X_min=max(floor(min(X)),1);
X_max=min(ceil(max(X)),3*N1);
Y_min=max(floor(min(Y)),1);
Y_max=min(ceil(max(Y)),3*M1);

if(X_min>N1+1)
    X_min=N1+1;
end
if(X_max<2*N1)
    X_max=2*N1;
end
if(Y_min>M1+1)
    Y_min=M1+1;
end
if(Y_max<2*M1)
    Y_max=2*M1;
end
if(size(fusion_image,3)==1)
    fusion_image=fusion_image(Y_min:Y_max,X_min:X_max);
    f_1=f_1(Y_min:Y_max,X_min:X_max);
    f_2=f_2(Y_min:Y_max,X_min:X_max);
elseif(size(fusion_image,3)==3)
    fusion_image=fusion_image(Y_min:Y_max,X_min:X_max,:);
    f_1=f_1(Y_min:Y_max,X_min:X_max,:);
    f_2=f_2(Y_min:Y_max,X_min:X_max,:);
end

% figure;
% imshow(fusion_image);
% title('融合后的图像');
% str=['.\save_image\','配准融合后图像.jpg'];
% imwrite(fusion_image,str,'jpg');

grid_num=10;
grid_size=floor(min(size(f_1,1),size(f_1,2))/grid_num);
[f1_1,f2_2,f_3]=mosaic_map(f_1,f_2,grid_size);

end


