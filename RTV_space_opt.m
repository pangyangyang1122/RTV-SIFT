function [RTV_scale_space] = RTV_space_opt(image_11,MaxLayer,sigma)

lambda = 0.004;
sharpness = 0.04;
dec=2.0;
sigma_iter = sigma;
x = image_11;
% figure,imshow(x)
[M, N] = size(x);
RTV_scale_space = zeros(M,N,MaxLayer);
% RTV_scale_space(:,:,1) = image_11;

for i = 1:1:MaxLayer
    [wx, wy] = computeTextureWeights(x, sigma_iter, sharpness); % ���x����һ����Vs
    x = solveLinearEquation(image_11, wx, wy, lambda); % ���x����һ����Vs
    x = adapthisteq(x); % ���ƶԱȶȵ�����Ӧֱ��ͼ���⻯
    % ÿ�ε���sigma�ͻ���Сһ��
    sigma_iter = sigma_iter/dec;  
    if sigma_iter < 0.5
        sigma_iter = 0.5;
    end
%     figure,imshow(x)
    
%     if i == 1
%         figure,imshow(x)
%         imwrite(x,['.\save_image\',num2str(i),'_rtv_harris_space.png']);
%     end
    RTV_scale_space(:,:,i) = x;
end
end