function [RTV_scale_space] = RTV_space_sar(image_11,MaxLayer,sigma)

lambda = 0.004;
sharpness = 0.02;
dec=2.0;
sigma_iter = sigma;
x = image_11;
% figure,imshow(x)
[M, N] = size(x);
RTV_scale_space = zeros(M,N,MaxLayer);
% RTV_scale_space(:,:,1) = image_11;

for i = 1:1:MaxLayer
    [wx, wy] = computeTextureWeights(x, sigma_iter, sharpness); % 这个x是上一代的Vs
    x = solveLinearEquation(image_11, wx, wy, lambda); % 这个x是下一代的Vs
    sigma_iter = sigma_iter/dec;  % 每次迭代sigma就会缩小一半
    if sigma_iter < 0.5
    sigma_iter = 0.5;
    end
%     figure,imshow(x)
    RTV_scale_space(:,:,i) = x;
%     str=['.\save_image\','RTV_space_',num2str(i),'.png'];
%     imwrite(x,str);
end
end