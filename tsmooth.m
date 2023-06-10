function S = tsmooth(I,lambda,sigma,sharpness,maxIter)
    % 载入参数
    if (~exist('lambda','var'))
       lambda=0.01;
    end   
    if (~exist('sigma','var'))
       sigma=3.0;
    end 
    if (~exist('sharpness','var'))
        sharpness = 0.02;
    end
    if (~exist('maxIter','var'))
       maxIter=4;
    end
    
    % 载入图像
    I = im2double(I);
    x = I;
    sigma_iter = sigma;
    lambda = lambda/2.0;
    dec=2.0;
    for iter = 1:maxIter
        [wx, wy] = computeTextureWeights(x, sigma_iter, sharpness); % 这个x是上一代的Vs
        x = solveLinearEquation(I, wx, wy, lambda); % 这个x是下一代的Vs
        % 每次迭代sigma就会缩小一半
        sigma_iter = sigma_iter/dec;  
        if sigma_iter < 0.5
            sigma_iter = 0.5;
        end 
    end
    S = x;      
end

function [retx, rety] = computeTextureWeights(fin, sigma,sharpness)

   fx = diff(fin,1,2);  % 列间一阶差分
   fx = padarray(fx, [0 1 0], 'post');  % 用0填充到矩阵最后一列
   fy = diff(fin,1,1);  % 行间一阶差分
   fy = padarray(fy, [1 0 0], 'post');  % 用0填充到矩阵最后一行
      
   vareps_s = sharpness;
   vareps = 0.001;

   wto = max(sum(sqrt(fx.^2+fy.^2),3)/size(fin,3),vareps_s).^(-1); 
   fbin = lpfilter(fin, sigma);
   gfx = diff(fbin,1,2);
   gfx = padarray(gfx, [0 1], 'post');
   gfy = diff(fbin,1,1);
   gfy = padarray(gfy, [1 0], 'post');     
   wtbx = max(sum(abs(gfx),3)/size(fin,3),vareps).^(-1); 
   wtby = max(sum(abs(gfy),3)/size(fin,3),vareps).^(-1);   
   retx = wtbx.*wto;
   rety = wtby.*wto;

   retx(:,end) = 0;
   rety(end,:) = 0;
   
end

function ret = conv2_sep(im, sigma)
  ksize = bitor(round(5*sigma),1);
  g = fspecial('gaussian', [1,ksize], sigma); % gaussian
  ret = conv2(im,g,'same');
  ret = conv2(ret,g','same');  
end

function FBImg = lpfilter(FImg, sigma)     
    FBImg = FImg;
    for ic = 1:size(FBImg,3)
        FBImg(:,:,ic) = conv2_sep(FImg(:,:,ic), sigma);
    end   
end

function OUT = solveLinearEquation(IN, wx, wy, lambda)
    [r,c,ch] = size(IN);
    k = r*c;
    dx = -lambda*wx(:); % 把矩阵拉成一行
    dy = -lambda*wy(:); % 把矩阵拉成一行
    B(:,1) = dx;
    B(:,2) = dy;
    d = [-r,-1];
    A = spdiags(B,d,k,k);

    e = dx;
    w = padarray(dx, r, 'pre'); w = w(1:end-r);
    s = dy;
    n = padarray(dy, 1, 'pre'); n = n(1:end-1);
    D = 1-(e+w+s+n); % D是个列向量
 
    A = A + A' + spdiags(D, 0, k, k);  % 这个A就是（1+lambda*L）
    if exist('ichol','builtin')
        % 如果存在cholesky分解，则用快速算法
        L = ichol(A,struct('michol','on'));    
        OUT = IN;
        for ii=1:ch
            tin = IN(:,:,ii);
            [tout, flag] = pcg(A, tin(:),0.1,100, L, L'); 
            OUT(:,:,ii) = reshape(tout, r, c);
        end    
    else
        % 如果不存在快速算法，就直接逆运算
        OUT = IN;
        disp('second');
        for ii=1:ch
            tin = IN(:,:,ii);
            tout = A\tin(:); % 这是A的逆*tin(:)
            OUT(:,:,ii) = reshape(tout, r, c);
        end    
    end
        
end