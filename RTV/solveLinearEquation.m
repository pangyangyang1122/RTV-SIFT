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