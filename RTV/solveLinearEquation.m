function OUT = solveLinearEquation(IN, wx, wy, lambda)
    [r,c,ch] = size(IN);
    k = r*c;
    dx = -lambda*wx(:); % �Ѿ�������һ��
    dy = -lambda*wy(:); % �Ѿ�������һ��
    B(:,1) = dx;
    B(:,2) = dy;
    d = [-r,-1];
    A = spdiags(B,d,k,k);

    e = dx;
    w = padarray(dx, r, 'pre'); w = w(1:end-r);
    s = dy;
    n = padarray(dy, 1, 'pre'); n = n(1:end-1);
    D = 1-(e+w+s+n); % D�Ǹ�������
 
    A = A + A' + spdiags(D, 0, k, k);  % ���A���ǣ�1+lambda*L��
    if exist('ichol','builtin')
        % �������cholesky�ֽ⣬���ÿ����㷨
        L = ichol(A,struct('michol','on'));    
        OUT = IN;
        for ii=1:ch
            tin = IN(:,:,ii);
            [tout, flag] = pcg(A, tin(:),0.1,100, L, L'); 
            OUT(:,:,ii) = reshape(tout, r, c);
        end    
    else
        % ��������ڿ����㷨����ֱ��������
        OUT = IN;
        disp('second');
        for ii=1:ch
            tin = IN(:,:,ii);
            tout = A\tin(:); % ����A����*tin(:)
            OUT(:,:,ii) = reshape(tout, r, c);
        end    
    end
        
end