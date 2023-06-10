function [retx, rety] = computeTextureWeights(fin, sigma,sharpness)

   fx = diff(fin,1,2);  % �м�һ�ײ��
   fx = padarray(fx, [0 1 0], 'post');  % ��0��䵽�������һ��
   fy = diff(fin,1,1);  % �м�һ�ײ��
   fy = padarray(fy, [1 0 0], 'post');  % ��0��䵽�������һ��
      
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