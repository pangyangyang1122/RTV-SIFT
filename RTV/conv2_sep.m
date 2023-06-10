function ret = conv2_sep(im, sigma)
  ksize = bitor(round(5*sigma),1);
  g = fspecial('gaussian', [1,ksize], sigma); % gaussian
  ret = conv2(im,g,'same');
  ret = conv2(ret,g','same');  
end