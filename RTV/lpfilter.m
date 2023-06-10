function FBImg = lpfilter(FImg, sigma)     
    FBImg = FImg;
    for ic = 1:size(FBImg,3)
        FBImg(:,:,ic) = conv2_sep(FImg(:,:,ic), sigma);
    end   
end