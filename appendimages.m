function  [fhand]=appendimages(image1, image2,correspond1,correspond2)
    image1 = im2double(image1);
    image2 = im2double(image2);
    temp1=size(image1,3);
    temp2=size(image2,3);
    if(temp1==1 && temp2==3)
        image2=rgb2gray(image2);
    elseif(temp1==3 && temp2==1)
        image1=rgb2gray(image1);
    end
    % Select the image with the fewest rows and fill in enough empty rows
    %   to make it the same height as the other image.
    rows1 = size(image1,1);
    rows2 = size(image2,1);

    col1=size(image1,2);
    col2=size(image2,2);
  
    if (rows1 < rows2)
         image1(rows1+1:rows2,1:col1,:) = 0;
    elseif(rows1 >rows2)
         image2(rows2+1:rows1,1:col2,:) = 0;
    end
    % Now append both images side-by-side.
    im3 = [image1 image2]; 

    %第一种方法删除白边
    % fhand=figure('Position', [100 100 size(im3,2) size(im3,1)]);
    % imshow(im3);

    %第二中方法删除白边
    fhand=figure;
    imshow(im3,'border','tight','initialmagnification','fit');
    %title(['左侧是参考图像---配对数目',num2str(size(correspond1,1)),'---右侧是待配准图像']);
    set (gcf,'Position',[0,0,size(im3,2) size(im3,1)]);
    axis normal;

    %第三种方法删除白边
    % set(0,'CurrentFigure',fhand);
    % set(gcf,'PaperPositionMode','auto');
    % set(gca,'position',[0,0,1,1]);
    % set(gcf,'position',[1,1,size(im3,2) size(im3,1)]);

    hold on;
    cols1 = size(image1,2);
    for i = 1: size(correspond1,1)
        line([correspond1(i,1) correspond2(i,1)+cols1],[correspond1(i,2) correspond2(i,2)], 'Color', 'y','LineWidth',1);
    end

    hold off;

end






