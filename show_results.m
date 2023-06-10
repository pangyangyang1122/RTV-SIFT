function S = show_results(im1,im2,cor1,cor2,show_print)
    uni1=cor1(:,[1,2,3,4]);
    num = size(cor1,1);
    [~,i,~]=unique(uni1,'rows','first');
    cor1=cor1(sort(i)',:);
%     cor1_x=cor1(:,1);cor1_y=cor1(:,2);
%     f=figure;colormap('gray');imshow(im1,'border','tight','initialmagnification','fit');
%     set (gcf,'Position',[0,0,size(im1,2) size(im1,1)]);
%     axis normal;hold on;
%     scatter(cor1_x,cor1_y,60,'r','filled');hold on;
%     str1=['.\save_image\','1_opt_correct.png'];
%     saveas(f,str1,'png');
    S = Scat(im1,cor1(:,[1,2]));
    % points_num_1 = size(cor1,1);
    % x_med = median(cor1_x);
    % y_med = median(cor1_y);
    % x_med = repmat(x_med,num,1);
    % dis_x = cor1_x - x_med;
    % dis_y = cor1_y - y_med;
    % distance = sum(sqrt(dis_x.^2+dis_y.^2))/num;

    if show_print ==1
        fprintf('参考图像对应点个数是%d\n', num);
    end

    uni1=cor2(:,[1,2,3,4]);
    [~,i,~]=unique(uni1,'rows','first');
    cor2=cor2(sort(i)',:);
%     cor2_x=cor2(:,1);cor2_y=cor2(:,2);
%     f=figure;colormap('gray');imshow(im1,'border','tight','initialmagnification','fit');
%     set (gcf,'Position',[0,0,size(im1,2) size(im1,1)]);
%     axis normal;hold on;
%     scatter(cor2_x,cor2_y,60,'r','filled');hold off;
%     str1=['.\save_image\','sar_corr_distribution','.png'];
%     saveas(f,str1,'png');
end

function S = Scat(im1,a)
    [h,w] = size(im1);
    a(:,1) = a(:,1)./w;
    a(:,2) = a(:,2)./h;
    M = size(a,1);
    if M<3
        S = 0;
    else
        distances = zeros(M,M);
        for i=1:1:M
            for j=1:1:M
                dis = sqrt((a(i,1)-a(j,1))^2+(a(i,2)-a(j,2))^2);
                distances(i,j) = dis;
            end 
        end

        dis_sort = sort(distances,2);
        S = median(dis_sort(:,2:M),2);
        S = mean(S,1);
    end
end

