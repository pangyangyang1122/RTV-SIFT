function [end_solution,point_1,point_2,rmse]=...
second_restriction(solution,loc1,loc2,im1,im2,des1,des2,Td,change_form,distance)
                                          
loc1=loc1(:,[1 2 3 4 5]);
loc2=loc2(:,[1 2 3 4 5]);

M_des1=size(des1,1);
M_des2=size(des2,1);
    
match=zeros(M_des1,2);
for k=1:1:M_des1
    %position error
    temp1_xy=loc1(k,[1,2]);
    temp1_xy=repmat(temp1_xy,M_des2,1);
    temp1_xy=temp1_xy';
    temp1_xy=[temp1_xy;ones(1,M_des2)];
    T_temp1_xy=solution*temp1_xy;

    temp2_xy=loc2(:,[1,2]);
    temp2_xy=temp2_xy';
    temp2_xy=[temp2_xy;ones(1,M_des2)];
    diff_xy=T_temp1_xy-temp2_xy;
    T_error=sqrt(sum(diff_xy.^2,1));
    % [ind_r,ind_c,v] = find(T_error<10);
    % N_indx = size(ind_c,2);
    % fprintf(' %d matches.\n', N_indx);
%     if N_indx==0
%         continue
%     end
    
    %%Scale ratio error
%     temp_scale=loc1(k,3);
%     temp_scale=repmat(temp_scale,N_indx,1);
%     ratio_scale=loc2(indx,3)./temp_scale;
%     error_S=abs(1-ratio_scale);
%     error_S=error_S';
    
    JD=(1+T_error).*distance(k,:);

    [vals,index] = sort(JD); 
%     if N_indx==1
%         match(k,1)=index(1);
%         Dk=vals(1);
%         match(k,2)=Dk;
    if(vals(1)<vals(2)*Td)
        match(k,1)=index(1);
        Dk=vals(1);
        match(k,2)=Dk;
    else
        match(k,1)=0;
        match(k,2)=0;
    end    
end

temp_match=match(:,1)';
num = sum(temp_match > 0);
fprintf('二次配准初始距离比Found %d matches.\n', num);

[~,point1,point2]=find(temp_match);
loc11=loc1(point1,[1,2,3,4,5]);
loc22=loc2(point2,[1,2,3,4,5]);
loc11=[loc11,point2'];
loc22=[loc22,point2'];

%% logic filter
% [loc11,loc22,]=logic_filter(loc11,loc22,delat_x,delat_y,unit_x,unit_y,S,RMO);
% fprintf('logic filter后Found %d matches.\n', size(loc11,1));

uni1=[loc11(:,[1 2]),loc22(:,[1 2])];
[~,i,~]=unique(uni1,'rows','first');
loc11=loc11(sort(i)',:);loc22=loc22(sort(i)',:);
appendimages(im2,im1,loc22,loc11);

[end_solution,rmse,loc11,loc22]=FSC(loc11,loc22,change_form,1);%%Sub-pixel accuracy
%fprintf('尺度，位置，方向联合约束Found %d matches.\n', size(loc11,1));

point_1=loc11;
point_2=loc22;

end













