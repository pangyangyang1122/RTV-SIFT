function [key_point_array]=find_scale_extreme...
    (sar_harris_function,...%�����߶ȼ���õ���SAR-HARRIS����
    threshold,...%�ǵ�����ֵ��Ĭ����0.1
    sigma,...%��һ��ĳ߶ȣ�Ĭ����2
    ratio,...%��������ĳ߶ȱ�
    gradient,...%�����ݶ�
    angle,...%����ĽǶ�
    first_layer)

%% �ú����������Ǹ������ɵ�SAR-Harris�߶ȿռ��SAR-Harris����Ѱ�Ҿֲ���ֵ��
%������˫���Բ�ֵ��ȷ��λ��ֵ��λ�ã��������sar_harris_function��sar_harris_function
%�߶ȿռ䣬key_point_array��һ���ṹ�����飬�������������Ϣ��
%����threshold��һ������SAR-HARRIS����ֵ
%sigma�ǵײ�߶�,ratio�ǳ߶ȱ���
%����harris_fun�Ǹ����Harris����

%% ��ʼ 
[M,N,num]=size(sar_harris_function);
BORDER_WIDTH=2;%�߽糣��
HIST_BIN=36;%ֱ��ͼ����������36��ֱ��ͼ��ÿ10��һ��
SIFT_ORI_PEAK_RATIO=0.9;
key_number=0;%�������������
%key_point_array=zeros(M*N*num,6);
key_point_array=zeros(M,6);
%% key_point_array��һ����ά�������ڱ���������������Ϣ������λ��x,y,���ڵĲ��
% �߶�scale,���ڵĲ�layer,������������ĽǶ�angle���ۼӵ��ݶ�gradient
for i=first_layer:1:num
    if i==1
        temp_current=sar_harris_function(:,:,i);
        temp_prev=sar_harris_function(:,:,i+1);
        temp_next=sar_harris_function(:,:,i+2);
    else
        if i==6
        temp_current=sar_harris_function(:,:,i);
        temp_prev=sar_harris_function(:,:,i-1);
        temp_next=sar_harris_function(:,:,i-2);
        else
            temp_prev=sar_harris_function(:,:,i-1);
            temp_next=sar_harris_function(:,:,i+1);
            temp_current=sar_harris_function(:,:,i);  % ��������˱������������ԱȵĽ������ʾ��harris����������������֮��ȡ�ü���ֵ
        end
    end
    gradient_current=gradient(:,:,i);%������������ڲ���ݶ�
    angle_current=angle(:,:,i);%������������ٲ�ĽǶ�
    
    for j=BORDER_WIDTH:1:M-BORDER_WIDTH%��
        for k=BORDER_WIDTH:1:N-BORDER_WIDTH%��
            temp=temp_current(j,k);
            if(temp>threshold &&...
                temp>temp_current(j-1,k-1) && temp>temp_current(j-1,k) && temp>temp_current(j-1,k+1) &&...
                temp>temp_current(j,k-1) && temp>temp_current(j,k+1) &&...
                temp>temp_current(j+1,k-1) && temp>temp_current(j+1,k) && temp>temp_current(j+1,k+1)&&...
                temp>temp_prev(j-1,k-1) && temp>temp_prev(j-1,k) && temp>temp_prev(j-1,k+1) &&...
                temp>temp_prev(j,k-1) && temp>temp_prev(j,k+1) &&...
                temp>temp_prev(j+1,k-1) && temp>temp_prev(j+1,k) && temp>temp_prev(j+1,k+1)&&...
                temp>temp_next(j-1,k-1) && temp>temp_next(j-1,k) && temp>temp_next(j-1,k+1) &&...
                temp>temp_next(j,k-1) && temp>temp_next(j,k+1) &&...
                temp>temp_next(j+1,k-1) && temp>temp_next(j+1,k) && temp>temp_next(j+1,k+1))
                
               scale=sigma*ratio^(i-1);
                %% ������������ֱ��ͼ�����ֵ���򣬼�������
                [hist,max_value]=calculate_oritation_hist(k,j,scale,...
                        gradient_current,angle_current,HIST_BIN);                
                mag_thr=max_value*SIFT_ORI_PEAK_RATIO; %%���������С   
                for kk=1:1:HIST_BIN %k�ǵ�ǰֱ��ͼ������
                    if(kk==1)
                        k1=HIST_BIN;
                    else
                        k1=kk-1;
                    end 
                    if(kk==HIST_BIN)%k2�ǵ�ǰֱ��ͼ�ұߵ�����
                        k2=1;
                    else
                        k2=kk+1;
                    end
                 %%��ǰ��ֵ����ǰ����ֵ�����Ҵ��������0.8��
                    if(hist(kk)>hist(k1) && hist(kk)>hist(k2)...
                         && hist(kk)>mag_thr)
                        bin=kk-1+0.5*(hist(k1)-hist(k2))/(hist(k1)+hist(k2)-2*hist(kk));
                        if(bin<0)
                            bin=HIST_BIN+bin;
                        elseif(bin>=HIST_BIN)
                            bin=bin-HIST_BIN;
                        end
                        %����������
                        key_number=key_number+1;
                        key_point_array(key_number,1)=k;%�����������꣬����x
                        key_point_array(key_number,2)=j;%�����������꣬����y
                        key_point_array(key_number,3)=sigma*ratio^(i-1);%���ڲ�ĳ߶�
                        key_point_array(key_number,4)=i;%���ڲ�
                        key_point_array(key_number,5)=(360/HIST_BIN)*bin;%0-360��
                        key_point_array(key_number,6)=hist(kk);%�ݶ�
                    end
                end
                
            end
        end
    end
end
key_point_array=key_point_array(1:key_number,:);
end

