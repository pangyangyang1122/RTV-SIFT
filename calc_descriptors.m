%%���������㼯�ϼ�������������
function [descriptors,locs]=calc_descriptors...
(   gradient,...%�߶ȿռ���ݶ�
    angle,...%�߶ȿռ�ĽǶ�
    key_point_array,....%���������
    sift_or_log_polar...
)
%% ������ʼ��
%sift����������
SIFT_DESCR_WIDTH=4;
SIFT_DESC_HIST_BINS=8;

%��������������������
LOG_POLAR_WIDTH=8;%�������������
LOG_POLAR__HIST_BINS=8;%�Ƕȷ����Ϊ8�����֣�����ÿ45��һ������

M=size(key_point_array,1);%%��ȡ���������
if(strcmp(sift_or_log_polar,'SIFT������'))
    d=SIFT_DESCR_WIDTH;
    n=SIFT_DESC_HIST_BINS;
    descriptors=zeros(M,d*d*n);%SIFT������
    temp=1;
elseif(strcmp(sift_or_log_polar,'����������������'))
    d=LOG_POLAR_WIDTH;
    n=LOG_POLAR__HIST_BINS;
    descriptors=zeros(M,(2*d+1)*n);%����������������
    temp=0;
end

%locs���������������Ϣ[x,y,�߶ȣ��㣬�Ƕȣ��ݶ�]�����locs��һ����СΪM*6�ľ���
locs=key_point_array;
for i=1:1:M
    x=key_point_array(i,1);%�������ˮƽ����
    y=key_point_array(i,2);%���������ֱ����
    scale=key_point_array(i,3);%���������ڵĳ߶�
    layer=key_point_array(i,4);%���������ڵĲ���
    main_angle=key_point_array(i,5);%�������������

    %�õ����������ڲ�������ݶȺͶ��нǶ�
    current_gradient=gradient(:,:,layer);
    current_angle=angle(:,:,layer);
    
    %% sift������
    if(temp==1)
        descriptors(i,:)=calc_sift_descriptor(current_gradient,current_angle,...
                                                x,y,scale,main_angle,d,n);
    elseif(temp==0)                                                                       
        descriptors(i,:)=calc_log_polar_descriptor(current_gradient,current_angle,...
                                                x,y,scale,main_angle,d,n);
    end
    
end


