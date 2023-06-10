function [nonelinear_space]=create_nonlinear_scale_space(image,sigma_1,sigma_2,...
                                                        ratio,layers,nbin,perc,...
                                                        which_diff,is_auto)
%�ú������������Գ߶ȿռ�
%image�������ԭʼͼ������Ӧ���Ǹ������͵����ݣ���Χ��0-1
%sigma_1�ǵ�һ���ͼ��ĳ߶ȣ�Ĭ����1.6���߶ȿռ��һ���ͼ����image������׼��
%��sigma_1�ĸ�˹�˲��õ�
%sigma_2��ÿ�μ�����һ��ͼ��֮ǰ����֮ǰ��ͼ��ĸ�˹ƽ����׼��,Ĭ����1����
%nbin�Ǽ���Աȶ�����ʱ����Ҫ�ĳ�����Ĭ����300
%perc�Ǽ���Աȶ����ӵ��ݶȰٷ�λ������Ĭ����0.7
%ratio���������ĳ߶ȱ�
%layers�ǹ����ĳ߶ȿռ�Ĳ���������û��ʹ���²�������
%which_diff������ʹ���ĸ��������㴫��ϵ��ȡֵ��1,2,3
%nonelinear_space�ǹ����ĳ߶ȿռ�ͼ��

%%
[M,N]=size(image);
nonelinear_space=zeros(M,N,layers);

%���ȶ�����ͼ����и�˹ƽ��
windows_size=2*round(2*sigma_1)+1;
W=fspecial('gaussian',[windows_size windows_size],sigma_1);
image=imfilter(image,W,'replicate');%base_image�ĳ߶���sigma_1
nonelinear_space(:,:,1)=image;%base_image��Ϊ�߶ȿռ�ĵ�һ��ͼ��

%��ȡ�˲�������
h=[-1,0,1;-2,0,2;-1,0,1];%����˲�ģ��

%����ÿ��ĳ߶�
sigma=zeros(1,layers);
for i=1:1:layers
    sigma(i)=sigma_1*ratio^(i-1);%ÿ��ĳ߶�
end

%% ���������Գ߶ȿռ�
for i=2:1:layers
    %֮ǰ��ķ�������ɢ��ĵ�ͼ��,�����ݶ�֮ǰ����ƽ����Ŀ����Ϊ����������
    prev_image=nonelinear_space(:,:,i-1);
    windows_size=2*round(2*sigma_2)+1;
    W=fspecial('gaussian',[windows_size,windows_size],sigma_2);
    prev_smooth=imfilter(prev_image,W,'replicate');
    
    %����֮ǰ�㱻ƽ��ͼ���x��y�����һ���ݶ�
    Lx=imfilter(prev_smooth,h,'replicate');
    Ly=imfilter(prev_smooth,h','replicate');
    %ÿ�ε���ʱ����Ҫ���¶Աȶ�����k
    if(strcmp(is_auto,'NO'))
        [k_percentile]=compute_k_percentile(Lx,Ly,perc,nbin);
    elseif(strcmp(is_auto,'YES'))
        [k_percentile]=compute_k_percentile_auto(Lx,Ly,perc);
    end
    if(which_diff==1)
        [diff_c]=pm_g1(Lx,Ly,k_percentile);
    elseif(which_diff==2)
        [diff_c]=pm_g2(Lx,Ly,k_percentile);
    else
        [diff_c]=weickert_diffusivity(Lx,Ly,k_percentile);
    end
    
    %���㵱ǰ��߶�ͼ��
    step=1/2*(sigma(i)^2-sigma(i-1)^2);%��������
    nonelinear_space(:,:,i)=AOS(prev_image,step,diff_c);
end
end


%% ��ɢϵ�����㺯��1
function [g1]=pm_g1(Lx,Ly,k)
%�ú�������PM����ϵ��g1,Lx��ˮƽ����ĵ�����Ly����ֱ����ĵ���
%k��һ���Աȶ����Ӳ�����k��ȡֵһ�����ͳ������
%g1=exp(-(Lx^2+Ly^2)/k^2)

g1=exp(-(Lx.^2+Ly.^2)/k^2);

end

%% ��ɢϵ�����㺯��2
function [g2]=pm_g2(Lx,Ly,k)
%�ú�������PM���̵���ɢϵ�����ڶ��ַ���
%Lx��Ly�ֱ���ˮƽ�������ֱ����Ĳ�֣�k�ǶԱȶ����Ӳ���
%g2=1/(1+(Lx^2+Ly^2)/(k^2)),����kֵ��ȷ��һ����ͨ��ͳ�Ʒ����õ�

g2=1./(1+(Lx.^2+Ly.^2)/(k^2));
end

%% ��ɢϵ�����㺯��3
function [g3]=weickert_diffusivity(Lx,Ly,k)
%�����������weickert����ϵ��
%Lx��Ly��ˮƽ�������ֱ�����һ�ײ���ݶȣ�k�ǶԱȶ�ϵ��
%k��ȡֵһ��ͨ��ͳ�Ʒ����õ�
%g3=1-exp(-3.315/((Lx^2+Ly^2)/k^4))

g3=1-exp(-3.315./((Lx.^2+Ly.^2).^4/k^8));
end















