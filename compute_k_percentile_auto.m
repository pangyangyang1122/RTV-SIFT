function [k_percentile]=compute_k_percentile_auto(gradient_x,gradient_y,perc)
%�ú�������һ���ԱȶȲ���k,����ԱȶȲ������ڼ�����ɢϵ��
%gradient_x��ˮƽ������ݶȣ�gradient_y����ֱ������ݶ�
%perc���ݶ�ֱ��ͼ�İٷ�λ����Ĭ��ȡֵ��0.7��k��ȡֵ��������ٷ�λ��ȷ��
%����ϵ����������k������������˶�����ͬ���ݶ�ֵ�����kֵ�ϴ��򴫵�ϵ��ֵ�ϴ�
%�����ɢ��ƽ�����أ���˿��Կ��������Ҫ����ϸ����Ҫ��С��kֵ
%�ú����Զ�����kֵ��������ﲻ��Ҫָ��bin�Ĵ�С

%ֱ��ͼ���
unit=0.005;

%���Ա߽�����ݶȵ����ֵ
gradient=sqrt(gradient_x.^2+gradient_y.^2);
[M,N]=size(gradient);
temp_gradient=gradient(2:M-1,2:N-1);

%���Ա߽����ֱ��ͼ
temp_gradient=temp_gradient(temp_gradient>0);
max_gradient=max(max(temp_gradient));
min_gradient=min(min(temp_gradient));
temp_gradient=round((temp_gradient-min_gradient)/unit);
nbin=round((max_gradient-min_gradient)/unit);
hist=zeros(1,nbin+1);
[M1,N1]=size(temp_gradient);
sum_pix=M1*N1;%���������ݶȸ���

%����ֱ��ͼ
for i=1:1:M1
    for j=1:1:N1
        hist(temp_gradient(i,j)+1)=hist(temp_gradient(i,j)+1)+1;
    end
end

%ֱ��ͼ�ٷ�λ
nthreshold=perc*sum_pix;
nelements=0;
temp_i=0;
for i=1:1:nbin+1
    nelements=nelements+hist(i);
    if(nelements>=nthreshold)
        temp_i=i;
        break;
    end
end
 
%   k_percentile=max_gradient*(temp_i)/(nbin+1);
    k_percentile=(temp_i-1)*unit+min_gradient;


end

















