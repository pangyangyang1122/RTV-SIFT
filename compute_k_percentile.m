function [k_percentile]=compute_k_percentile(gradient_x,gradient_y,perc,nbin)
%�ú�������һ���ԱȶȲ���k,����ԱȶȲ������ڼ�����ɢϵ��
%gradient_x��ˮƽ������ݶȣ�gradient_y����ֱ������ݶ�
%perc���ݶ�ֱ��ͼ�İٷ�λ����Ĭ��ȡֵ��0.7��k��ȡֵ��������ٷ�λ��ȷ��
%nbin��ֱ��ͼbin�ĸ�����Ĭ����300
%����ϵ����������k������������˶�����ͬ���ݶ�ֵ�����kֵ�ϴ��򴫵�ϵ��ֵ�ϴ�
%�����ɢ��ƽ�����أ���˿��Կ��������Ҫ����ϸ����Ҫ��С��kֵ

%���Ա߽�����ݶȵ����ֵ
gradient=sqrt(gradient_x.^2+gradient_y.^2);
[M,N]=size(gradient);
temp_gradient=gradient(2:M-1,2:N-1);
max_gradient=max(max(temp_gradient));

%���Ա߽����ֱ��ͼ
temp_gradient=temp_gradient(temp_gradient>0);
temp_gradient=floor(temp_gradient*nbin/max_gradient);
temp_gradient(temp_gradient==nbin)=nbin-1;
[M1,N1]=size(temp_gradient);
sum_pix=M1*N1;%���������ݶȸ���
hist=zeros(1,nbin);%ֱ��ͼ��ʼ��

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
for i=1:1:nbin
    nelements=nelements+hist(i);
    if(nelements>=nthreshold)
        temp_i=i;
        break;
    end
end
 
if(temp_i==0)
    k_percentile=0.03;
else
    k_percentile=max_gradient*(temp_i)/nbin;
end

end















