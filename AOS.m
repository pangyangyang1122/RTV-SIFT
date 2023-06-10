function [U]=AOS(U_prev,step,diff_c)
%�ú���ʵ�ּ��Է����㷨
%U_prev��ǰһ��ĸ�������ͼ��
%step�����õĲ������ӣ��������ĳ߶���sigma(i),��ǰһ��ĳ߶���sigma(i-1)
%��˲���������1/2*(sigma(i)^2-sigma(i-1)^2)
%diff_c�Ǹ���ǰһ��߶�ͼ���˹�˲������õ�����ɢϵ����СҲ��M*N
%U�ǵ�ǰ��ĸ������Գ߶�ͼ�񣬴�С��U_prevһ��

[U1]=AOS_row(U_prev,step,diff_c);
[U2]=AOS_col(U_prev,step,diff_c);
U=1/2*(U1+U2);

end


function [U1]=AOS_row(U1_prev,step,diff_c)
%�ú������з��������ɢ
%%
[M,N]=size(U1_prev);
U1=zeros(M,N);
%����ÿһ��
for i=1:1:M
    d=U1_prev(i,:);%������
    
    %���Ǿ���ĶԽ��߲���
    a=diff_c(i,:);
    a(2:N-1)=2*a(2:N-1);
    a(1:N-1)=a(1:N-1)+diff_c(i,2:N);
    a(2:N)=a(2:N)+diff_c(i,1:N-1);
    a=-1/2*a;
    
    %���Ǿ�������
    b=diff_c(i,1:N-1)+diff_c(i,2:N);
    b=1/2*b;
    
    %���Ǿ�������,�þ���Գƣ����c=b
    c=b;
    
    %�������ԽǷ�����Ľ�
    a=1-2*step*a;
    b=-2*step*b;
    c=-2*step*c;
    x=thomas_algorith(a,b,c,d);
    U1(i,:)=x;
end

end
    
function [U2]=AOS_col(U2_prev,step,diff_c)
%�ú������з������ɢ
%%
[M,N]=size(U2_prev);
U2=zeros(M,N);
%����ÿһ��
for i=1:1:N
    d=U2_prev(:,i);%������
    %���Ǿ���ĶԽ��߲���
    a=diff_c(:,i);
    a(2:M-1)=2*a(2:M-1);
    a(1:M-1)=a(1:M-1)+diff_c(2:M,i);
    a(2:M)=a(2:M)+diff_c(1:M-1,i);
    a=-1/2*a;
    
    %���Ǿ�������
    b=diff_c(1:M-1,i)+diff_c(2:M,i);
    b=1/2*b;
    
    %���Ǿ�������,�þ���Գƣ����c=b
    c=b;
    
    %�������ԽǷ�����Ľ�
    a=1-2*step*a;
    b=-2*step*b;
    c=-2*step*c;
    x=thomas_algorith(a',b',c',d');
    U2(:,i)=x';
end

end

