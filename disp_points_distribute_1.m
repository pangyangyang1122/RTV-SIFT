function button=disp_points_distribute_1(locs_1,locs_2,cor2,cor1,Mmax)
%�ú�����ʾ��ʼ��⵽�ĵ���ÿ��ķֲ��������ȷ��Եĵ��ڸ���ķֲ�
%locs_1�ǲο�ͼ��㣬locs_2�Ǵ���׼ͼ���
%cor2�ǲο�ͼ����ȷ��⵽�ĵ㣬cor2�Ǵ���׼ͼ����ȷ��⵽�ĵ�
%Mmax�ǲ���

%% ��ʼ������
dis_num1=zeros(1,Mmax);%����ο�ͼ��ÿ���⵽��������
dis_num2=zeros(1,Mmax);%�������׼ͼ��ÿ���⵽��������
dis_num11=zeros(1,Mmax);%����ο�ͼ�����ÿ����ȷ��Ե�
dis_num22=zeros(1,Mmax);%�������׼ͼ�����ÿ����ȷ��Ե�
for i=1:1:Mmax
    %���м�⵽�ĵ�
    temp1=find(locs_1(:,4)==i);
    temp2=find(locs_2(:,4)==i);
    dis_num1(1,i)=size(temp1,1);
    dis_num2(1,i)=size(temp2,1);
    %��ȷ��Եĵ�
    temp1=find(cor2(:,4)==i);
    temp2=find(cor1(:,4)==i);
    dis_num11(1,i)=size(temp1,1);
    dis_num22(1,i)=size(temp2,1);
end
horz=1:1:Mmax;

%% ��ʾ��figure
%�ο�ͼ���ʼ��ֲ�
button=figure;
subplot(2,2,1);
bar(horz,dis_num1);
axis([0 Mmax+2 0 1.1*max(dis_num1)]);
xlabel('���');
ylabel('����');
title(['�ο���ʼ��ֲ�',num2str(size(locs_1,1)),'��']);
height=max(dis_num1)*0.1;
for i=1:size(horz,2)
    text(horz(1,i),dis_num1(1,i)+height,num2str(dis_num1(1,i)));
end
%����׼ͼ���ʼ��ֲ�
subplot(2,2,2);
bar(horz,dis_num2);
axis([0 Mmax+2 0 1.1*max(dis_num2)]);
xlabel('���');
ylabel('����');
title(['����׼��ʼ��ֲ�',num2str(size(locs_2,1)),'��']);
height=max(dis_num2)*0.1;
for i=1:size(horz,2)
    text(horz(1,i),dis_num2(1,i)+height,num2str(dis_num2(1,i)));
end

%�ο�ͼ����ȷ��Ե�ֲ�
subplot(2,2,3);
bar(horz,dis_num11);
xlabel('���');
ylabel('����');
axis([0 Mmax+2 0 1.1*max(dis_num11)]);
title(['�ο���ȷ��ֲ�',num2str(size(cor2,1)),'��']);
height=max(dis_num11)*0.1;
for i=1:size(horz,2)
    text(horz(1,i),dis_num11(1,i)+height,num2str(dis_num11(1,i)));
end

%����׼ͼ����ȷ��Ե�ֲ�
subplot(2,2,4);
bar(horz,dis_num22);
xlabel('���');
ylabel('����');
axis([0 Mmax+2 0 1.1*max(dis_num22)]);
title(['����׼��ȷ��ֲ�',num2str(size(cor1,1)),'��']);
height=max(dis_num22)*0.1;
for i=1:size(horz,2)
    text(horz(1,i),dis_num22(1,i)+height,num2str(dis_num22(1,i)));
end
end



