function [CC] = transform_out(data,y,x) 
%transform�Ǻ����ĵ������ƣ�����Ҫ�����ǽ�MATLAB������ģ��2D��������תΪ��ΪSurfer��.grd��ʽ��
%������Ҫ����Ĳ���Ϊ��f ��ת�������ݣ�x ԭʼ���ݵĺ�������������y Ϊ��������������
%%
%��MATLAB��2ά����Ϊ����A*B�е�A��������Ҳ����������B�������Ҳ����������
%MATLAB��2ά�����ǰ���Y*X�ķ�ʽ�����о���ģ�
%Surfer��.grd�����ǰ�ĳһ���������ΰ��з����γ�����һ�����ݣ�
x=x;y=y;
M=length(x);N=length(y);
for k=1:N
BB(1+((k-1)*M):k*M)=data(k,:)';
XX(1+((k-1)*M):k*M)=x;
YY(1+((k-1)*M):k*M)=y(k);
end
BB=BB';
XX=XX';
YY=YY';
CC=[XX YY BB];
