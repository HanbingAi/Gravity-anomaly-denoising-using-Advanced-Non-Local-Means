function [CC] = transform_out(data,y,x) 
%transform是函数的调用名称，其主要功能是将MATLAB的理论模型2D网格数据转为成为Surfer的.grd格式。
%其中需要输入的参数为：f 待转换的数据；x 原始数据的横坐标列向量；y 为纵坐标行向量；
%%
%以MATLAB的2维矩阵为例，A*B中的A代表纵轴也就是行数，B代表横轴也就是列数；
%MATLAB的2维矩阵是按照Y*X的方式来排列矩阵的；
%Surfer的.grd数据是把某一行数据依次按列放置形成最后的一列数据；
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
