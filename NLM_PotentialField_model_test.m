clc;
clear all;
close all;
%%
DATA_1=load('ground_truth.txt');DATA_2=load('Noise.txt');DATA_3=load('irregular.txt');
DATA_4=load('directional.txt');DATA_5=load('non_uniform.txt');
%
[D_1,X,Y]=transform_in(DATA_1,126,126);[D_2,X,Y]=transform_in(DATA_2,126,126);
[D_3,X,Y]=transform_in(DATA_3,126,126);[D_4,X,Y]=transform_in(DATA_4,126,126);
[D_5,X,Y]=transform_in(DATA_5,126,126);
%%
%
INPUT=D_5;
%
figure
imagesc(INPUT)
colormap('hsv')
colorbar
% INPUT=noise_new;
%
ds=3;% block size for calculate weight
Ds=[6:2:20];% search block
sigma=0.001*mean(mean(INPUT));
h=[100:100:1000].*sigma;
%
%%
for i=1:length(Ds)
    for j=1:length(h)
       N_processed(i,j)= sqrt(mean(mean((D_1-NLM_II(INPUT,Ds(i),ds,h(j))).^2)));
    end
end
%
[R,C]=find(N_processed==min(min(N_processed)));
%
OUTPUT=NLM_II(INPUT,Ds(R(1)),ds,h(C(1)));
%
sigma=0.005*mean(mean(OUTPUT));
h=[100:100:1000].*sigma;
%
%%
for i=1:length(Ds)
    for j=1:length(h)
       N_processed(i,j)= sqrt(mean(mean((D_1-NLM_II(OUTPUT,Ds(i),ds,h(j))).^2)));
    end
end
%
[R,C]=find(N_processed==min(min(N_processed)));
%
OUTPUT=NLM_II(OUTPUT,Ds(R(1)),ds,h(C(1)));
RMS_1=min(min(N_processed));
%
figure
subplot(2,2,1)
imagesc(INPUT-OUTPUT)
colormap('hsv')
colorbar
subplot(2,2,2)
% figure
% imagesc(OUTPUT)
contourf(OUTPUT,15)
colormap('hsv')
colorbar
subplot(2,2,3)
imagesc(D_1-OUTPUT)
colormap('hsv')
colorbar
subplot(2,2,4)
imagesc(INPUT)
colormap('hsv')
colorbar

% figure
% contourf(D_1,15)
% colormap('hsv')
% colorbar
% figure
% contourf(OUTPUT,15)
% colormap('hsv')
% colorbar



DT1=transform_out(OUTPUT,YY,XX);
DT2=transform_out(D_1-OUTPUT,YY,XX);
DT3=transform_out(INPUT-OUTPUT,YY,XX);

[Ds(R(1)),h(C(1))]

figure
imagesc(Ds,h,N_processed)
% shading interp
h_c=colorbar('southoutside');
hold on
% scatter(Ds(R(1)),h(C(1)),40,'markerfacecolor','r','markeredgecolor','k')
colormap('bone')
axis tight
set(gca,'FontSize',12,'fontname','times new roman')
set(gca,'Ydir','normal')
h_x=xlabel('Ds');      
h_y=ylabel('h');   
set(h_x,'fontname','times new roman','fontsize',12);
set(h_y,'fontname','times new roman','fontsize',12);
set(gcf,'unit', 'normalized', 'position',[0.1 0.1 0.5 0.5])
set(gcf,'color',[1 1 1])

%% Wavelet
% [c,s]=wavedec2(INPUT,3,'sym4');  
% OUTPUT_wavelet=wrcoef2('a',c,s,'sym4');
% %
% figure
% imagesc(OUTPUT_wavelet);          
% colormap('hsv')
% colorbar
% %
% RMS_2=sqrt(mean(mean((D_1-OUTPUT_wavelet).^2)));
%% Gaussian Blur:
sigma1=1;sigma2=1;
Lx=3;Ly=3;
w_x=-Lx:1:Lx;w_y=-Ly:1:Ly;
N_enlarged=padarray(INPUT,[Lx,Ly],'symmetric','both');
[X,Y]=size(INPUT);
%%
for i=1:2*Lx+1
    for j=1:2*Ly+1
        W(i,j)=(1/(2*pi*sigma1*sigma2))*exp(-w_x(i)^2/(2*sigma1^2))*exp(-w_y(j)^2/(2*sigma2^2));
    end
end
%%
N_W=W./sum(sum(W));
%%
for i=Lx+1:X+Lx
    for j=Ly+1:Y+Ly
        RESULT(i,j)=sum(sum(N_W.*N_enlarged(i-Lx:i+Lx,j-Ly:j+Ly)));
    end
end
%%
OUTPUT_Gaussian=RESULT(Lx+1:X+Lx,Ly+1:Y+Ly);
RMS_2=sqrt(mean(mean((D_1-OUTPUT_Gaussian).^2)));
%
figure
imagesc(OUTPUT_Gaussian)
colorbar;
colormap('hsv')
%% DCT
GRAVITY=dct2(INPUT);
[m n]=size(GRAVITY);
%
QC1=zeros(m,n);
for i=1:15
for j=1:15
    QC1(i,j)=GRAVITY(i,j); 
end
end
%
OUTPUT_DCT=idct2(QC1);
RMS_3=sqrt(mean(mean((D_1-OUTPUT_DCT).^2)));
%
figure
imagesc(OUTPUT_DCT);
colormap('hsv')
colorbar
%% SVD
[u,s,v]=svd(INPUT);
%
s1=s;s0=s;
s0(abs(s1)<10)=0;
%
OUTPUT_SVD=u*s0*v';
RMS_4=sqrt(mean(mean((D_1-OUTPUT_SVD).^2)));
figure
imagesc(OUTPUT_SVD);
colormap('hsv')
colorbar
%%
RMS=[RMS_1 RMS_2 RMS_3 RMS_4];

Dt1=transform_out(OUTPUT_Gaussian,YY,XX);
Dt2=transform_out(D_1-OUTPUT_Gaussian,YY,XX);
Dt3=transform_out(INPUT-OUTPUT_Gaussian,YY,XX);

Dt4=transform_out(OUTPUT_DCT,YY,XX);
Dt5=transform_out(D_1-OUTPUT_DCT,YY,XX);
Dt6=transform_out(INPUT-OUTPUT_DCT,YY,XX);

Dt7=transform_out(OUTPUT_SVD,YY,XX);
Dt8=transform_out(D_1-OUTPUT_SVD,YY,XX);
Dt9=transform_out(INPUT-OUTPUT_SVD,YY,XX);