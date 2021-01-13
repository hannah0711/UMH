
tic
clear all
close all
clc
fc=5.25e9;
wc=2*pi*fc;
c=3e8;
rol=c/fc/2;
M=6;
B=80e6;%
delta_f=3.125e5;
N=B/delta_f;
N=N-11;

K=floor(M/3);%floor((N-1)/2)-1;
 P=floor(N/3);%floor((N)/2)-1;

L=4; 
SNR=10;
% t=3e-8:-2.5000e-009:2e-8;%
%1/delta_f=3.2000e-006

 
 load('a1_var.mat')
tau_MAX=100;
step_tau= 0.01;
tau_vector_points=0:step_tau:tau_MAX;
tau_vector_points=tau_vector_points*1e-9+51.766e-9;
theta_MAX=90;
step_theta=0.1;
theta_vector_points_degrees=[-90:step_theta:theta_MAX];
theta_vector_points=pi*theta_vector_points_degrees/180;


 
   Tau_estimates={};
   theta_estimate={};
   
TempDelay=zeros(4,size(ff,3)); 
DoaALL = zeros(4,size(ff,3));
for jj=1:size(ff,3) 
for kk=1:size(ff,1)
% HData=CSIraw{kk,1}.csi;
% H=squeeze(HData(:,1,:));
H1=ff(kk,:,jj,:);
H=squeeze(H1);
S_time=ones(1,K);
N_FFT=N;
noise_variance=1/10^(SNR/10);
rho_g=1;
R=100;
rho_Lc=800000000000000;
[g_prime_bar, S_frequency, X_hat]=Importance_function_gridOFDM(wc,H,S_time,theta_vector_points,step_theta,tau_vector_points,step_tau,M,N_FFT,N,noise_variance,rho_g,delta_f,rol,c);
% % fff
% ct_theta= UMPfunction(H,L,K,P,M,N,delta_f,c,fc,rol);
% ct_theta(:,1)-51.766e-9 
%   [ii,jjj]=max(sum(g_prime_bar,1));
% TempDelay(kk,jj)=tau_vector_points(jjj)-51.766e-9  ;
  [ii,kkkk]=max(sum(g_prime_bar,2));
DoaALL(kk,jj)=theta_vector_points_degrees(kkkk)  ;
% fff
% if TempDelay<0
%      Tau_estimates{kk,1} = nan;
%   theta_estimate{kk,1} = nan;
% else
%  
  [  theta_est_circular, tau_est_circular, tau_est_circular_corrected]=IS_variable_generationOFDM(g_prime_bar,theta_vector_points,step_theta,tau_vector_points,step_tau,L,R,S_frequency, X_hat,noise_variance,rho_g,rho_Lc,M,S_time,N_FFT,H,wc,delta_f,rol,c);

%   [ii,jjj]=max(sum(g_prime_bar,1));
TempDelay(kk,jj)=tau_est_circular_corrected(1) ;
%   [ii,kkkk]=max(sum(g_prime_bar,2));
DoaALL(kk,jj)=theta_est_circular(1)  ;
% end
 jj
end
end
% fff
temp=TempDelay  ;
zhatInitVect=zeros(2,size(ff,3));
zhatVect=zeros(2,size(ff,3));
for hhh=1:size(ff,3)
[zhatInit,zhat]= positionesttimation(temp,c,hhh,DoaALL)
zhatInitVect(:,hhh)=zhatInit;
zhatVect(:,hhh)=zhat;
end
bs1=[146.2 , -172.6].'*1e-2;
bs2=[841.6 , -213.2].'*1e-2;
bs3=[907.7 , 338.2].'*1e-2;
bs4=[19.1 , 333.4].'*1e-2;
MU=[442.2 , 179.4].'*1e-2;
scatter(bs1(1,:),bs1(2,:),'LineWidth',1.5)
grid on
hold on
scatter(bs2(1,:),bs2(2,:),'LineWidth',1.5)
hold on
scatter(bs3(1,:),bs3(2,:),'LineWidth',1.5)
hold on
scatter(bs4(1,:),bs4(2,:),'LineWidth',1.5)
hold on
scatter(MU(1,:),MU(2,:),'LineWidth',1.5)
hold on
scatter(zhatInitVect(1,:),zhatInitVect(2,:),'x')
hold on
scatter(zhatInitVect(1,:),zhatInitVect(2,:),'<')
% ff