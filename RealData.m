
tic
clear all
close all
clc
fc=5.2e9;
wc=2*pi*fc;
c=3e8;
rol=c/fc/2;
M=2;%天线个数
B=20e6;%带宽Hz
delta_f=3.125e5;
N=B/delta_f;
N=N-8;


% K=floor(M/3);%floor((N-1)/2)-1;
%  P=floor(N/3);%floor((N)/2)-1;
 K=1;
 P=floor(N/2);%floor((N)/2)-1;
 
 
L=1;%有效路径
SNR=10;
% t=3e-8:-2.5000e-009:2e-8;%路径时延
%1/delta_f=3.2000e-006

 
load('C:\Users\souheibbenamor\Dropbox\170913\170913\csi_AP_ss_4m.mat')
tau_MAX=20;
step_tau= 0.01;
tau_vector_points=-1:step_tau:tau_MAX;
tau_vector_points=tau_vector_points*1e-6;
theta_MAX=90;
step_theta=0.1;
theta_vector_points_degrees=[-90:step_theta:theta_MAX];
theta_vector_points=pi*theta_vector_points_degrees/180;


 
   Tau_estimates={};
   theta_estimate={};
   
   
for kk=1:length(CSIraw)
HData=CSIraw{4,1}.csi;
H=squeeze(HData(:,1,:));


S_time=ones(1,K);
N_FFT=N;
noise_variance=1/10^(SNR/10);
rho_g=10;
R=4000;
rho_Lc=80000;
ct_theta= UMPfunction(H,L,K,P,M,N,delta_f,c,fc,rol);
ct_theta(:,1)
 ct_theta(:,2)

[g_prime_bar, S_frequency, X_hat]=Importance_function_gridOFDM(wc,H,S_time,theta_vector_points,step_theta,tau_vector_points,step_tau,M,N_FFT,N,noise_variance,rho_g,delta_f,rol,c);

 [ii,jj]=max(sum(g_prime_bar,1));
TempDelay=tau_vector_points(jj) 
ff
if TempDelay<0
     Tau_estimates{kk,1} = nan;
  theta_estimate{kk,1} = nan;
else
 
  [  theta_est_circular, tau_est_circular, tau_est_circular_corrected]=IS_variable_generationOFDM(g_prime_bar,theta_vector_points,step_theta,tau_vector_points,step_tau,L,R,S_frequency, X_hat,noise_variance,rho_g,rho_Lc,M,S_time,N_FFT,H,wc,delta_f,rol,c);

  Tau_estimates{kk,1} = tau_est_circular_corrected;
  theta_estimate{kk,1} = theta_est_circular;
end
Tau_estimates{kk,1} = tau_est_circular_corrected
  theta_estimate{kk,1} = theta_est_circular
  ff

end
% ff