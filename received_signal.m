%%---------closely spaced reinject circular angles to correct the delay
% addpath('jade','IS')
% rmpath('jade','IS')
clc
clear all
close all
R=2000;
rho_g=10;
rho_Lc=8000000;
SNR_dB=[-10:5: 30];
Mc=500;
k=1;
P=5;% number of sensors
Q=2;
chann_coeff=[1  1];
tau_vector=k*[2.5 5 ];% these delays are multiples of Ts (the sampling period of the known waveform) in order to support fractional delys of the sampling period (for observation)
theta_vector_degrees= [75  90]; % the reflection angles in degrees
theta_vector=theta_vector_degrees/180*pi; % the same reflection angles in radians
Np=length(tau_vector);

N=1;
beta =0.25;
P1=2;
L =128/P1;

[H,g] = genspacetimechan(P,P1,L,beta,theta_vector_degrees,tau_vector,chann_coeff,N);
%   plot(g)
%    figure
%    plot(abs(fftshift(fft(H(1,:)))))
%    ff
[M11,PL] = size(H);
    if length(g) < PL, g = [g zeros(1,PL-length(g))]; end	% zero pad
    if length(g) > PL, 
	H = [H zeros(M11,length(g)-PL)]; 		% zero pad
	[M11,PL] = size(H);
    end	% zero pad
   recieved_signal_matrix_noise_free= H;
   S_time=g;
   M=PL; 
%    plot(g,'r')
%    hold on 
%    plot(H(1,:))
%----------- Jade err init

   tau_estimate_MC_SNR_matrix=zeros(length(SNR_dB),Q,Mc);
   tau_estimate_MC_SNR_matrix_corrected=zeros(length(SNR_dB),Q,Mc);
   theta_estimate_MC_SNR_matrix=zeros(length(SNR_dB),Q,Mc);
   
   tau_estimate_jade_MC_SNR_matrix=zeros(length(SNR_dB),Q,Mc);
   theta_estimate_jade_MC_SNR_matrix=zeros(length(SNR_dB),Q,Mc);
   
    tau_estimate_UMP_MC_SNR_matrix=zeros(length(SNR_dB),Q,Mc);    
   theta_estimate_UMP_MC_SNR_matrix=zeros(length(SNR_dB),Q,Mc);
   %--------------rmse JADE estimator init
rmse_tau_JADE=zeros(length(SNR_dB),2);
rmse_theta_JADE=zeros(length(SNR_dB),2);

rmse_tau_UMP=zeros(length(SNR_dB),2);
rmse_theta_UMP=zeros(length(SNR_dB),2);

   m1 = 3;	% stacking parameters in the algorithm
m2 = 1;
for snr_ind=1:length(SNR_dB)
     tau_est_err_Jade=0;
theta_est_err_Jade=0;
for jj=1:Mc
 

   
    SNR=10^(SNR_dB(snr_ind)/10);

Power_path1= sum(abs(S_time).^2)/M;
noise_variance=abs(chann_coeff(1,1))^2*Power_path1/SNR;


noise_matrix=zeros(P,M);
for ii=1:P
    noise_matrix(ii,:)=sqrt(noise_variance)*(randn(1,M)+1j*randn(1,M))/sqrt(2);
end

recieved_signal_matrix=recieved_signal_matrix_noise_free+noise_matrix;

%------------------------ define the grid points of the delays and the angles---------
tau_MAX=15;
step_tau=0.01;
tau_vector_points=0:step_tau:tau_MAX;

theta_MAX=180;
step_theta=0.1;
theta_vector_points_degrees=[0:step_theta:theta_MAX];
theta_vector_points=pi*theta_vector_points_degrees/180;
%------------------------ end of define the grid points of the delys and the angles---------

m1 = 3;	% stacking parameters in the algorithm
m2 = 1;
r=length(tau_vector);
[theta_jade,tau_jade] = jade(recieved_signal_matrix,S_time,r,P1,m1,m2); %joint angle-delay estimation 
tau_jade=tau_jade.';
theta_jade=(theta_jade).';

%------------------------ IS 
N_FFT=M; % do not forget to comment this line for Npoint FFT using the higher N_FFT

[g_prime_bar, S_frequency, X_hat]=Importance_function_grid(recieved_signal_matrix,S_time,theta_vector_points,step_theta,tau_vector_points,step_tau,P,N_FFT,M,noise_variance,rho_g);

[  theta_est_circular, tau_est_circular, tau_est_circular_corrected]=IS_variable_generation(g_prime_bar,theta_vector_points,step_theta,tau_vector_points,step_tau,length(tau_vector),R,S_frequency, X_hat,noise_variance,rho_g,rho_Lc,P,M,S_time,N_FFT,recieved_signal_matrix);
% ff
% %----------------------end of IS



tau_est_circular_corrected_edited=zeros(1,Q);
tau_est_circular_edited=zeros(1,Q);
theta_est_circular_edited=zeros(1,Q);
for qq=1:length(tau_est_circular)
   
 tau_est_circular_corrected_edited(1,qq)=tau_est_circular_corrected(1,qq);
tau_est_circular_edited(1,qq)=tau_est_circular(1,qq);
theta_est_circular_edited(1,qq)=theta_est_circular(1,qq);
  
end
  tau_estimate_MC_SNR_matrix(jj,:,snr_ind)=tau_est_circular_corrected_edited;
  tau_estimate_MC_SNR_matrix_corrected(jj,:,snr_ind)=tau_est_circular_corrected_edited;
  theta_estimate_MC_SNR_matrix(jj,:,snr_ind)=theta_est_circular_edited;


%   tau_estimate_MC_SNR_matrix(jj,:,snr_ind)=tau_est_circular;
%   tau_estimate_MC_SNR_matrix_corrected(jj,:,snr_ind)=tau_est_circular_corrected;
%   theta_estimate_MC_SNR_matrix(jj,:,snr_ind)=theta_est_circular;
 
   
   
%   tau_est_err_Jade=tau_est_err_Jade+(sort(tau_jade)-sort(tau_vector)/k).^2;
% theta_est_err_Jade=theta_est_err_Jade+(sort(theta_jade)-sort(theta_vector_degrees)).^2;
 

end
% rmse_tau_JADE(snr_ind,:)=rmse_tau_JADE(snr_ind,:)+(tau_est_err_Jade);
% rmse_theta_JADE(snr_ind,:)=rmse_theta_JADE(snr_ind,:)+(theta_est_err_Jade);
save tau_estimate_MC_SNR_matrix_4_5    tau_estimate_MC_SNR_matrix
save tau_est_circular_corrected_MC_SNR_matrix_4_5    tau_estimate_MC_SNR_matrix_corrected
save theta_estimate_MC_SNR_matrix_90_91     theta_estimate_MC_SNR_matrix
save tau_estimate_jade_MC_SNR_matrix_25_5    tau_estimate_jade_MC_SNR_matrix
save theta_estimate_jade_MC_SNR_matrix_75_90     theta_estimate_jade_MC_SNR_matrix
end
% rmse_mc_rmse_rmse_tau_JADE=sqrt(rmse_tau_JADE/Mc);
% rmse_mc_rmse_rmse_theta_JADE=sqrt(rmse_theta_JADE/Mc);
[ CRLB_delay,CRLB_AoA1,]=CRLB(SNR_dB,P);
% YMatrix1=[rmse_mc_rmse_rmse_tau_circular(:,1) rmse_tau_iterative(:,1) rmse_mc_rmse_rmse_tau_JADE(:,1) sqrt(CRLB_delay(:,1))*2/2];
% YMatrix2=[rmse_mc_rmse_rmse_tau_circular(:,2) rmse_tau_iterative(:,2) rmse_mc_rmse_rmse_tau_JADE(:,2) sqrt(CRLB_delay(:,1))*2/2];
% YMatrix3=[rmse_mc_rmse_theta_conditioned(:,1) rmse_theta_iterative(:,1) rmse_mc_rmse_rmse_theta_JADE(:,1) 180*sqrt(CRLB_AoA1(:,1))/pi];
% YMatrix4=[rmse_mc_rmse_theta_conditioned(:,2) rmse_theta_iterative(:,2) rmse_mc_rmse_rmse_theta_JADE(:,2) 180*sqrt(CRLB_AoA1(:,1))/pi];
% 
Q=2;
snr_db=-10:5:30;
Mc=500;
rmse_tau=zeros(length(snr_db),Q);
rmse_theta=zeros(length(snr_db),Q);
for ii=1:length(snr_db)
rmse_tau(ii,:)=sqrt(sum((tau_estimate_MC_SNR_matrix_corrected(1:Mc,:,ii)/2-repmat([2.5 5],Mc,1)).^2)/Mc)
end
for ii=1:length(snr_db)
rmse_theta(ii,:)=sqrt(sum((theta_estimate_MC_SNR_matrix(1:Mc,:,ii)-repmat([75 90],Mc,1)).^2)/Mc)
end

figure (5)
semilogy(-10:5:30,rmse_tau)
figure(2)
semilogy(-10:5:30,rmse_theta)


