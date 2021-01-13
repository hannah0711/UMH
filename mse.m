%  clear all
 close all
SNR=[-10:10:30]
B=80e6;
Mc=10
t = [ 3 8 ]/B;

theta=[10, 45]; 
load('tau_estimate_jade_MC_SNR_matrix_25_5.mat')
load('tau_estimate_MC_SNR_matrix_4_5.mat')
load('tau_est_circular_corrected_MC_SNR_matrix_4_5.mat')
load('theta_estimate_jade_MC_SNR_matrix_75_90.mat')
load('theta_estimate_MC_SNR_matrix_90_91.mat')
for snr_ind=1:length(SNR)
    
     
TauMse(snr_ind,:)=sqrt(sum((squeeze(tau_estimate_MC_SNR_matrix(:,:,snr_ind))-repmat(t,Mc,1,1)).^2,1)/Mc);
 TauCorrectedMse(snr_ind,:)= sqrt(sum((squeeze(tau_estimate_MC_SNR_matrix_corrected(:,:,snr_ind))-repmat(t,Mc,1,1)).^2,1)/Mc);
 ThetaMSE(snr_ind,:)=sqrt(sum((squeeze( real(theta_estimate_MC_SNR_matrix(:,:,snr_ind)))-repmat(theta,Mc,1,1)).^2,1)/Mc);


  %TauMseUMP(snr_ind,:)= sqrt(sum((squeeze(real(tau_estimate_UMP_MC_SNR_matrix(:,:,snr_ind)))-repmat(t,Mc,1,1)).^2,1)/Mc);
   %ThetaMSEUMP(snr_ind,:) =sqrt(sum((squeeze(real(theta_estimate_UMP_MC_SNR_matrix(:,:,snr_ind)))-repmat(theta,Mc,1,1)).^2,1)/Mc);
%     
%     
%     
    
end

figure 
semilogy(SNR,TauCorrectedMse*1e9)
hold on
%semilogy(SNR,TauMseUMP*1e9)

figure 
semilogy(SNR,ThetaMSE)
hold on
%semilogy(SNR,ThetaMSEUMP)


