function [  theta_est_circular, Delay_est_conditioned, tau_est_circular_corrected]=IS_variable_generationOFDMV3(g_prime_bar,theta_vector_points,step_theta,tau_vector_points,~,Q,R,~, ~,noise_variance,rho_g,rho_Lc,M,S_time,N_FFT,H,wc,delta_f,rol,c)

P=M;
M=N_FFT;

fc=5.25e9;
t = [210 240 525]/fc;
%--------argmax



% Qhat=PathDetection(g_prime_bar,step_theta)

%--------argmax
 pdf_tau=step_theta*sum(g_prime_bar,1);
 [~,locs]=findpeaks(pdf_tau,'SORTSTR','ascend');
 max_tau=locs(end-Q+1:end);
 max_points=sort(tau_vector_points(max_tau));  
 pdf_tetha=g_prime_bar(:,sort(max_tau))./repmat(pdf_tau(sort(locs(end-Q+1:end))),length(theta_vector_points),1);
% plot(theta_vector_points*180/pi,pdf_tetha)
 tetha_max=zeros(1,Q);
 for gg=1:Q
      [~,loc]=findpeaks(pdf_tetha(:,gg),'SORTSTR','ascend');
 max_tetha=loc(end);
     tetha_max(1,gg)=theta_vector_points(max_tetha);
 end
 tetha_max=sort(tetha_max);
% tetha_max*180/pi 
%-------------------end of argmax
%------------local intervals
theta_vector_points=[];
tau_vector_points=[];
step_theta=0.01;
step_tau=0.0001;
 
tetha_max_degrees=tetha_max*180/pi;
for hh=1:Q
    theta_vector_points=[theta_vector_points tetha_max_degrees(1,hh)-2:step_theta:tetha_max_degrees(1,hh)+2];
    tau_vector_points=[tau_vector_points max_points(1,hh)-0.14*1e-8:step_tau*1e-8:max_points(1,hh)+0.14*1e-8];
 

end

 
theta_vector_pointsrad=theta_vector_points/180*pi;
N=N_FFT;
[g_prime_bar, S_frequency, X_hat]=Importance_function_gridOFDM(wc,H,S_time,theta_vector_pointsrad,step_theta,tau_vector_points,step_tau,P,N_FFT,N,noise_variance,rho_g,delta_f,rol,c);

 %-------------------------------------
%---------local pdf
 size_avoiding_interval_tau=ceil(0.14/step_tau);
size_avoiding_interval_theta=ceil(1.6/step_theta);
pdf_tau=step_theta*sum(g_prime_bar,1);
tetha_maxdeg=tetha_max_degrees;
 theta_vector_pointsdeg=theta_vector_points;
 [theta_est_matrix_conditioned,tau_est_matrix]= VariableGeneration(max_points,tetha_maxdeg,R,Q,step_theta,step_tau,g_prime_bar,tau_vector_points,theta_vector_pointsdeg,size_avoiding_interval_theta,size_avoiding_interval_tau);

 
Delay_est_conditioned=1/R*sum((tau_est_matrix),1);

% tau_est_circular = DelayWC(R,theta_est_conditioned,tau_est_matrix,tau_vector_points,noise_variance,Q,P,rho_g,rho_Lc,wc,delta_f,rol,c,N,X_hat,S_frequency)



theta_est_circular = ThetaWC(R,Delay_est_conditioned,theta_est_matrix_conditioned,theta_vector_pointsrad,noise_variance,Q,P,rho_g,rho_Lc,wc,delta_f,rol,c,N,X_hat,S_frequency);

 

%%--------- generate local delays

max_points=Delay_est_conditioned ;
size_avoiding_interval_tau=ceil(0.1/step_tau);
size_avoiding_interval_theta=ceil(0.3/step_theta);
tetha_maxdeg = theta_est_circular;

[theta_est_matrix_conditioned,tau_est_matrix]= VariableGeneration(max_points,tetha_maxdeg,R,Q,step_theta,step_tau,g_prime_bar,tau_vector_points,theta_vector_pointsdeg,size_avoiding_interval_theta,size_avoiding_interval_tau);

tau_est_circular_corrected = DelayWC(R,theta_est_circular,tau_est_matrix,tau_vector_points,noise_variance,Q,P,rho_g,rho_Lc,wc,delta_f,rol,c,N,X_hat,S_frequency);

theta_est_circular = ThetaWC(R,tau_est_circular_corrected,theta_est_matrix_conditioned,theta_vector_pointsrad,noise_variance,Q,P,rho_g,rho_Lc,wc,delta_f,rol,c,N,X_hat,S_frequency);

max_points=tau_est_circular_corrected ;
size_avoiding_interval_tau=ceil(0.02/step_tau);
size_avoiding_interval_theta=ceil(0.3/step_theta);
tetha_maxdeg = theta_est_circular;

[~,tau_est_matrix]= VariableGeneration(max_points,tetha_maxdeg,R,Q,step_theta,step_tau,g_prime_bar,tau_vector_points,theta_vector_pointsdeg,size_avoiding_interval_theta,size_avoiding_interval_tau);

tau_est_circular_corrected = DelayWC(R,theta_est_circular,tau_est_matrix,tau_vector_points,noise_variance,Q,P,rho_g,rho_Lc,wc,delta_f,rol,c,N,X_hat,S_frequency);


