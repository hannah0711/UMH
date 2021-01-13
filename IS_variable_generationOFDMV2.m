function [  theta_est_circular, tau_est_circular, tau_est_circular_corrected]=IS_variable_generationOFDMV2(g_prime_bar,theta_vector_points,step_theta,tau_vector_points,~,Q,R,~, ~,noise_variance,rho_g,rho_Lc,M,S_time,N_FFT,H,wc,delta_f,rol,c,B)

P=M;
M=N_FFT;

fc=5.25e9;
B=80e6;
t = [2 3 ]/B;
%--------argmax



% Qhat=PathDetection(g_prime_bar,step_theta)

%--------argmax

%Doing marginal pdf of delay
 pdf_tau=step_theta*sum(g_prime_bar,1);
 
%finding peaks of path delay (initial estimate)
 [~,locs]=findpeaks(pdf_tau,'SORTSTR','ascend');
 
%Assigning values
 max_tau=locs(end-Q+1:end);
 max_points=sort(tau_vector_points(max_tau)); 
%Conditional pdf of angle repected to delay
 pdf_tetha=g_prime_bar(:,sort(max_tau))./repmat(pdf_tau(sort(locs(end-Q+1:end))),length(theta_vector_points),1);
% plot(theta_vector_points*180/pi,pdf_tetha)
 tetha_max=zeros(1,Q);
 %Find peaks for each conditional
 for gg=1:Q
      [~,loc]=findpeaks(pdf_tetha(:,gg),'SORTSTR','ascend');
 max_tetha=loc(end);
     tetha_max(1,gg)=theta_vector_points(max_tetha);
 end
 tetha_max=sort(tetha_max);
% tetha_max*180/pi 
%-------------------end of argmax

%Setting local intervals for pairs of DOA and delay%%%%%%%%%%%%%%%
theta_vector_points=[];
tau_vector_points=[];
step_theta=0.01;
step_tau=0.0001;

 
tetha_max_degrees=tetha_max*180/pi;
for hh=1:Q
%     theta_vector_points=[theta_vector_points tetha_max_degrees(1,hh)-2:step_theta:tetha_max_degrees(1,hh)+2];
%     tau_vector_points=[tau_vector_points max_points(1,hh)-0.6*1e-8:step_tau*1e-8:max_points(1,hh)+0.6*1e-8];
     theta_vector_points=[theta_vector_points tetha_max_degrees(1,hh)-0.2:step_theta:tetha_max_degrees(1,hh)+0.2];
    tau_vector_points=[tau_vector_points max_points(1,hh)-0.3/B:step_tau*1e-8:max_points(1,hh)+0.3/B];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
theta_vector_pointsrad=theta_vector_points/180*pi;
N=N_FFT;

%Evaluation of importance function over a local interval
[g_prime_bar, S_frequency, X_hat]=Importance_function_gridOFDM(wc,H,S_time,theta_vector_pointsrad,step_theta,tau_vector_points,step_tau,P,N_FFT,N,noise_variance,rho_g,delta_f,rol,c);

 %-------------------------------------
%---------local pdf
%  size_avoiding_interval_tau=ceil(0.5/step_tau);
% size_avoiding_interval_theta=ceil(0.2/step_theta);
 size_avoiding_interval_tau=ceil(0.25/step_tau);
size_avoiding_interval_theta=ceil(0.2/step_theta);
pdf_tau=step_theta*sum(g_prime_bar,1);
tetha_maxdeg=tetha_max_degrees;
 theta_vector_pointsdeg=theta_vector_points;
 
 %Variable generation, use g bar to generate samples 
 %Section 7
 [theta_est_matrix_conditioned,tau_est_matrix]= VariableGeneration(max_points,tetha_maxdeg,R,Q,step_theta,step_tau,g_prime_bar,tau_vector_points,theta_vector_pointsdeg,size_avoiding_interval_theta,size_avoiding_interval_tau);

%  max_points-t, linear
theta_est_conditioned=1/R*sum((theta_est_matrix_conditioned),1);

tau_est_circular = DelayWC(R,theta_est_conditioned,tau_est_matrix,tau_vector_points,noise_variance,Q,P,rho_g,rho_Lc,wc,delta_f,rol,c,N,X_hat,S_frequency);


% tau_est_circular-t
theta_est_circular = ThetaWC(R,tau_est_circular,theta_est_matrix_conditioned,theta_vector_pointsrad,noise_variance,Q,P,rho_g,rho_Lc,wc,delta_f,rol,c,N,X_hat,S_frequency);

 

%%--------- generate local delays

max_points=tau_est_circular ;
size_avoiding_interval_tau=ceil(0.0125/step_tau);
size_avoiding_interval_theta=ceil(0.01/step_theta);
tetha_maxdeg = theta_est_circular;


[theta_est_matrix_conditioned,tau_est_matrix]= VariableGeneration(max_points,tetha_maxdeg,R,Q,step_theta,step_tau,g_prime_bar,tau_vector_points,theta_vector_pointsdeg,size_avoiding_interval_theta,size_avoiding_interval_tau);

tau_est_circular_corrected = DelayWC(R,theta_est_conditioned,tau_est_matrix,tau_vector_points,noise_variance,Q,P,rho_g,rho_Lc,wc,delta_f,rol,c,N,X_hat,S_frequency);

% theta_est_circular = ThetaWC(R,tau_est_circular_corrected,theta_est_matrix_conditioned,theta_vector_pointsrad,noise_variance,Q,P,rho_g,rho_Lc,wc,delta_f,rol,c,N,X_hat,S_frequency);

% tau_est_circular_corrected-t