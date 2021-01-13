function [g_prime_bar, S_frequency, X_hat]=Importance_function_gridOFDM(wc,recieved_signal_matrix,S_time,theta_vector_points,step_theta,tau_vector_points,step_tau,P,N_FFT,M,~,rho_g,delta_f,rol,c,B)

% N_FFT=M; % do not forget to comment this line for Npoint FFT using the higher N_FFT

N=M;
% M
S_frequency=S_time; % compute the FFT of the known transmitted signal
X_hat=recieved_signal_matrix;% verify the entries of the received matrix!!!!!!!
% size(X_hat)
% P
% N
X_hat_vector=reshape(X_hat,P*N,1);
 


V_tau= exp(-1j*kron((2*pi*([1:N]-floor(N/2)-1)*delta_f),ones(1,P)).'*tau_vector_points);
% V_theta=exp( 1j*sin(theta_vector_points.')*kron((wc+2*pi*([1:N]-floor(N/2)-1)*delta_f),[0:P-1])*rol/c);
V_theta=exp( -1j*sin(theta_vector_points.')*kron((wc*ones(1,N)),[0:P-1])*rol/c);


%the same I (periodigram) in the paper, formula 39.
I_telde_matrix=V_theta*diag(conj (X_hat_vector))*V_tau;
I_matrix=abs(I_telde_matrix).^2/max(max(abs(I_telde_matrix).^2));% try the difference

% rho_prime=rho_g/(noise_variance*P*norm(S_frequency)^2);
rho_prime=rho_g;% do not forget to consider the effect of noise variance and the signal energy by using rho_prime
I_matrix_exponential=exp(rho_prime*I_matrix);

%Formula 40
g_prime_bar=I_matrix_exponential/(step_tau*step_theta*sum(sum(I_matrix_exponential)));
% 
% figure
% 
% plot(tau_vector_points*B,sum(g_prime_bar,1));
% % % hold on
% % % grid on
% figure
% %  hold on
% plot(theta_vector_points*180/pi,sum(g_prime_bar,2));
% hold on
% grid on
% % figure
% % plot(theta_vector_points*180/pi,sum(g_prime_bar,2));
% 
% [ii,jj]=max(sum(g_prime_bar,1));
% ffff=tau_vector_points(jj)-3.2000e-06
% ffff*3e8 
% ff








