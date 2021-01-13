function tau_est_circular = DelayWC(R,theta_est_conditioned,tau_est_no_interpolation_matrix,tau_vector_points,noise_variance,Q,P,rho_g,rho_Lc,wc,delta_f,rol,c,N,X_hat,S_frequency)


L_c=zeros(R,1);
g=zeros(R,1);

for hh=1:R
theta_est_conditioned_rad=theta_est_conditioned*pi/180;
% theta_est_conditioned_rad=[ 90 91 ]*pi/180;
    [L_c(hh),g(hh)]=weighting_FactorOFDM(X_hat,S_frequency,theta_est_conditioned_rad,tau_est_no_interpolation_matrix(hh,:),noise_variance,Q,P,rho_g,rho_Lc,wc,delta_f,rol,c,N);
end
wFactor=exp(L_c-g-max(L_c-g));
% L_c
% g
% max(L_c-g)
tau_est_circular=max(tau_vector_points)/(2*pi)*angle(real(wFactor)'*exp(1j*2*pi*tau_est_no_interpolation_matrix/max(tau_vector_points)));
tau_est_circular(find(tau_est_circular<0)) = max(tau_vector_points)+ tau_est_circular(find(tau_est_circular<0));