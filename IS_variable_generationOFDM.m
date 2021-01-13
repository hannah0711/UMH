function [  theta_est_circular, tau_est_circular, tau_est_circular_corrected]=IS_variable_generationOFDM(g_prime_bar,theta_vector_points,step_theta,tau_vector_points,~,Q,R,~, ~,noise_variance,rho_g,rho_Lc,M,S_time,N_FFT,H,wc,delta_f,rol,c)

P=M;
M=N_FFT;


% plot(pdf_tau)
%--------argmax
 pdf_tau=step_theta*sum(g_prime_bar,1);
%  figure(1)
%     plot(tau_vector_points,pdf_tau)

%  [vall,locs]=findpeaks(pdf_tau,'SORTSTR','ascend');
%  energie_peak=abs(vall).^2;
%  ener=zeros(length(energie_peak),1);
%   loc_var=0;
%   eta=1-1e-4;
%   Q=1;
%  for iii=length(energie_peak):-1:1
%     
%      loc_var=loc_var+(energie_peak(iii))/sum(energie_peak);
%      if (loc_var<eta)
%          Q=Q+1;
%      end
%   ener(length(energie_peak)-iii+1,1)=loc_var; 
%  end
%  Q
% plot(tau_vector_points,pdf_tau)
% figure
%--------argmax
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
normalizatioFactor=1e-7;
tetha_max_degrees=tetha_max*180/pi;
for hh=1:Q
    theta_vector_points=[theta_vector_points tetha_max_degrees(1,hh)-1.2:step_theta:tetha_max_degrees(1,hh)+1.2];
    tau_vector_points=[tau_vector_points max_points(1,hh)-0.15*1e-8:step_tau*1e-8:max_points(1,hh)+0.15*1e-8];
%     tau_vector_points=[tau_vector_points max_points(1,hh)-0.2:step_tau:max_points(1,hh)];%%  check this

end

 
theta_vector_points=theta_vector_points/180*pi;
N=N_FFT;
[g_prime_bar, S_frequency, X_hat]=Importance_function_gridOFDM(wc,H,S_time,theta_vector_points,step_theta,tau_vector_points,step_tau,P,N_FFT,N,noise_variance,rho_g,delta_f,rol,c);

% [g_prime_bar, S_frequency, X_hat]=Importance_function_grid(recieved_signal_matrix,S_time,theta_vector_points,step_theta,tau_vector_points,step_tau,P,N_FFT,M,noise_variance,rho_g);
%-------------------------------------
%---------local pdf
 size_avoiding_interval_tau=ceil(0.12/step_tau);
size_avoiding_interval_theta=ceil(1/step_theta);
pdf_tau=step_theta*sum(g_prime_bar,1);

% plot(pdf_tau)

%-----------------------
cdf_tau=step_tau*[0 cumsum(pdf_tau(1:end-1))];% try the other way when the integral is computed using the second points of the discretization intervals !!!
  
tau_est_matrix=zeros(R,Q);
theta_est_matrix_conditioned=zeros(R,Q);
tau_est_no_interpolation_matrix=zeros(R,Q);
L_c=zeros(R,1);
g=zeros(R,1);


for rr=1:R
 tau_est_no_interpolation=zeros(1,Q);
 tau_est=zeros(1,Q);
 theta_est_conditioned=zeros(1,Q);
for qq=1:Q
%     max_points
%     tau_vector_points.'
   [~, min_pos1]=min(abs(tau_vector_points-max_points(1,qq)));
%    min_pos1
%    size_avoiding_interval_tau
%    size(cdf_tau)
%    test1=cdf_tau(min_pos1+size_avoiding_interval_tau)-cdf_tau(min_pos1-size_avoiding_interval_tau)
%    test2=cdf_tau(min_pos1-size_avoiding_interval_tau)
   u_tau=rand*(cdf_tau(min_pos1+size_avoiding_interval_tau)-cdf_tau(min_pos1-size_avoiding_interval_tau))+cdf_tau(min_pos1-size_avoiding_interval_tau);
   [~,min_pos]=min(abs(cdf_tau-u_tau));
   tau_est_no_interpolation(qq)=tau_vector_points(min_pos);
    if (max(cdf_tau(1:min_pos)-u_tau)<0)
       a=(cdf_tau(min_pos+1)-cdf_tau(min_pos))/step_tau;
       b=cdf_tau(min_pos)-a*tau_vector_points(min_pos);
       tau_est(qq)=(u_tau-b)/a;
    else
       a=(cdf_tau(min_pos)-cdf_tau(min_pos-1))/step_tau;
       b=cdf_tau(min_pos-1)-a*tau_vector_points(min_pos-1);
       tau_est(qq)=(u_tau-b)/a;
    end % endif

  
    
    % --------------- generating the angle with the condionned pdf-----------
    pdf_theta_conditioned=g_prime_bar(:,min_pos)/pdf_tau(min_pos);
    cdf_theta_conditioned=step_theta*[0 cumsum(pdf_theta_conditioned(1:end-1)).'];% try the other way when the integral is computed using the second points of the discretization intervals.. try also the wrong way !!!
     [~, min_pos3]=min(abs(theta_vector_points-tetha_max(1,qq)));

u_theta_conditioned=rand*(cdf_theta_conditioned(min_pos3+size_avoiding_interval_theta)-cdf_theta_conditioned(min_pos3-size_avoiding_interval_theta))+cdf_theta_conditioned(min_pos3-size_avoiding_interval_theta);

    [~,min_pos]=min(abs(cdf_theta_conditioned-u_theta_conditioned));
    theta_est_conditioned(qq)=theta_vector_points(min_pos);

      % --------------- end of generating the angle with the condionned pdf-----------

     
end % end qq

tau_est_matrix(rr,:)=(tau_est);
tau_est_no_interpolation_matrix(rr,:)=(tau_est_no_interpolation);
theta_est_matrix_conditioned(rr,:)=(theta_est_conditioned)*180/pi;

% [L_c(rr),g(rr)]=weighting_Factor(X_hat,S_frequency,theta_est_conditioned,tau_est,noise_variance,P,Q,M,rho_g,rho_Lc);

end % end rr
theta_est_conditioned=1/R*sum((theta_est_matrix_conditioned),1);

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
% wFactor
% 
% tau_est_no_interpolation_matrix
% 
% tau_est_circular
% ff





%%------- circular doa using circular tau----
for hh=1:R
theta_est_matrix_conditioned_rad=theta_est_matrix_conditioned*pi/180;
    [L_c(hh),g(hh)]=weighting_FactorOFDM(X_hat,S_frequency,theta_est_matrix_conditioned_rad(hh,:),tau_est_circular,noise_variance,Q,P,rho_g,rho_Lc,wc,delta_f,rol,c,N);
end
% theta_est_matrix_conditioned
wFactor=exp(L_c-g-max(L_c-g));
theta_vector_points=theta_vector_points*180/pi;
% theta_est_circular=max(theta_vector_points)/(2*pi)*angle(real(wFactor)'*exp(1j*2*pi*theta_est_matrix_conditioned/max(theta_vector_points)));
% theta_est_circular(find(theta_est_circular<0)) = max(theta_vector_points)+ theta_est_circular(find(theta_est_circular<0));
theta_est_circular=max(abs(theta_vector_points))*angle(real(wFactor)'*exp(1j*theta_est_matrix_conditioned/max(abs(theta_vector_points))));
% theta_est_circular(find(theta_est_circular<0)) = max(theta_vector_points)+ theta_est_circular(find(theta_est_circular<0));
% tau_est_circular
% theta_est_conditioned
% theta_est_matrix_conditioned
% tetha_max_degrees
% theta_est_circular
% wFactor
% ff
% tau_est_circular
%%---------- end of  circular doa using circular tau---

%%--------------circular delay using circular doa-----

%%--------- generate local delays

max_points=tau_est_circular ;
 size_avoiding_interval_tau=ceil(0.02/step_tau);

for rr=1:R
 tau_est_no_interpolation=zeros(1,Q);
 tau_est=zeros(1,Q);
for qq=1:Q

   [~, min_pos1]=min(abs(tau_vector_points-max_points(1,qq)));
%       min_pos1
%    size_avoiding_interval_tau
%    min_pos1-size_avoiding_interval_tau
%    size(cdf_tau)
%    test1=cdf_tau(min_pos1+size_avoiding_interval_tau)
%    cdf_tau(min_pos1-size_avoiding_interval_tau)
%    test2=cdf_tau(min_pos1-size_avoiding_interval_tau)
   u_tau=rand*(cdf_tau(min_pos1+size_avoiding_interval_tau)-cdf_tau(min_pos1-size_avoiding_interval_tau))+cdf_tau(min_pos1-size_avoiding_interval_tau);
   [~,min_pos]=min(abs(cdf_tau-u_tau));
   tau_est_no_interpolation(qq)=tau_vector_points(min_pos);
    if (max(cdf_tau(1:min_pos)-u_tau)<0)
       a=(cdf_tau(min_pos+1)-cdf_tau(min_pos))/step_tau;
       b=cdf_tau(min_pos)-a*tau_vector_points(min_pos);
       tau_est(qq)=(u_tau-b)/a;
    else
       a=(cdf_tau(min_pos)-cdf_tau(min_pos-1))/step_tau;
       b=cdf_tau(min_pos-1)-a*tau_vector_points(min_pos-1);
       tau_est(qq)=(u_tau-b)/a;
    end % endif

end % end qq

tau_est_no_interpolation_matrix(rr,:)=(tau_est_no_interpolation);


%%----------- end of generate local delays


theta_est_circular_rad=theta_est_circular*pi/180;
% theta_est_circular_rad=[20 30]*pi/180;
    [L_c(rr),g(rr)]=weighting_FactorOFDM(X_hat,S_frequency,theta_est_circular_rad,tau_est_no_interpolation_matrix(rr,:),noise_variance,Q,P,rho_g,rho_Lc,wc,delta_f,rol,c,N);
end
wFactor=exp(L_c-g-max(L_c-g));
tau_est_circular_corrected=max(tau_vector_points)/(2*pi)*angle(real(wFactor)'*exp(1j*2*pi*tau_est_no_interpolation_matrix/max(tau_vector_points)));
tau_est_circular_corrected(find(tau_est_circular_corrected<0)) = max(tau_vector_points)+ tau_est_circular_corrected(find(tau_est_circular_corrected<0));

%%---------------end of circular delay using circular doa

% tau_est_no_interpolation_matrix
% [val_test,pos_test]=min(sum((tau_est_no_interpolation_matrix-repmat([8 10],R,1)).^2,2))
% tau_est_no_interpolation_matrix(pos_test,:)
% pos_tau_circular=find(wFactor)
% max_points
% tau_est_no_interpolation_=1/R*sum(tau_est_no_interpolation_matrix)
% tau_est_linear=1/R*sum(tau_est_matrix);

