

function [FIM]=CRLBOFDM(NoiseVariance,delay,doa,wc,delta_f,rol,c,P,N,L)



X = exp(-1j*kron(wc*sin(doa),[0:P-1])*rol/c).';
Z = exp(-1j*kron(delay,(2*pi*([1:N]-floor(N/2)-1)*delta_f))).';
Ad = eye(L);%complex gain matrix
% Ad=pinv(X)*H*pinv(Z.');
A = kron(eye(2),Ad);
S =   kr(Z,X);
X1 = ((-1j*kron(wc*cos(doa),[0:P-1])*rol/c).*exp(-1j*kron(wc*sin(doa),[0:P-1])*rol/c)).';
Z1 = ((-1j*kron(ones(L,1),(2*pi*([1:N]-floor(N/2)-1)*delta_f))).*exp(-1j*kron(delay,(2*pi*([1:N]-floor(N/2)-1)*delta_f)))).';
derivativeStheta =  kr(Z,X1);
derivativeStau =  kr(Z1,X);
D = [derivativeStheta,derivativeStau];
II=eye(size(S*pinv(S)));
FIM = (2/NoiseVariance)*(real(A'*D'*(II-S*pinv(S))*D*A));




 