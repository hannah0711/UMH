function [zhatInit,zhat]= positionesttimation(delay,c,jj,DoaALL)
% % Init step
QDelay=cov((delay(2:end,:).'))*c^2 ;
QDoa=diag([var(DoaALL(1,:)),var(DoaALL(2,:)),var(DoaALL(3,:)),var(DoaALL(4,:))]);
bs1=[146.2 , -172.6]*1e-2;
bs2=[841.6 , -213.2]*1e-2;
bs3=[907.7 , 338.2]*1e-2;
bs4=[19.1 , 333.4]*1e-2;

x1=bs1(1);
x2=bs2(1);
x3=bs3(1);
x4=bs4(1);

y1=bs1(2);
y2=bs2(2);
y3=bs3(2);
y4=bs4(2);

% d21=sqrt((bs1-bs2)*(bs1-bs2).');
% d31=sqrt((bs1-bs3)*(bs1-bs3).');
% d41=sqrt((bs1-bs4)*(bs1-bs4).');

 d21=c*(delay(2,jj)-delay(1,jj)) ;
d31=c*(delay(3,jj)-delay(1,jj))  ;
d41=c*(delay(4,jj)-delay(1,jj)) ;

A=[[x1-x2;x1-x3;x1-x4] [y1-y2;y1-y3;y1-y4]]; 
b1= 0.5*[x1^2-x2^2 + y1^2-y2^2 + d21^2; x1^2-x3^2 + y1^2-y3^2 + d31^2; x1^2-x4^2 + y1^2-y4^2 + d41^2] ;
b2=[d21,d31,d41].' ;
% zhatInit=inv(A.'*inv(QDelay)*A)*A.'*inv(QDelay)*(b1+b2*d1hat);
% inv(QDelay)
% zhatInit1=inv(A.'*A)*A.'*(b1+b2*d1hat) 
aVect=inv(A.'*inv(QDelay)*A)*A.'*inv(QDelay)*(b1) ;
cVect=inv(A.'*inv(QDelay)*A)*A.'*inv(QDelay)*(b2) ;

% aVect=inv(A.'*A)*A.'*(b1) 
% cVect=inv(A.'*A)*A.'*(b2) 

 b = 2*cVect(1,1)*(aVect(1,1)-x1)+2*cVect(2,1)*(aVect(2,1)-y1) ;
 a = cVect(1,1)^2+cVect(2,1)^2-1 ;
 cc = (aVect(1,1) - x1)^2+(aVect(2,1) - y1)^2 ;
 
d1hat1=(-b+sqrt(b^2-4*a*cc))/(2*a) ;
d1hat2=(-b-sqrt(b^2-4*a*cc))/(2*a) ;
if(d1hat1>0)
  d1hat=d1hat1;
else
  d1hat=d1hat2;
end
%  delay(1,1)
% d1hat2=c*delay(1,1)
% d1hat2
% d1hat
zhatInit=aVect+cVect*d1hat;
% 
% % Hybrid Algorithm
d1primehat=sqrt((bs1-zhatInit.')*(bs1-zhatInit.').') ;
 d2hat=sqrt((bs2-zhatInit.')*(bs2-zhatInit.').') ;
d3hat=sqrt((bs3-zhatInit.')*(bs3-zhatInit.').') ;
d4hat=sqrt((bs4-zhatInit.')*(bs4-zhatInit.').') ;


QDoadhat=QDoa.*diag([d1primehat,d2hat,d3hat,d4hat].^2);
d12hat=d2hat-d1primehat;
d13hat=d3hat-d1primehat;
d14hat=d4hat-d1primehat;
b2hat=[d12hat,d13hat,d14hat].';
temppost1=( bs1 -zhatInit.')/d1hat;
temppost2=( bs2 -zhatInit.')/d2hat;
temppost3=( bs3 -zhatInit.')/d3hat;
temppost4=( bs4 -zhatInit.')/d4hat;
thathat(1,1)=45-DoaALL(1,jj);
thathat(2,1)=135-DoaALL(2,jj);
thathat(3,1)=-135-DoaALL(3,jj);
thathat(4,1)=-45-DoaALL(4,jj);

thetahatVEct=[-sind(thathat) cosd(thathat)];

G=[temppost1-temppost2;temppost1-temppost3;temppost1-temppost4;thetahatVEct];
BSMU=sum(([bs1;bs2;bs3;bs4]-repmat(zhatInit.',4,1)).*[-sind(thathat),cosd(thathat)],2);
h=[b2-b2hat;BSMU];
Q=  blkdiag(QDelay,QDoadhat);
zhat=zhatInit+inv(G.'*inv(Q)*G)*G.'*inv(Q)*h;
