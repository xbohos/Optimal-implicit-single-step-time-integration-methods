function [u1,dotu1,dot2u1]=OESS1(M,C,K,F,lou,dt,u,dotu,dot2u)
%2023.10.17
alpha_m=(2-lou)/(1+lou);
alpha_c=(2*lou^2-7*lou+11)/(3*(1+lou)*(3-lou));
alpha_k=1/(1+lou);
beta=1/(1+lou)^2;
gamma=1/2+alpha_m-alpha_k;
M1=alpha_m.*M+alpha_c*gamma*dt.*C+alpha_k*beta*dt^2.*K;
M2=(1-alpha_m).*M+(alpha_k-alpha_c)*dt.*C+alpha_c*(1-gamma)*dt.*C+alpha_k*(1/2-beta)*dt^2.*K;
C1=C+K*alpha_k*dt;
K1=K;
dot2u1=M1\(alpha_k*F(:,2)+(1-alpha_k)*F(:,1)-C1*dotu-K1*u-M2*dot2u);
u1=u+dt*dotu+(1/2-beta)*dt^2*dot2u+beta*dt^2*dot2u1;
dotu1=dotu+(1-gamma)*dt*dot2u+gamma*dt*dot2u1;