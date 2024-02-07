function [u,dotu,dot2u]=OESS2(M,C,K,F,lou,dt)
%2023.10.17
alpha_m=(2-lou)/(1+lou);
alpha_c=(2*lou^2-7*lou+11)/(3*(1+lou)*(3-lou));
alpha_k=1/(1+lou);
beta=1/(1+lou)^2;
gamma=1/2+alpha_m-alpha_k;
t_end=size(F,2);
u=zeros(size(M,1),t_end);dotu=u;
dot2u=M\(F(:,1)-C*dotu(:,1)-K*u(:,1));
M1=alpha_m*M+alpha_c*gamma*dt*C+alpha_k*beta*dt^2*K;
C1=C+K*alpha_k*dt;
K1=(1-alpha_k)*K+alpha_k*K;
for i=2:t_end
      dot2u(:,i)=M1\(alpha_k*F(i)+(1-alpha_k)*F(i-1)-C1*dotu(:,i-1)-K1*u(:,i-1));
      u(:,i)=u(:,i-1)+dt*dotu(:,i-1)+(1/2-beta)*dt^2*dot2u(:,i-1)+beta*dt^2*dot2u(:,i);
      dotu(:,i)=dotu(:,i-1)+(1-gamma)*dt*dot2u(:,i-1)+gamma*dt*dot2u(:,i);
end