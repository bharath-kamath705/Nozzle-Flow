%---------AE5327 Project----------------
%Numerical solution of quasi 1D C-D nozzle
%Method used: Macormac's method
%note that ALL variables are non dimensional
%-------------------------------------------

clc
clear all
close all

gamma=1.4;


%CD nozzle Area profile
i_tot=101; %total spatial grid points
L=1;
x=linspace(0,L,i_tot);

A=zeros(1,i_tot);
for i=1:i_tot
    if x(i)<=0.5
        A(i)=5.95-9.9*x(i);
    else
        A(i)=0.313+1.374*x(i);
    end
end

%solution parameters
dx=x(2);
n_tot=5000;
cfl=0.8;


%initial conditions
rho=1-0.5*x;
T=1-0.3*x;
V=(0.1+1.2*x).*sqrt(T);

%initial mach and pressure (resulting from the above conditions)
M_i=V(1,:)./sqrt(T(1,:));
p_i=rho(1,:).*T(1,:);


%predictor step
for n=1:n_tot-1
    a=sqrt(T(n,:));
    del_t=cfl*dx*(a+V(n,:)).^-1;
    
    %note:expressions broken into parts to avoid paranthetical errors
   
    %predictor step (note functions 'fd' and 'bd' defined at the end)
        drho_dtp=-rho(n,:).*fd(V(n,:));
        drho_dtp=drho_dtp -rho(n,:).*V(n,:).*fd(log(A));
        drho_dtp=drho_dtp-V(n,:).*fd(rho(n,:));
        drho_dtp=inv(dx)*drho_dtp;
        
        dV_dtp=-V(n,:).*fd(V(n,:));
        dV_dtp=dV_dtp-inv(gamma).*fd(T(n,:));
        dV_dtp=dV_dtp -inv(gamma)*T(n,:).*rho(n,:).^-1 .*fd(rho(n,:));
        dV_dtp=dV_dtp*inv(dx);
        
        dT_dtp=fd(V(n,:))+V(n,:).*fd(log(A));
        dT_dtp=dT_dtp*-(gamma-1).*T(n,:);
        dT_dtp=dT_dtp-V(n,:).*fd(T(n,:));
        dT_dtp=dT_dtp*inv(dx);
        
        rho(n+1,:)=rho(n,:)+drho_dtp*min(del_t);
        V(n+1,:)=V(n,:)+dV_dtp*min(del_t);
        T(n+1,:)=T(n,:)+dT_dtp*min(del_t);
       
        %corrector step
        drho_dtc=-rho(n+1,:).*bd(V(n+1,:));
        drho_dtc=drho_dtc -rho(n+1,:).*V(n+1,:).*bd(log(A));
        drho_dtc=drho_dtc-V(n+1,:).*bd(rho(n+1,:));
        drho_dtc=inv(dx)*drho_dtc;
        
        dV_dtc=-V(n+1,:).*bd(V(n+1,:));
        dV_dtc=dV_dtc-inv(gamma).*bd(T(n+1,:));
        dV_dtc=dV_dtc -inv(gamma)*T(n+1,:).*rho(n+1,:).^-1 .*bd(rho(n+1,:));
        dV_dtc=dV_dtc*inv(dx);
        
        dT_dtc=bd(V(n+1,:))+V(n+1,:).*bd(log(A));
        dT_dtc=dT_dtc*-(gamma-1).*T(n+1,:);
        dT_dtc=dT_dtc-V(n+1,:).*bd(T(n+1,:));
        dT_dtc=dT_dtc*inv(dx);
        
        %average 
        rho(n+1,:)=rho(n,:)+0.5*(drho_dtp+drho_dtc)*min(del_t);
        V(n+1,:)=V(n,:)+0.5*(dV_dtp+dV_dtc)*min(del_t);
        T(n+1,:)=T(n,:)+0.5*(dT_dtp+dT_dtc)*min(del_t);
        
        %applying boundary conditions
        %inlet
        V(n+1,1)=2*V(n+1,2)-V(n+1,3);
        rho(n+1,1)=1;
        T(n+1,1)=1;
        %outlet
        V(n+1,i_tot)=2*V(n+1,i_tot-1)-V(n+1,i_tot-2);
        rho(n+1,i_tot)=2*rho(n+1,i_tot-1)-rho(n+1,i_tot-2);
        T(n+1,i_tot)=2*T(n+1,i_tot-1)-T(n+1,i_tot-2);
    
end

%final results
M=V(end,:)./sqrt(T(end,:));
p=rho(end,:).*T(end,:);

%------------Analytic solution----------------------------
for i=1:i_tot
    
    %Find analytic mach number using area M relations
    fun = @(M_a,A_a) (A_a*M_a)^2-( (2/(gamma+1))*(1+0.5*M_a^2*(gamma-1)))^((gamma+1)/(gamma-1));
    A_a = A(i);
    AM_fun = @(M_a) fun(M_a,A_a);
    
    if i<i_tot/2
        M_a(i) = fzero(AM_fun,0.2);
        
    else
        M_a(i) = fzero(AM_fun,10);
    end
    
    if isnan(M_a(i))
        M_a(i)=1;
    end
    
    %use isentropic relations
    %find analytic pressure based on analytic mach number
    p_a(i)=(1+0.5*(gamma-1)*M_a(i)^2)^(-gamma/(gamma-1));
end
%-----------------------------------------------------------------

set(groot,'defaultLineLineWidth',2)

figure(1)
plot(x,0.5*A,x,-0.5*A),hold on
%plot(x,zeros(1,i_tot),'-o','MarkerSize',7)
ylabel('A/A_t','FontSize',15),xlabel('x/L','FontSize',15)

figure(2)
hold on,plot(x,M_i,x,M,'k'),plot(x,M_a,'--','MarkerSize',19)
ylabel('M','FontSize',15),xlabel('x/L','FontSize',15)
legend({'Initial condition',strcat(num2str(n_tot),' iterations'),'Exact solution'}... 
    ,'FontSize',16,'Location','North')

figure(3)
plot(x,p_i,x,p,'k'),hold on,plot(x,p_a,'--','MarkerSize',19)
ylabel('p/p_0','FontSize',15),xlabel('x/L','FontSize',15)
legend({'Initial condition',strcat(num2str(n_tot),' iterations'),'Exact solution'}... 
    ,'FontSize',16)

%Function that returns forward difference
function z = fd(x)
y=[x(2:end) 0];
z=y-x;
z(end)=0;
end

%Function that returns backward difference
function z = bd(x)
y=[0 x(1:end-1)];
z=x-y;
z(1)=0;
end
