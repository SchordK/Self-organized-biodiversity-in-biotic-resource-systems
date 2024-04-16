% Implements the chasing pair & intraspecific interference model as described in:
% This script runs simulations while scanning the parameters for species abundance.
% Click the plot to show the (stochastic) trajectories for species abundances vs time.

%eg. Fig.2A
clear
tic
clc
N=2;
 % set the model parameters
Alpha = 1.25;
par.a = 0.05;par.u = 0.0625; 
par.d= 2; par.p=0;
par.w = 0.32;par.k= 0.22;
par.R00 =0.15; par.K0 = 10000; 
par.D =[0.055;0.058];
par.v = 0.2; 


%define time mesh
t0=1e5;
tspan = 0:1:t0;
tmesh=linspace(0,t0,t0);

% initial species abundances
x0= zeros(1,N);
Y0= zeros(1,N);
C0= 1000*ones(1,N);
R0=2000;
y0 = [x0 Y0 C0 R0];%N=2

[t,y] = ode45(@(t,y) odefcn(t,y,N,par),tspan,y0);  % ODEs simulation
tra= SSA(par,tmesh,y0,N);  % SSA  simulation

figure;
plot(t,y(:,1+2*N:3*N));hold on %Ci_ODEs results
plot(t,y(:,1+3*N),'r','linewidth',3); %R_ODEs results
c=tra(1:N,:)+2*tra(1+N:2*N,:)+tra(1+2*N:3*N,:);r=sum(tra(1:N,:))+tra(1+3*N,:);
plot(tmesh,r,'r','linewidth',3) ;hold on  %R_SSA results
plot(tmesh,c) ;hold on  %Ci_SSA results
xlabel('Time')
ylabel('Consumers and Resource')
set(gca,'XScale','log')
set(gca,'YScale','log')
% axis([1,t0,1e0,1*1e3]);






toc
