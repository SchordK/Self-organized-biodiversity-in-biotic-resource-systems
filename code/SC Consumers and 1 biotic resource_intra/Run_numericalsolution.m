% Implements the chasing pair & intraspecific interference model as described in:
% This script runs simulations while scanning the parameters for species abundance.
% Click the plot to show the (stochastic) trajectories for species abundances vs time.

clear
tic
clc

%eg. Fig.4
N=20; %the type of consumer species

% set the model parameters
Alpha = 1.25;
par.a = 0.01;par.u = Alpha*par.a; 
par.d= 0.8;
par.w = 0.16;par.k= 0.18;
par.p = 0.0; 
par.R00 =0.16;
par.K0 = 180000; 
rng(100)
%par.a=rand(1,N)*0.0+0.3;
par.D=normrnd(1.0,0.05,1,N)*0.015;
%par.v = 0.5; % 分岔参数--oscillation
par.v = 0.05; %分岔参数--stable


% initial species abundances
x0= zeros(1,N);
Y0= zeros(1,N);
C0= 50*ones(1,N);
R0=80;
y0 = [x0 Y0 C0 R0];

%define time mesh 
t0=1*1e5;
tspan = 0:1:t0;
tmesh=linspace(0,t0,t0);

%ODEs simulation
[t,y] = ode45(@(t,y) odefcn(t,y,N,par),tspan,y0);
%ode=[t,y(:,1+2*N:3*N),y(:,1+3*N)];


figure
plot(t,y(:,1+2*N:3*N));hold on
plot(t,y(:,1+3*N),'r','linewidth',3);
set(gca,'XScale','log')
set(gca,'YScale','log')
%axis([1,t0,1e-1,1*1e5]);

%SSA simulation
tra= SSA(par,tmesh,y0,N);
c=tra(1:N,:)+2*tra(1+N:2*N,:)+tra(1+2*N:3*N,:);r=sum(tra(1:N,:))+tra(1+3*N,:);
%ssa=[tmesh',c',r'];

figure;
plot(tmesh,r,'r','linewidth',2) ;hold on   %R
plot(tmesh,c);hold on    %C1
xlabel('Time')
ylabel('Consumers and Resource')
set(gca,'XScale','log')
set(gca,'YScale','log')
%axis([1,t0,1e0,1*1e2]);
toc
