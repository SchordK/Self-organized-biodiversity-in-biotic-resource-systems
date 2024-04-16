% Implements the chasing pair & intraspecific interference model as described in:
% This script runs simulations while scanning the parameters for species abundance.
% Click the plot to show the (stochastic) trajectories for species abundances vs time.

clear
tic
clc
%eg. Fig.S5C-D                  
N=200;%the type of consumer species
sn=1;%seed
s=randi([100,1000],1,sn);

% set the model parameters
Alpha = 1.25;
par.a1 = 0.2*ones(1,N);   par.k1= 0.3*ones(1,N);   par.d1= 0.4*ones(1,N);
par.a2 = 0.2*ones(1,N);   par.k2=0.3*ones(1,N);    par.d2= 0.4*ones(1,N);
par.a3 = 0.2*ones(1,N);   par.k3= 0.3*ones(1,N);   par.d3= 0.4*ones(1,N);
par.w1 = 0.3*ones(1,N);   par.w2= 0.3*ones(1,N);   par.w3= 0.3*ones(1,N);
par.u = Alpha*par.a1; par.v = 0.2*ones(1,N);par.p = zeros(1,N);
par.R01 =0.85; par.K01 = 4000; 
par.R02 =0.95; par.K02 = 2000; 
par.R03 =0.9; par.K03 = 3000; 
rng(110)
par.D = rand(1,N)*0.03+0.02  ;

% initial species abundances
x0= zeros(1,3*N);
Y0= zeros(1,N);
C0= 10*ones(1,N);
R0=30*ones(1,3);
y0 = [x0 Y0 C0 R0];

%define time mesh
t0=1e5;
tspan = 0:1:t0;
tmesh=linspace(0,t0,1*t0);


%ODEs simulation
[t,y] = ode45(@(t,y) odefcn(t,y,N,par),tspan,y0);

figure;
hold on
plot(t,y(:,1+4*N:5*N),'linewidth',1);  %Ci
plot(t,y(:,1+5*N),'r','linewidth',3);  %R1
plot(t,y(:,2+5*N),'b','linewidth',3);  %R2
plot(t,y(:,3+5*N),'m','linewidth',3);  %R3
set(gca,'XScale','log')
set(gca,'YScale','log')
%axis([1,t0,1e-2,1*1e5]);



%SSA simulation
ss=1;
tra= SSA(par,tmesh,y0,N,s(ss));  %stochastic simulation eq. 1,2
tt{ss}=tra;

c=tra(1:N,:)+tra(1+N:2*N,:)+tra(1+2*N:3*N,:)+2*tra(1+3*N:4*N,:)+tra(1+4*N:5*N,:);
r1=sum(tra(1:N,:))+tra(1+5*N,:);r2=sum(tra(1+N:2*N,:))+tra(1+5*N,:);r3=sum(tra(1+2*N:3*N,:))+tra(1+5*N,:);
%ssa=[tmesh',c',r1',r2',r3'];
ttc{ss}=tra(1:N,t0)+tra(1+N:2*N,t0)+tra(1+2*N:3*N,t0)+2*tra(1+3*N:4*N,t0)+tra(1+4*N:5*N,t0);
ttr1{ss}=sum(tra(1:N,:))+tra(1+5*N,:);
ttr2{ss}=sum(tra(1+N:2*N,:))+tra(1+5*N,:);
ttr3{ss}=sum(tra(1+2*N:3*N,:))+tra(1+5*N,:);


figure;
plot(tmesh,r1,'r','linewidth',3) ;hold on  %R1
plot(tmesh,r2,'b','linewidth',3) ;hold on  %R2
plot(tmesh,r3,'m','linewidth',3) ;hold on  %R3
plot(tmesh,c) ;hold on  %Ci
xlabel('Time')
ylabel('Species abundances')
set(gca,'XScale','log')
set(gca,'YScale','log')





toc
