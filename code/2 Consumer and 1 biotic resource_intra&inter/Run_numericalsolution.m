% Implements the chasing pair & intraspecific interference model as described in:
% This script runs simulations while scanning the parameters for species abundance.
% Click the plot to show the (stochastic) trajectories for species abundances vs time.

clear
tic
clc

%eg. Fig.S4g-h

% set the model parameters
 Alpha = 1.1; par.a1 = 0.05; par.a2 = 0.05;
par.k1 = 0.11; par.k2 = 0.11;
par.d1 = 0.2; par.d2 = 0.2;
par.u1 =Alpha*par.a1; par.u2 = Alpha*par.a2;
par.p1 = 0.00; par.p2 = 0.00;
par.b = par.u1; 
par.w1 = 0.12; par.w2 = 0.12;
par.D1 = 0.008; par.D2 = 0.0079;
par.v1 = 0.2;par.v2 = 0.2;
par.K0=200000;
par.e = 0.150;
par.R00 = 0.10; %oscillation
%par.R00 = 0.20; %stable

 %define time mesh
t0=1e5;
tspan = 0:0.1:t0;
tmesh=linspace(0,t0,t0);

% initial species abundances
init=[80;30;30;0;0;0;0;0];
y0 = [0.0 0.0 0.0 0.0 0.0 30 30 80];

%ODEs simulation
[t,y] = ode45(@(t,y) odefcn(t,y,par),tspan,y0);

%SSA simulation
tra = SSA(par,tmesh,init);

figure;
plot(t,y(:,8),'r','linewidth',2);hold on
plot(t,y( :,6),'g','linewidth',2);
plot(t,y( :,7),'b','linewidth',2);
a=tra(1,:)+tra(4,:)+tra(5,:);b=tra(2,:)+tra(4,:)+2*tra(6,:)+tra(8,:);c=tra(3,:)+tra(5,:)+2*tra(7,:)+tra(8,:);
plot(tmesh,a,'r','linewidth',1)   %R
plot(tmesh,b,'g','linewidth',1)   %C1
plot(tmesh,c,'b','linewidth',1)   %C2
set(gca,'XScale','log');
set(gca,'YScale','log');


toc


