% Implements the chasing pair & intraspecific interference model as described in:
% This script runs simulations while scanning the parameters for species abundance.
% Click the plot to show the (stochastic) trajectories for species abundances vs time.

clear
tic
clc

%eg. Fig.S1b,d


% set the model parameters
alpha=1.05;
par.a1 = 0.06; par.a2 = 0.060; 
par.d1 = 0.2; par.d2 = 0.2; 
par.k2 = 0.1; par.k1 = 0.1;
par.K0 = 100; par.w1 = 0.06; par.w2 = 0.06;
par.D1 = 0.0009; par.D2 = 0.00075;
par.R00 = 0.05;par.R0 = 0.05;
par.e = 0.05; par.b = alpha*par.a1;

%define time mesh
t0=1*1e5;
na=0.0*t0;
tspan =0:1:t0;
tmesh=linspace(0,t0,t0);

% initial species abundances
init=[30;10;8;0;0;0];
y0 = [0.0 0.0 0.0 10 8 80.];

%ODEs simulation
[t,y] = ode45(@(t,y) odefcn(t,y,par),tspan,y0);

%SSA simulation
tra = SSA(par,tmesh,init);
mm1=[y(:,4),y(:,5),y(:,6)];

figure;
plot(t,y(:,4),'g','linewidth',1);hold on
plot(t,y(:,5),'b','linewidth',1);
plot(t,y(:,6),'r','linewidth',1);
a=tra(1,:)+tra(4,:)+tra(5,:);b=tra(2,:)+tra(4,:)+tra(6,:);c=tra(3,:)+tra(5,:)+tra(6,:);
plot(tmesh,a,'r','linewidth',2)   %R
plot(tmesh,b,'g','linewidth',2)   %C1
plot(tmesh,c,'b','linewidth',2)   %C2
SSA=[tmesh',b',c',a'];
%plot(t,R,'g','linewidth',1);hold on
legend('C_1','C_2','R','fontname','Times New Roman','FontSize',20','linewidth',1)
xlabel('{Time}','fontname','Times New Roman','FontSize',24');
ylabel('{Species abundances}','fontname','Times New Roman','FontSize',24);
set(gca,'FontName','Times New Roman','FontSize',24,'linewidth',2);
set(gca,'XScale','log');
set(gca,'YScale','log');
%axis([-1000,t0,-0.5,20]);
%title('{Biotic}','fontname','Times New Roman','FontSize',24);
%print(gcf,'-dpng','-r1000','figCase3Biotic,TimeseriesR=0.048.png')
%print(gcf,'-dpng','-r1000','figCase3Biotic,TimeseriesR=0.04.png')







