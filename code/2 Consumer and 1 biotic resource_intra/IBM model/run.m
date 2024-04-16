% Implements the chasing pair & intraspecific interference model as described in:
% This script runs simulations while scanning the parameters for species abundance  and position.
% Click the plot to show the (stochastic) trajectories for species abundances vs time and position.
%eg. Fig.3dclear
clc
tic
%define time mesh
t0=1e5;
tmesh=linspace(0,t0,t0);
T=tmesh(1:end-1);
tspan = [0 t0];

%define initial condition
init=[100;60;60;0;0;0;0];
y0 = [0.0 0.0 0.0 0.0 60 60 100];

% set the model parameters
%define reaction rate
alpha=1.25;
par.a=0.1;   %a1 d1 k1
par.a1=par.a*alpha;par.k1=0;  %a1' d1' k1'
par.aa=0.1;   %a2 d2 k2
par.aa1=par.a*alpha;par.kk1=0;    %a2' d2' k2'
par.R00=0.7;par.K0=1000;
par.DD = 0.02;par.D =0.023;
par.w = 0.25; par.ww = 0.25;
par.k = 0.1; par.kk = 0.1;
par.d1 = 0.5; par.dd1 = 0.5; %d'
par.d=0.3;par.dd=0.3;

width=50;   %x=[-width,width]
data=cell(3,2);


par.r=5;
time=[0 1e1 1e2 1e3 1e4 1e5 1e6];
[position ibm]=simulation_model_movie(init,width,tmesh,par,time);
a=ibm(1,:)+ibm(4,:)+ibm(5,:);b=ibm(2,:)+ibm(4,:)+2*ibm(6,:);c=ibm(3,:)+ibm(5,:)+2*ibm(7,:);
%IBM=[tmesh',b',c',a'];


figure%('visible','off')
hold on
plot(T,a,'r');   %R
plot(T,b,'g')   %C1
plot(T,c,'b')   %C2
legend('Resource','Consumer1','Consumer2')
xlabel('Time')
ylabel('species abundance')
% set(gca,'XScale','log')
% set(gca,'YScale','log')
axis([1,inf,0,inf])





toc
%save databiotic