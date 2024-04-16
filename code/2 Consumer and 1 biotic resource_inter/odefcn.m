function dydt = odefcn(t,y,par)
%Chasing pair & interspecific interference model 
% returns a function that can be used to simulate population dynamics
% the 'type' variable determines of species or pair abundance
dydt = zeros(6,1);
dydt(1) = par.a1*(y(4)-y(1)-y(3))*(y(6)-y(1)-y(2))-(par.d1+par.k1)*y(1);    %x1_chasing pair abundance
dydt(2) = par.a2*(y(5)-y(2)-y(3))*(y(6)-y(1)-y(2))-(par.d2+par.k2)*y(2);    %x2_chasing pair abundance
dydt(3) = par.b*(y(4)-y(1)-y(3))*(y(5)-y(2)-y(3))-par.e*y(3);  %z_interspecific interference pair abundance
dydt(4) = par.w1*par.k1*y(1)-par.D1*y(4);  %C1_species abundance
dydt(5) = par.w2*par.k2*y(2)-par.D2*y(5);  %C2_species abundance
dydt(6) = par.R00*y(6)*(1-y(6)/par.K0)-par.k1*y(1)-par.k2*y(2);  %R_species abundance
end






