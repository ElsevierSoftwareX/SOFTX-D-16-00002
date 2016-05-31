% Moose implementation of example created by Patricio E. Valenzuela
%
% The example requires MOOSE2 and yalmip version 20130201 or above.

% Author: Patricio E. Valenzuela
% Copyright (c) 2015 Mariette Annergren & Christian A. Larsson

clear variables
close all
clc
%% SETUP THE PROBLEM 
% ---- SYSTEM
theta0     = [2, 0.5];
A0         = 1;
B0         = [0 theta0(1)];
F0         = [1 theta0(2)];
C0         = 1;
D0         = 1;
lambda     = 1;
Ts         = 1;
trueSystem = idpoly(A0,B0,C0,D0,F0,lambda,Ts);
K          = tf(-0.2,1,Ts); % Controller
G          = tf(B0,F0,Ts);
Gc         = minreal(G/(1+G*K));
% ----- SYSTEM IDENTIFICATION MODEL AND EXPERIMENT
model      = trueSystem;
Nident     = 5e2;
% --- SPECTRUM CONSTRAINTS
% - Excitation
Gr_ub = tf([0 0.5], [1 -0.5], Ts);
Gr_lb = tf([0 0.1], [1 -0.5], Ts);
% - Input
Gu_ub = tf([0 10], [1 -0.8], Ts);
Gu_lb = tf([0 0.1], [1 -0.8], Ts);
% - Output
Gy_ub = tf([0 10], [1 -0.8], Ts);
Gy_lb = tf([0 0.1], [1 -0.8], Ts);
% --- QUALITY CONSTRAINTS
T       = minreal(G*K/(1+G*K));
gamma_q = 0.5;
con(1)  = oidQualityConstraint(T,3,[],2);
con(2)  = oidQualityConstraint(T,gamma_q,[],inf);

%% INPUT DESIGN USING MOOSE
% The optimal input design problem is specified. A 10 parameter FIR spectrum is used.
prob = oidProblem(model,Nident,'FIR',10,K,'fixed');

% The contraints are added to the input design problem 
prob.spectrum.excitation.ub = Gr_ub;
prob.spectrum.excitation.lb = Gr_lb;
prob.spectrum.input.ub      = Gu_ub;
prob.spectrum.input.lb      = Gu_lb;
prob.spectrum.output.ub     = Gy_ub;
prob.spectrum.output.lb     = Gy_lb;
prob.constraints{1}         = con(1);
prob.constraints{2}         = con(2);

% The input design problem is solved and the
% corresponding optimal spectral factor is obtained.
[optH,info,iF,optVal,signalPow,c,con] = solve(prob,[1,0.5,2]);

%% PLOT RESULTS
figure;
suptitle('Excitation spectrum constraint')
hold on
bodemag(optH*optH',Gr_ub*Gr_ub',Gr_lb*Gr_lb')
legend('\Phi_r', 'upper bound','lower bound')
axis([0.5*1e-2 pi -60 60 ])

figure;
suptitle('Input spectrum constraint')
hold on
bodemag(G*optH*optH'*G',Gu_ub*Gu_ub',Gu_lb*Gu_lb')
legend('\Phi_u', 'upper bound','lower bound')
axis([0.5*1e-2 pi -100 100 ])

figure;
suptitle('Output spectrum constraint')
hold on
bodemag(Gc*optH*optH'*Gc',Gy_ub*Gy_ub',Gy_lb*Gy_lb')
legend('\Phi_y', 'upper bound','lower bound')
axis([0.5*1e-2 pi -100 100 ])

figure;
suptitle('Quality constraint')
hold on
bodemag(T,G)
legend('T', 'G')
axis([0.5*1e-2 pi -100 100 ])


