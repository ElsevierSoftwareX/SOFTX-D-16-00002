% Moose implementation of example based on the example in Section 3.1.4 and Section 5.3
% in Henrik Jansson, "Experiment design with applications in identification for control", 
% PhD thesis. 
% 
% The example requires MOOSE2 and yalmip version 20130201 or above.

% Author: Mariette Annergren & Christian A. Larsson
% Copyright (c) 2015 Mariette Annergren & Christian A. Larsson

clear variables
close all
clc
%% SETUP THE PROBLEM 
% ----- SYSTEM
Ts         = 0.05;  
a          = [-1.99185 2.20265 -1.84083 0.89413];
b          = [0.10276 0.18123];
theta0     = [b a];
A0         = [1 a];
B0         = [0 0 0 b];
F0         = 1;
C0         = 1;
D0         = 1;
lambda     = 0.01;            
trueSystem = idpoly(A0,B0,C0,D0,F0,lambda,Ts);
% ----- SYSTEM IDENTIFICATION MODEL AND EXPERIMENT
model      = trueSystem;
G0         = tf(B0,A0,Ts,'Variable','z^-1');        
Nident     = 100;
alpha      = 0.95;
% ----- ELLIPSOIDAL QUALITY CONSTRAINT
w0_lowBW   = 3;      % low bandwidth of weighting function T
w0_highBW  = 8;      % high bandwidth of weighting function T
xi         = 0.7;
T_lowBW    = tf(w0_lowBW^2,[1 2*xi*w0_lowBW w0_lowBW^2],'Variable','s');     % weighting function
T_highBW   = tf(w0_highBW^2,[1 2*xi*w0_highBW w0_highBW^2],'Variable','s');  % weighting function
V_lowBW    = minreal(c2d(T_lowBW,Ts));
V_highBW   = minreal(c2d(T_highBW,Ts));
Wn         = 1; 
Wd         = 1; 
Xn         = absorbDelay(-G0); 
Yn_lowBW   = V_lowBW;
Yn_highBW  = V_highBW;
Yd         = 1; 
Xd         = 0; 
Kn         = 0; 
Kd         = 0; 
gamma      = 0.1;        
wSamp      = linspace(1e-3,0.99*pi,100);
%% INPUT DESIGN USING MOOSE 
% The optimal input design problem is specified. A 30 parameter FIR spectrum is used.
optInputDesignFreq_lowBW                   = oidProblem(model,Nident,'FIR',30);

% The contraints are added to the input design problem 
optInputDesignFreq_lowBW.constraints{1}    = oidEllipsoidalQualityConstraint(Kd,Kn,Wd,Wn,Xd,Xn,Yd,Yn_lowBW,gamma^2,alpha,wSamp);

% The input design problem is solved for the minimum input power and the
% corresponding optimal spectral factor is obtained.
[optHFreq_lowBW,~,~,~,optPowerFreq_lowBW]  = solve(optInputDesignFreq_lowBW, [1 0 0]);

% The optimal input design problem is specified. A 30 parameter FIR spectrum is used.
optInputDesignFreq_highBW                  = oidProblem(model,Nident,'FIR',30);

% The contraints are added to the input design problem 
optInputDesignFreq_highBW.constraints{1}   = oidEllipsoidalQualityConstraint(Kd,Kn,Wd,Wn,Xd,Xn,Yd,Yn_highBW,gamma^2,alpha,wSamp);

% The input design problem is solved for the minimum input power and the
% corresponding optimal spectral factor is obtained.
[optHFreq_highBW,~,~,~,optPowerCon_highBW] = solve(optInputDesignFreq_highBW, [1 0 0]);

%% PLOT RESULTS
figure;
suptitle('Low band width of weighting function T')
hold on
bodemag(V_lowBW,G0,optHFreq_lowBW*optHFreq_lowBW')
legend('T','G0','\Phi_u')
axis([1e-3 80 -120 60 ])
figure;
suptitle('High bandwidth of weighting function T')
hold on
bodemag(V_highBW,G0,optHFreq_highBW*optHFreq_highBW')
legend('T','G0','\Phi_u')
axis([1e-3 80 -120 60 ])














