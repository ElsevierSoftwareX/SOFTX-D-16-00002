% Moose implementation of example based on the example in Section 3.1.4 in Henrik Jansson, 
% "Experiment design with applications in identification for control", PhD thesis. 
% (Note that no ellipsoidal quality constraint is imposed! The variance of Delta is bounded instead.)
% 
% The example requires MOOSE2 and yalmip version 20130201 or above.

% Author: Mariette Annergren & Christian A. Larsson
% Copyright (c) 2015 Mariette Annergren & Christian A. Larsson

clear variables
close all
clc
%% SETUP THE PROBLEM 
% ----- SYSTEM
Ts         = 1; 
a          = [-1.99185 2.20265 -1.84083 0.89413];
b          = [0.10276 0.18123];
theta0     = [b a];
A0         = [1 a];
B0         = [0 0 0 b];
F0         = 1;
C0         = 1;
D0         = 1;
lambda     = 0.05;             
trueSystem = idpoly(A0,B0,C0,D0,F0,lambda,Ts);
% ----- SYSTEM IDENTIFICATION MODEL AND EXPERIMENT
model      = trueSystem;
G0         = tf(B0,A0,Ts,'Variable','z^-1');       
Nident     = 500;
% ----- QUALITY CONSTRAINT
w0         = 0.4;
xi         = 0.7;
T          = tf(w0^2,[1 2*xi*w0 w0^2],'Variable','s'); % weighting function 
V          = minreal(c2d(T,Ts));
gamma      = 0.1;         
wSamp      = linspace(1e-3,pi,100);
% ----- INPUT SPECTRUM CONSTRAINT
freqResp   = db2mag([55*ones(1,50) 5*ones(1,50) 40*ones(1,200) -21*ones(1,50)]);
freq       = [linspace(1e-3,0.393-0.04,50) linspace(0.393-0.04,0.393+0.04,50) linspace(0.393+0.04,2,200) linspace(2,pi,50)];
conSys     = struct('B',freqResp,'w',freq);

%% INPUT DESIGN USING MOOSE (with frequency grid and spectrum constraint)
% The optimal input design problem is specified. A 40 parameter FIR spectrum is used.
optInputDesignFreq                        = oidProblem(model,Nident,'FIR',40);

% The contraints are added to the input design problem 
optInputDesignFreq.constraints{1}         = oidQualityConstraint(V,gamma^2,wSamp,inf);
optInputDesignFreq.spectrum.excitation.ub = conSys;

% The input design problem is solved for the minimum input power and the
% corresponding optimal spectral factor is obtained.
[optHFreq,~,~,optPowerFreq]               = solve(optInputDesignFreq,[1 0 0]);

%% INPUT DESIGN USING MOOSE (with frequency grid and no spectrum constraint)
% The optimal input design problem is specified. A 20 parameter FIR spectrum is used.
optInputDesignFreqNoCon                = oidProblem(model,Nident,'FIR',20);

% The contraints are added to the input design problem 
optInputDesignFreqNoCon.constraints{1} = oidQualityConstraint(V,gamma^2,wSamp,inf);

% The input design problem is solved for the minimum input power and the
% corresponding optimal spectral factor is obtained.
[optHFreqNoCon,~,~,optPowerFreqNoCon]  = solve(optInputDesignFreqNoCon, [1 0 0]);

%% PLOT RESULTS
figure;
suptitle('With spectrum constraint')
hold on
bodemag(V,G0,optHFreq*optHFreq')
plot((optInputDesignFreq.spectrum.excitation.ub.w),mag2db(squeeze(optInputDesignFreq.spectrum.excitation.ub.B)),'k')
legend('V','G0','\Phi_u', 'upper bound')
axis([0.5*1e-2 pi -60 60 ])
figure;
suptitle('Without spectrum constraint')
hold on
bodemag(V,G0,optHFreqNoCon*optHFreqNoCon')
legend('V','G0','\Phi_u')
axis([1e-3 pi -60 30 ])















