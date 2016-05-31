% Moose implementation of input design.
%
% The example illustrates the syntax used to set different criteria 
% in optimal experiment design (SISO case).
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
theta0     = [0.95 1];
A0         = [1 theta0(1)];
B0         = theta0(2);
F0         = 1;
C0         = 1;
D0         = 1;
lambda     = 1;             % Variance of measurement noise
trueSystem = idpoly(A0,B0,C0,D0,F0,lambda,Ts);
% ----- SYSTEM IDENTIFICATION MODEL AND EXPERIMENT
model       = trueSystem;
Nident      = 500;
alpha       = 0.95;
%% INPUT DESIGN USING MOOSE -- A-Optimality
% The optimal input design problem is specified. A 20 parameter FIR spectrum is used.
optInputDesignA                                                 = oidProblem(model,Nident,'FIR',20);

% Input spectrum constraint (input constraint <=> excitation constraint if open-loop)
conPow                                                          = 1; 
optInputDesignA.spectrum.input.power.ub                         = conPow;  

% The input design problem is solved for the specific optimality criterion and the
% corresponding optimal spectral factor is obtained.
[optHA,yalmipInfoA,infoMatrixA,objValueA,signalPowerA,spectrumParamA,conA] = solve(optInputDesignA,'A');
%% INPUT DESIGN USING MOOSE -- E-Optimality
% The optimal input design problem is specified. A 20 parameter FIR spectrum is used.
optInputDesignE                                                 = oidProblem(model,Nident,'FIR',20);

% Input spectrum constraint (input constraint <=> excitation constraint if open-loop)
conPow                                                          = 1; 
optInputDesignE.spectrum.input.power.ub                         = conPow;  

% The input design problem is solved for the specific optimality criterion and the
% corresponding optimal spectral factor is obtained.
[optHE,yalmipInfoE,infoMatrixE,objValueE,signalPowerE,spectrumParamE,conE] = solve(optInputDesignE,'E');
%% INPUT DESIGN USING MOOSE -- D-Optimality
% The optimal input design problem is specified. A 20 parameter FIR spectrum is used.
optInputDesignD                                                 = oidProblem(model,Nident,'FIR',20);

% Input spectrum constraint (input constraint <=> excitation constraint if open-loop)
conPow                                                          = 1; 
optInputDesignD.spectrum.input.power.ub                         = conPow;  

% The input design problem is solved for the specific optimality criterion and the
% corresponding optimal spectral factor is obtained.
[optHD,yalmipInfoD,infoMatrixD,objValueD,signalPowerD,spectrumParamD,conD] = solve(optInputDesignD,'D');

%% PLOT RESULTS
figure;
hold on
ellipse(Nident/chi2inv(alpha,2)*(infoMatrixA),theta0,'b','empty');
ellipse(Nident/chi2inv(alpha,2)*(infoMatrixE),theta0,'r','empty');
ellipse(Nident/chi2inv(alpha,2)*(infoMatrixD),theta0,'g');
xlabel('\theta_1')
ylabel('\theta_2')
legend('Identification set - A','Identification set - E','Identification set - D','True parameters')

P                = bodeoptions; P.MagUnits = 'abs'; P.FreqScale = 'linear';
wOut             = linspace(1e-3,pi);

% Check of infoMatrix
A_objectiveValueFromYalmip_explicit = objValueA
A_objectiveValueFromYalmip_implicit = trace(inv(infoMatrixA))
E_objectiveValueFromYalmip_explicit = objValueE
D_objectiveValueFromYalmip_explicit = objValueD

% Check of input
m                  = size(tf(model));

uA                 = lsim(optHA,randn(10000,1));
[mag,~, wout]      = bode(minreal(optHA*optHA'),wOut);
A_powerInputOpt    = value(signalPowerA.input.tot)
A_powerInputCheck1 = var(uA)
A_powerInputCheck2 = trapz(wOut,mag)/pi

uE                 = lsim(optHE,randn(10000,1));
[mag,~, wout]      = bode(minreal(optHE*optHE'),wOut);
E_powerInputOpt    = value(signalPowerE.input.tot)
E_powerInputCheck1 = var(uE)
E_powerInputCheck2 = trapz(wOut,mag)/pi

uD                 = lsim(optHD,randn(10000,1));
[mag,~, wout]      = bode(minreal(optHD*optHD'),wOut);
D_powerInputOpt    = value(signalPowerD.input.tot)
D_powerInputCheck1 = var(uD)
D_powerInputCheck2 = trapz(wOut,mag)/pi







