% Moose implementation of input design.
%
% The example illustrates the syntax used to set power constraints.
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
theta0     = [-0.5 1];
A0         = [1 theta0(1)];
B0         = theta0(2);
F0         = 1;
C0         = 1;
D0         = 1;
lambda     = 1;           
trueSystem = idpoly(A0,B0,C0,D0,F0,lambda,Ts);
% ----- SYSTEM IDENTIFICATION MODEL AND EXPERIMENT
model      = trueSystem;
Nident     = 500;
alpha      = 0.95;
%% INPUT DESIGN USING MOOSE
% The optimal input design problem is specified. A 20 parameter FIR spectrum is used.
optInputDesign                                                  = oidProblem(model,Nident,'FIR',20);

% Excitation power constraint
% conPow                                                          = 0.1;
% optInputDesign.spectrum.excitation.power.lb                     = conPow; 
% optInputDesign.spectrum.excitation.power.ub                     = conPow; 

% Input spectrum constraint (input constraint <=> excitation constraint if open-loop)
% conPow                                                          = 0.1;
% optInputDesign.spectrum.input.power.lb                          = conPow; 
% optInputDesign.spectrum.input.power.ub                          = conPow;  

% Output spectrum constraint
conPow                                                          = 10;
optInputDesign.spectrum.output.power.lb                         = conPow; 
% optInputDesign.spectrum.output.power.ub                         = conPow;  

% The input design problem is solved for the minimum input power and the
% corresponding optimal spectral factor is obtained.
[optH,yalmipInfo,infoMatrix,objValue,signalPower,spectrumParam] = solve(optInputDesign, [1 0 0]);

%% MONTE-CARLO STUDY
NMC        = 100;
thetaHatMC = zeros(2,NMC);
h          = waitbar(0,'Running Monte-Carlo study...');
for tk = 1:NMC
   e                = sqrt(lambda)*randn(Nident,1);
   u                = lsim(optH,randn(Nident,1));
   y                = sim(trueSystem,[u e]);
   modelEst         = arx(iddata(y,u,Ts),model,...
                        'MaxIter',100,'Tolerance',1e-10);
   theta            = getpvec(modelEst);
   thetaHatMC(:,tk) = [theta(1); theta(2)]; 
   waitbar(tk/NMC,h);
end
close(h)
%% PLOT RESULTS
figure;
hold on
ellipse(Nident/chi2inv(alpha,2)*(infoMatrix),theta0,'k');
plot(thetaHatMC(1,:),thetaHatMC(2,:),'xk');
xlabel('\theta_1')
ylabel('\theta_2')
legend('Identification set','True parameters', 'Estimates')

P                 = bodeoptions; P.MagUnits = 'abs'; P.FreqScale = 'linear';
wOut              = linspace(1e-3,pi);

% Check of excitation
r                 = lsim(optH,randn(10000,1));
[mag,~, wout]     = bode(minreal((optH)*(optH)'),wOut);
powerExOpt        = signalPower.excitation
powerExCheck1     = var(r)
powerExCheck2     = trapz(wOut,mag)/pi

% Check of input
m                 = size(tf(model));
Su                = eye(m)/(eye(m) - tf(model)*optInputDesign.controller.K);
u                 = lsim(Su*optH,randn(10000,1));
[mag,~, wout]     = bode(minreal((Su*optH)*(optH*Su)'),wOut);
powerInputOpt     = value(signalPower.input.tot)
powerInputCheck1  = var(u)
powerInputCheck2  = trapz(wOut,mag)/pi

% Check of output
H                    = tf(model);
Sy                   = eye(m)/(eye(m) - tf(model)*optInputDesign.controller.K)*tf(model);
y_r                  = lsim(minreal(Sy*optH),randn(10000,1));
[mag,~, wout]        = bode(minreal((Sy*optH)*(Sy*optH)'),wOut);
powerOutputOpt_r     = value(signalPower.output.r)
powerOutput_r_Check1 = var(y_r)
powerOutput_r_Check2 = trapz(wOut,mag)/pi
y_e                  = lsim(minreal(eye(m)/(eye(m) - tf(model)*optInputDesign.controller.K)*H),sqrt(lambda)*randn(10000,1));
powerOutputOpt_e     = value(signalPower.output.e)
powerOutput_e_Check1 = var(y_e)







