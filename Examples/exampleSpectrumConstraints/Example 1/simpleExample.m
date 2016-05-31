% Moose implementation of input design.
%
% The example illustrates the syntax used to set signal spectrum constraints (SISO case).
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
lambda     = 1;             % Variance of measurement noise
trueSystem = idpoly(A0,B0,C0,D0,F0,lambda,Ts);
% ----- SYSTEM IDENTIFICATION MODEL AND EXPERIMENT
model       = trueSystem;
Nident      = 500;
alpha       = 0.95;
%% INPUT DESIGN USING MOOSE
% The optimal input design problem is specified. A 20 parameter FIR spectrum is used.
optInputDesign                                                  = oidProblem(model,Nident,'FIR',20);

% Excitation spectrum constraint
%-- constraint is LTI
conSys                                                          = tf(0.1,[1 0.1],1,'Variable','z^-1');
% Note that the LTI lower / upper bound constraint is treated as Phi_r >= / <=  conSys*conSys'.
%-- constraint is STRUCT
[magconSys,~,wSys]                                              = bode(conSys,1e-3:0.01:pi);
conSys                                                          = struct('B',magconSys,'w',wSys);
%-- constraint is DOUBLE
conSys                                                          = 0.1;
optInputDesign.spectrum.excitation.lb                           = conSys; 
% optInputDesign.spectrum.excitation.ub                           = conSys; 

% Input spectrum constraint (input constraint <=> excitation constraint if open-loop)
%-- constraint is LTI
% conSys                                                          = tf(0.1,[1 0.1],1,'Variable','z^-1');
% Note that the LTI lower / upper bound constraint is treated as Phi_u >= / <=  conSys*conSys'.
%-- constraint is STRUCT
% [magconSys,~,wSys]                                              = bode(conSys,1e-3:0.001:pi);
% conSys                                                          = struct('B',magconSys,'w',wSys);
%-- constraint is DOUBLE
% conSys                                                          = 0.1;
% optInputDesign.spectrum.input.lb                                = conSys; 
% optInputDesign.spectrum.input.ub                                = conSys; 

% Output spectrum constraint
%-- constraint is LTI
% conSys                                                          = tf(0.1,[1 0.1],1,'Variable','z^-1');
% Note that the LTI lower / upper bound constraint is treated as Phi_y >= / <=  conSys*conSys'.
%-- constraint is STRUCT
% [magconSys,~,wSys]                                              = bode(conSys,1e-3:0.001:pi);
% conSys                                                          = struct('B',magconSys,'w',wSys);
%-- constraint is DOUBLE
% conSys                                                          = 0.1;
% optInputDesign.spectrum.output.lb                               = conSys; 
% optInputDesign.spectrum.output.ub                               = conSys; 

% The input design problem is solved for the minimum input power and the
% corresponding optimal spectral factor is obtained.
[optH,yalmipInfo,infoMatrix,objValue,signalPower,spectrumParam] = solve(optInputDesign,[1 0 0]);
%% MONTE-CARLO STUDY
NMC        = 100;
thetaHatMC = zeros(2,NMC);
h          = waitbar(0,'Running Monte-Carlo study...');
for tk = 1:NMC
   e                = sqrt(lambda)*randn(Nident,1);
   u                = lsim(optH,randn(Nident,1));
   y                = sim(trueSystem,[u e]);
   modelEst         = arx(iddata(y,u,Ts),model,'MaxIter',100,'Tolerance',1e-10);
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
legend('Application set','Identification set','True parameters', 'Estimates')

P                = bodeoptions; P.MagUnits = 'abs'; P.FreqScale = 'linear';
wOut             = linspace(1e-3,pi);

% Check of excitation
figure;
hold on
if ~isempty(optInputDesign.spectrum.excitation.lb)
    if isa(optInputDesign.spectrum.excitation.lb,'lti')
        bodemag((optInputDesign.spectrum.excitation.lb*optInputDesign.spectrum.excitation.lb'),wOut,'r-^',P)
    elseif isa(optInputDesign.spectrum.excitation.lb,'struct')
        plot(squeeze(optInputDesign.spectrum.excitation.lb.w),squeeze(optInputDesign.spectrum.excitation.lb.B),'r-^')
    else
        plot(wOut,optInputDesign.spectrum.excitation.lb*ones(length(wOut),1),'r-^')
    end
end
if ~isempty(optInputDesign.spectrum.excitation.ub)
    if isa(optInputDesign.spectrum.excitation.ub,'lti')
        bodemag((optInputDesign.spectrum.excitation.ub*optInputDesign.spectrum.excitation.ub'),wOut,'r-V',P)
    elseif isa(optInputDesign.spectrum.excitation.ub,'struct')
        plot(squeeze(optInputDesign.spectrum.excitation.ub.w),squeeze(optInputDesign.spectrum.excitation.ub.B),'r-V')
    else
        plot(wOut,optInputDesign.spectrum.excitation.ub*ones(length(wOut),1),'r-V')
    end
    magTrue = freqresp(minreal(optH*optH'),wOut);
    plot(wOut,squeeze(magTrue),'-*b')
    legend('lower bound', 'upper bound','optimal \Phi_r')
else
    magTrue = freqresp(minreal(optH*optH'),wOut);
    plot(wOut,squeeze(magTrue),'-*b')
    legend('lower bound','optimal \Phi_r')
end
title('\Phi_r')
r                = lsim(optH,randn(10000,1));
[mag,~, wout]    = bode(minreal((optH)*(optH)'),wOut);
powerExOpt       = signalPower.excitation
powerExCheck1    = var(r)
powerExCheck2    = trapz(wOut,mag)/pi

% Check of input
figure;
hold on
G                = minreal(tf(model));
m                = size(G);
Su               = minreal(eye(m)/(eye(m) - optInputDesign.controller.K*G));
if ~isempty(optInputDesign.spectrum.input.lb)
    if isa(optInputDesign.spectrum.input.lb,'lti')
        bodemag((optInputDesign.spectrum.input.lb*optInputDesign.spectrum.input.lb'),wOut,'r-^',P)
    elseif isa(optInputDesign.spectrum.input.lb,'struct')
        plot(squeeze(optInputDesign.spectrum.input.lb.w),squeeze(optInputDesign.spectrum.input.lb.B),'r-^')
    else
        plot(wOut,optInputDesign.spectrum.input.lb*ones(length(wOut),1),'r-^')
    end
end
if ~isempty(optInputDesign.spectrum.input.ub)
    if isa(optInputDesign.spectrum.input.ub,'lti')
        bodemag((optInputDesign.spectrum.input.ub*optInputDesign.spectrum.input.ub'),wOut,'r-V',P)
    elseif isa(optInputDesign.spectrum.input.ub,'struct')
        plot(squeeze(optInputDesign.spectrum.input.ub.w),squeeze(optInputDesign.spectrum.input.ub.B),'r-V')
    else
        plot(wOut,optInputDesign.spectrum.input.ub*ones(length(wOut),1),'r-V')
    end
end
bodemag(minreal(Su*(optH*optH')*Su'),wOut,'-*b',P);
if isempty(optInputDesign.spectrum.input.lb) && isempty(optInputDesign.spectrum.input.ub)
    legend('Optimal \Phi_u')
elseif ~isempty(optInputDesign.spectrum.input.lb) && isempty(optInputDesign.spectrum.input.ub)
    legend('lower bound','Optimal \Phi_u')
elseif isempty(optInputDesign.spectrum.input.lb) && ~isempty(optInputDesign.spectrum.input.ub)
    legend('upper bound','Optimal \Phi_u')
elseif ~isempty(optInputDesign.spectrum.input.lb) && ~isempty(optInputDesign.spectrum.input.ub)
    legend('lower bound','upper bound','Optimal \Phi_u')
end
title('\Phi_u')
u                   = lsim(Su*optH,randn(10000,1));
[mag,~, wout]       = bode(minreal((Su*optH)*(optH*Su)'),wOut);
powerInputOpt_r     = value(signalPower.input.r)
powerInput_r_Check1 = var(u)
powerInput_r_Check2 = trapz(wOut,mag)/pi

% Check of output
figure;
hold on
Sy = minreal(eye(m)/(eye(m) - optInputDesign.controller.K*G)*G);
if ~isempty(optInputDesign.spectrum.output.lb)
    if isa(optInputDesign.spectrum.output.lb,'lti')
        bodemag((optInputDesign.spectrum.output.lb*optInputDesign.spectrum.output.lb'),wOut,'r-^',P)
    elseif isa(optInputDesign.spectrum.output.lb,'struct')
        plot(squeeze(optInputDesign.spectrum.output.lb.w),squeeze(optInputDesign.spectrum.output.lb.B),'r-^')
    else
        plot(wOut,optInputDesign.spectrum.output.lb*ones(length(wOut),1),'r-^')
    end
end
if ~isempty(optInputDesign.spectrum.output.ub)
    if isa(optInputDesign.spectrum.output.ub,'lti')
        bodemag((optInputDesign.spectrum.output.ub*optInputDesign.spectrum.output.ub'),wOut,'r-V',P)
    elseif isa(optInputDesign.spectrum.output.ub,'struct')
        plot(squeeze(optInputDesign.spectrum.output.ub.w),squeeze(optInputDesign.spectrum.output.ub.B),'r-V')
    else
        plot(wOut,optInputDesign.spectrum.output.ub*ones(length(wOut),1),'r-V')
    end
end
bodemag(minreal(Sy*(optH*optH')*Sy'),wOut,'-*b',P);
if isempty(optInputDesign.spectrum.output.lb) && isempty(optInputDesign.spectrum.output.ub)
    legend('Optimal \Phi_y')
elseif ~isempty(optInputDesign.spectrum.output.lb) && isempty(optInputDesign.spectrum.output.ub)
    legend('lower bound','Optimal \Phi_y')
elseif isempty(optInputDesign.spectrum.output.lb) && ~isempty(optInputDesign.spectrum.output.ub)
    legend('upper bound','Optimal \Phi_y')
elseif ~isempty(optInputDesign.spectrum.output.lb) && ~isempty(optInputDesign.spectrum.output.ub)
    legend('lower bound','upper bound','Optimal \Phi_y')
end
title('\Phi_y')
y_r                  = lsim(minreal(Sy*optH),randn(10000,1));
[mag,~, wout]        = bode(minreal((Sy*optH)*(Sy*optH)'),wOut);
powerOutputOpt_r     = value(signalPower.output.r)
powerOutput_r_Check1 = var(y_r)
powerOutput_r_Check2 = trapz(wOut,mag)/pi








