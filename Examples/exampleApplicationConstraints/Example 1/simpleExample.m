% Moose implementation of input design for an FIR system with two parameters.
%
% The true system is an FIR system given by
% y(t) = theta(1)*u(t-1) + theta(2)u(t-2) + e(t)
% where e(t) is zero-mean, white, Gaussian noise with variance lambda.
% 
% The true system is unknown and needs to be estimated in an identification 
% experiment. The application cost only emphasizes estimating theta(2) with 
% high accuracy, see Vapp.m. Therefore, the system identification experiment 
% is designed such that only small deviations from theta0(2) are allowed.
%
% The example requires MOOSE2 and yalmip version 20130201 or above.

% Author: Mariette Annergren & Christian A. Larsson
% Copyright (c) 2015 Mariette Annergren & Christian A. Larsson

clear variables
close all
clc
%% SETUP THE PROBLEM 
oneParam   = 0; % Boolean: 1 - only estimate theta(1), 0 - estimate theta 

% ----- SYSTEM
Ts         = 1;                     
theta0     = [10 -9];
A0         = 1;
B0         = theta0;
F0         = 1;
C0         = 1;
D0         = 1;
lambda     = 1;             
trueSystem = idpoly(A0,B0,C0,D0,F0,lambda,Ts);
if oneParam; trueSystem.Structure.b.Free = [oneParam 0]; end

% ----- APPLICATION REGION
% The application region is defined to allow small deviations in the estimate of
% theta(2) and to not consider the theta(1), see Vapp.m.
if oneParam; VappHessian = 0.01; else VappHessian = hessian(@(x)Vapp(x),theta0); end
gamma       = 100;          % Allowed application degradation
alpha       = 0.95;         % Confidence level of estimated model

% ----- SYSTEM IDENTIFICATION MODEL AND EXPERIMENT
model       = trueSystem;
Nident      = 500;
%% INPUT DESIGN USING MOOSE
% The optimal input design problem is specified. A 20 parameter FIR
% spectrum is used.
optInputDesign                = oidProblem(model,Nident,'FIR',20);

% The application contraint is added to the input design problem
optInputDesign.constraints{1} = oidApplicationConstraint(VappHessian,gamma,alpha);

% The input design problem is solved for the minimum input power and the
% corresponding optimal spectral factor is obtained.
[optH, info, iF, optValue]          = solve(optInputDesign,[1 0 0]);
%% MONTE-CARLO STUDY
NMC = 100;
if oneParam; thetaHatMC = zeros(1,NMC); else thetaHatMC = zeros(2,NMC); end
h   = waitbar(0,'Running Monte-Carlo study...');
for tk = 1:NMC
   e     = sqrt(lambda)*randn(Nident,1);
   u     = lsim(optH,randn(Nident,1));
   y     = lsim(trueSystem,u) + e;
   theta = getpvec(oe(iddata(y,u,Ts),model,'MaxIter',100,'Tolerance',1e-10));
   if oneParam; thetaHatMC(tk) = theta(1); else thetaHatMC(:,tk) = [theta(1); theta(2)]; end
   waitbar(tk/NMC,h);
end
close(h)
%% PLOT RESULTS
if oneParam
    figure;
    hold on
    appDeviation = sqrt(2/(gamma*VappHessian));
    xApp         = linspace(theta0(1)-appDeviation,theta0(1)+appDeviation,100);
    idDeviation  = sqrt(chi2inv(alpha,1)/(Nident*optInputDesign.informationMatrix));
    xId          = linspace(theta0(1)-idDeviation,theta0(1)+idDeviation,100);
    appLine      = plot(xApp,theta0(2)*ones(length(xApp),1),'r','LineWidth',4);
    idLine       = plot(xId,theta0(2)*ones(length(xId),1),'k','LineWidth',2);
    trueParam    = plot(theta0(1),theta0(2),'xb','MarkerSize',10);
    estParam     = plot(thetaHatMC(:),theta0(2)*ones(NMC,1),'xk','MarkerSize',10);    
    legend('Application set','Identification set', 'True parameters', 'Estimates')   
    xlim([theta0(1)-appDeviation-0.1 theta0(1)+appDeviation+0.1])
    ylim([-9.1 -8.9])
    xlabel('\theta_1')
    ylabel('\theta_2')
    xlabel('\theta_1')
    ylabel('\theta_2')
    title('Estimating one parameter')
else
    figure;
    hold on
    ellipse(gamma/2*VappHessian,theta0,'r');
    ellipse(Nident/chi2inv(alpha,2)*(optInputDesign.informationMatrix),theta0,'k','empty');
    plot(thetaHatMC(1,:),thetaHatMC(2,:),'xk');
    xlim([9.9 10.1])
    ylim([-9.1 -8.9])
    xlabel('\theta_1')
    ylabel('\theta_2')
    xlabel('\theta_1')
    ylabel('\theta_2')
    legend('Application set', 'True parameters','Identification set', 'Estimates')
    title('Estimating two parameters')
end
