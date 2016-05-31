% Moose implementation of input design for a system to be controlled by
% proportional state feedback control.
%
% The true system is a first order output-error system given by
%  x(t+1) = a0*x(t) + b0*u(t)
%  y(t)   = c0*x(t) + e(t)
%  u(t)   = f*x(t)
%  where e(t) is zero-mean, white, Gaussian noise with variance le.
% 
% The optimal controller is given by f0 = a0/b0. The true system is unknown and
% needs to be estimated in an identification experiment. The estimated system
% results in the controller f = a/b. For good control performance, f should be
% close to f0. Therefore, the system identification experiment is designed such
% that only small deviations from f0 are allowed.
%
% Note that the application is proportional state feedback control, but the
% identification is to be performed in open-loop.
%
% The problem appears as Example 2 in:
%  "Application-oriented Input Design in System Identification: Optimal input
%   design for control", to be submitted.
%
% The example requires MOOSE2 and yalmip version 20130201 or above.

% Author: Mariette Annergren & Christian A. Larsson
% Copyright (c) 2015 Mariette Annergren & Christian A. Larsson

clear variables
close all
clc
%% SETUP THE PROBLEM 
% ----- SYSTEM
Ts                  = 1; 
a0                  = -0.95;
b0                  = 1;
c0                  = 1;
le                  = 1;
trueSystem          = ss(a0,b0,c0,0,Ts);

% ----- OPTIMAL CONTROLLER
f0                  = a0/b0;

% ----- APPLICATION REGION
% The application region is defined to allow small deviations in the estimate of
% f0. One such region is an ellipse with infinite semi-axis in the direction 
% [f; 1] and finite semi-axis in the direction [1; -f]
e                   = diag([0,1]);     % Eigenvalues of the matrix defining the application region
v                   = [f0, 1; 1 -f0];  % Eigenvectors of the application region matrix
application         = v*e/v;
gamma               = 100;             % The allowed application degradation
alpha               = 0.95;            % Confidence level for the estimated model

% ----- SYSTEM IDENTIFICATION MODEL AND EXPERIMENT
model               = idpoly(trueSystem);
model.NoiseVariance = le;
Nident              = 500;
%% INPUT DESIGN USING MOOSE
% The optimal input design problem is specified. A 10 parameter Moving Average
% spectrum is used.
% The input design uses a transfer function model instead of state space model.
% The parameters of the two models are related as
%  [theta_tf(1); theta_tf(2)] = [theta_ss(1); -theta_ss(2)]
% The application cost Hessian needs to be transformed accordingly in the design
%  VappHessian_tf = [1, 0; 0 -1]'*VappHessian_ss*[1, 0; 0 -1] 
%                 =: T*VappHessian_ss*T
T                               = [1,0;0,-1];
optInputDesign                  = oidProblem(model,Nident,'MA',10);

% The application contraint is added to the input design problem
optInputDesign.constraints{1}   = oidApplicationConstraint(T'*application*T,gamma,alpha);

% The input design problem is solved for the minimum input power and the
% corresponding optimal spectral factor is obtained.
optH                            = solve(optInputDesign,[1 0 0]);

% For comparison, a white design that satisfies the application specifications
% is also tried. This can be achieved with an MA (FIR) spectrum with one parameter.
whiteInputDesign                = oidProblem(model,Nident,'MA',1);
whiteInputDesign.constraints{1} = oidApplicationConstraint(T'*application*T,gamma,alpha);
whiteH                          = solve(whiteInputDesign,[1 0 0]);

%% COMPARISON OF THE TWO SOLUTIONS
e      = randn(Nident,1);
whiteU = lsim(whiteH,e);
optU   = lsim(optH,e);

% Some normalization for fair comparison
umin   = min(min(whiteU),min(optU)); 
umax   = max(max(whiteU),max(optU));
if abs(umin)<abs(umax)
   optUNorm   = optU/abs(umin);
   whiteUNorm = whiteU/abs(umin);
else
   optUNorm   = optU/abs(umax); 
   whiteUNorm = whiteU/abs(umax); 
end

%% PLOT RESULTS
figure;
subplot(1,2,1)
plot(1:Nident,optUNorm)
xlim([1,Nident]); ylim([-1,1]);
xlabel('Time'); ylabel('Input signal')
title('Optimal input signal')

subplot(1,2,2)
plot(1:Nident,whiteUNorm);
xlim([1,Nident]); ylim([-1,1]);
xlabel('Time'); ylabel('Input signal')
title('White input signal')


