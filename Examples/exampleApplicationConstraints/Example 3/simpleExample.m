% Illustration of concepts in the application-oriented input design framework for
% model based control design. 
%
% The true system is a first order output-error system given by
%  x(t+1) = a0*x(t) + b0*u(t)
%  y(t)   = c0*x(t) + e(t)
%  where e(t) is zero-mean, white, Gaussian noise with variance le.
%
% An MPC controller is designed for reference tracking with cost function
%  J(t) = (Y(t) - R(t))'*(Y(t) - R(t))
% where Y(t) contains predicted outputs over the MPC prediction horizon. 
%
% For the design of the MPC, the system needs to be estimated. An application
% cost is chosen to reflect the reference tracking capabilities of the MPC. The
% application set is characterized using the ellipsoidal approximation and the
% scenario approach.
%
% The impact of the MPC design on the application set is illustrated by 
% considering two cases:
%  1. MPC without constraints.
%  2. MPC with constraints on the input, which reduces the controller bandwidth.
%
% The entire Input Design Framework (IDF) is then illustrated in the steps:
%  1. Short initial identification of the true system
%  2. Input design based on the model estimated in step 1.
%  3. Longer identification using the optimally designed inputs.
% A Monte Carlo study is performed for the IDF on the two MPC cases.
% 
% The problem is presented in:
%  "Application-oriented Input Design in System Identification: Optimal input
%   design for control", to be submitted.
% 
% NB! The code takes long time to execute due to extensive numerical
% calculations.
%
% The example requires MOOSE2 and yalmip version 20130201 or above.

% Author: Mariette Annergren & Christian A. Larsson
% Copyright (c) 2015 Mariette Annergren & Christian A. Larsson

clear variables
close all
clc
%% SETUP THE PROBLEM
% ----- SYSTEM
Ts                         = 1;
theta0                     = [0.6; 0.9];
le                         = 1;
trueSystem                 = ss(theta0(2),1,theta0(1),0,Ts);

% ----- CONTROLLER PARAMETERS
% The MPC is implemented in the function simpleMPC.m
lambda                     = 0;
Ny                         = 5;
Nu                         = 5;
Ub                         = [1, -1];
Yb                         = [1, -1];

% ----- APPLICATION REGION
gamma                      = 200;
Napp                       = 10;
Nscen                      = 5000;
Nnum                       = 51;
alpha                      = 0.99;

% -------- Reference trajectory
r                          = zeros(Napp,1);
r(1:3)                     = 1;
r(4:7)                     = -1;

% -------- Nominal responses
[~,yUC]                    = Vapp(theta0,theta0,0*r,r,Ny,Nu,[inf,-inf],[inf,-inf]);
[~,yC]                     = Vapp(theta0,theta0,0*r,r,Ny,Nu,Ub,Yb);

% -------- Hessian of application cost for the two cases
VappHessUC                 = hessian(@(x)Vapp(x,theta0,yUC,r,Ny,Nu,[inf, -inf],[inf, -inf]),theta0);
VappHessC                  = hessian(@(x)Vapp(x,theta0,yC,r,Ny,Nu,Ub,Yb),theta0);

% -------- Scenario calculations for the two cases
scenariosUC                = 0.2-0.4*rand(2,Nscen) + repmat(theta0,1,Nscen);
scenariosC                 = 0.4-0.8*rand(2,Nscen) + repmat(theta0,1,Nscen);
VappScenariosUC            = zeros(1,Nscen);
VappScenariosC             = zeros(1,Nscen);
h                          = waitbar(0,'Running scenario evaluation of Vapp (based on true system)...');
for ik = 1:Nscen
    VappScenariosUC(end,ik) =...
      Vapp(scenariosUC(:,ik),theta0,yUC,r,Ny,Nu,[inf,-inf],[inf,-inf]);
    VappScenariosC(end,ik)  =...
      Vapp(scenariosC(:,ik),theta0,yC,r,Ny,Nu,Ub,Yb);  
    waitbar(ik/Nscen,h);
end
close(h)
VappScenariosMatrixUC      = [scenariosUC; VappScenariosUC];
VappScenariosMatrixC       = [scenariosC; VappScenariosC];
VappScenariosMatrixUC(2,:) = -VappScenariosMatrixUC(2,:);
VappScenariosMatrixC(2,:)  = -VappScenariosMatrixC(2,:);

% ------- Numerical evaluation of Vapp for contours for the two cases
VappNumUC                  = zeros(Nnum);
VappNumC                   = zeros(Nnum);
delta                      = linspace(-0.4,0.4,Nnum);
% h                          = waitbar(0,'Running contour evaluation of Vapp (based on true system)...');
% for ik = 1:Nnum
%    for jk = 1:Nnum
%       VappNumUC(ik,jk) = ...
%          Vapp(theta0+[delta(ik);delta(jk)],theta0,yUC,r,Ny,Nu,[inf, -inf],[inf,-inf]);
%       VappNumC(ik,jk)  = ...
%          Vapp(theta0+[delta(ik);delta(jk)],theta0,yC,r,Ny,Nu,Ub,Yb);
%    end
%    waitbar(ik/Nnum,h);
% end
% close(h)
load VappNumUC
load VappNumC
%% IDF PROCEDURE USING MOOSE
% ----- SYSTEM IDENTIFICATION EXPERIMENT
Ninitial                                 = 100;
Nident                                   = 400;

% ----- INITIAL IDENTIFICATION USING WHITE INPUT SIGNAL
u                                        = randn(Ninitial,1);
y                                        = lsim(trueSystem,u) + sqrt(le)*randn(Ninitial,1);
Z                                        = iddata(y,u,Ts);
modelEstInitial                          = oe(Z,[1,1,1],'MaxIter',30,'Tolerance',1e-4);
thetaHat                                 = getpvec(modelEstInitial); thetaHat(2) = -thetaHat(2);

% ----- CONSTRUCTION OF APPLICATION REGION (BOTH ELLIPSOIDAL APPROXIMATION AND SCENARIO APPROACH)
% -------- Initial model responses
[~,yInitialUC]                           = Vapp(thetaHat,thetaHat,0*r,r,Ny,Nu,[inf,-inf],[inf,-inf]);
[~,yInitialC]                            = Vapp(thetaHat,thetaHat,0*r,r,Ny,Nu,Ub,Yb);

% -------- Hessian of application cost for the two cases based on the intial model
VappHessInitialUC                        = hessian(@(x)...
                                    Vapp(x,thetaHat,yInitialUC,r,Nu,Ny,[inf,-inf],[inf,-inf]),thetaHat);
VappHessInitialC                         = hessian(@(x)...
                                    Vapp(x,thetaHat,yInitialC,r,Nu,Ny,Ub,Yb),thetaHat);

% -------- Scenario calculations for the two cases based on the intial model
scenariosInitialUC                       = 0.2-0.4*rand(2,Nscen) + repmat(thetaHat,1,Nscen);
scenariosInitialC                        = 0.4-0.8*rand(2,Nscen) + repmat(thetaHat,1,Nscen);
VappScenariosInitialUC                   = zeros(1,Nscen);
VappScenariosInitialC                    = zeros(1,Nscen);
h                                        = waitbar(0,'Running scenario evaluation of Vapp (based on initial model)...');
for ik = 1:Nscen
   VappScenariosInitialUC(end,ik)        =...
                                    Vapp(scenariosInitialUC(:,ik),thetaHat,yInitialUC,r,Ny,Nu,[inf,-inf],[inf,-inf]);
   VappScenariosInitialC(end,ik)         =...
                                    Vapp(scenariosInitialC(:,ik),thetaHat,yInitialC,r,Ny,Nu,Ub,Yb);  
   waitbar(ik/Nscen,h);
end
close(h)

VappScenariosMatrixInitialUC             = [scenariosInitialUC; VappScenariosInitialUC];
VappScenariosMatrixInitialC              = [scenariosInitialC; VappScenariosInitialC];
VappScenariosMatrixInitialUC(2,:)        = -VappScenariosMatrixInitialUC(2,:);
VappScenariosMatrixInitialC(2,:)         = -VappScenariosMatrixInitialC(2,:);

% ----- INPUT DESIGN USING INITIAL MODEL
% ------ Input design using ellipsoidal approximation
% The input design uses a transfer function model instead of state space model.
% The parameters of the two models are related as
%  [theta_tf(1); theta_tf(2)] = [theta_ss(1); -theta_ss(2)]
% The application cost Hessians need to be transformed accordingly in the design
%  VappHessian_tf = [1, 0; 0 -1]'*VappHessian_ss*[1, 0; 0 -1] 
%                 =: T*VappHessian_ss*T
T                                        = [1,0;0,-1];
optInputDesignUC                         = oidProblem(modelEstInitial,Nident,'MA',40);
optInputDesignUC.constraints{1}          = ...
                                    oidApplicationConstraint(T'*VappHessInitialUC*T,gamma,alpha);
optHUC                                   = solve(optInputDesignUC,[1 0 0]);

optInputDesignC                          = oidProblem(modelEstInitial,Nident,'MA',40);
optInputDesignC.constraints{1}           = ...
                                    oidApplicationConstraint(T'*VappHessInitialC*T,gamma,alpha);
optHC                                    = solve(optInputDesignC,[1 0 0]);

% ------ Input design using scenario approach
optInputDesignScenariosUC                = oidProblem(modelEstInitial,Nident,'MA',40);
optInputDesignScenariosUC.constraints{1} = ...
                                    oidApplicationConstraint(VappScenariosMatrixInitialUC(:,VappScenariosMatrixInitialUC(end,:)<=1/gamma),gamma,alpha);
optHScenariosUC                          = solve(optInputDesignScenariosUC, [1 0 0]);

optInputDesignScenariosC                 = oidProblem(modelEstInitial,Nident,'MA',40);
optInputDesignScenariosC.constraints{1}  = ...
                                    oidApplicationConstraint(VappScenariosMatrixInitialC(:,VappScenariosMatrixInitialC(end,:)<=1/gamma),gamma,alpha);
optHScenariosC                           = solve(optInputDesignScenariosUC, [1 0 0]);
%% MONTE-CARLO STUDY
NMC                 = 100;
thetaHatUC          = zeros(2,NMC);
thetaHatC           = zeros(2,NMC);
thetaHatScenariosUC = zeros(2,NMC);
thetaHatScenariosC  = zeros(2,NMC);
h                   = waitbar(0,'Running Monte-Carlo study...');
for tk = 1:NMC
   e                         = sqrt(le)*randn(Nident,1);
   u                         = lsim(optHUC,randn(Nident,1));
   y                         = lsim(trueSystem,u) + e;
   theta                     = getpvec(oe(iddata(y,u,Ts),[1,1,1],...
                                'MaxIter',100,'Tolerance',1e-10));
   thetaHatUC(:,tk)          = [theta(1); -theta(2)];
   
   u                         = lsim(optHC,randn(Nident,1));
   y                         = lsim(trueSystem,u) + e;
   theta                     = getpvec(oe(iddata(y,u,Ts),[1,1,1],...
                                'MaxIter',100,'Tolerance',1e-10));
   thetaHatC(:,tk)           = [theta(1); -theta(2)];
   
   u                         = lsim(optHScenariosUC,randn(Nident,1));
   y                         = lsim(trueSystem,u) + e;
   theta                     = getpvec(oe(iddata(y,u,Ts),[1,1,1],...
                                'MaxIter',100,'Tolerance',1e-10));
   thetaHatScenariosUC(:,tk) = [theta(1); -theta(2)];
   
   u                         = lsim(optHScenariosC,randn(Nident,1));
   y                         = lsim(trueSystem,u) + e;
   theta                     = getpvec(oe(iddata(y,u,Ts),[1,1,1],...
                                'MaxIter',100,'Tolerance',1e-10));
   thetaHatScenariosC(:,tk)  = [theta(1); -theta(2)];
   
   waitbar(tk/NMC,h);
end
close(h)
%% PLOT RESULTS
% ----- APPLICATION REGIONS
% ------- Unconstrained case
figure;
contour(theta0(1)+delta,theta0(2)+delta,VappNumUC',(1:10)/gamma); hold on
ellipse(gamma/2*VappHessUC,theta0,'r','empty');
plot(scenariosUC(1,VappScenariosUC<=(1/gamma)),scenariosUC(2,VappScenariosUC<=(1/gamma)),'x','Color',[1,0,0]);
ellipse(gamma/2*VappHessInitialUC,theta0,'b','empty');
plot(scenariosInitialUC(1,VappScenariosInitialUC<=(1/gamma))-thetaHat(1)+theta0(1),scenariosInitialUC(2,VappScenariosInitialUC<=(1/gamma))-thetaHat(2)+theta0(2),...
   'x','Color',[0,0,1]);
xlim([0.2 1])
xlabel('\theta_1')
ylabel('\theta_2')
legend('Application cost', 'Ellipsoidal approximation (true)', 'Scenarios (true)', 'Ellipsoidal approximation (initial)', 'Scenarios (initial)')
title('Unconstrained MPC')

% ------- Constrained case
figure;
contour(theta0(1)+delta,theta0(2)+delta,VappNumC',(1:10)/gamma); hold on
ellipse(gamma/2*VappHessC,theta0,'r','empty');
plot(scenariosC(1,VappScenariosC<=(1/gamma)),scenariosC(2,VappScenariosC<=(1/gamma)),'x','Color',[1,0,0]);
ellipse(gamma/2*VappHessInitialC,theta0,'b','empty');
plot(scenariosInitialC(1,VappScenariosInitialC<=(1/gamma))-thetaHat(1)+theta0(1),scenariosInitialC(2,VappScenariosInitialC<=(1/gamma))-thetaHat(2)+theta0(2),...
   'x','Color',[0,0,1]);
xlim([0.2 1])
xlabel('\theta_1')
ylabel('\theta_2')
legend('Application cost', 'Ellipsoidal approximation (true)', 'Scenarios (true)', 'Ellipsoidal approximation (initial)', 'Scenarios (initial)')
title('Constrained MPC')

% ----- SYSTEM IDENTIFICATION RESULTS USING INTIAL MODEL IN INPUT DESIGN
% ------- Unconstrained case
figure; hold on;
ellipse(gamma/2*VappHessUC,theta0,'r','empty');
ellipse(gamma/2*VappHessInitialUC,theta0,'b','empty');
ellipse(Nident/chi2inv(alpha,2)*(T\optInputDesignUC.informationMatrix/T),theta0,'k');
plot(thetaHatScenariosUC(1,:),thetaHatScenariosUC(2,:),'og');
plot(thetaHatUC(1,:),thetaHatUC(2,:),'xk');
xlim([0.5 0.7])
ylim([0.75,1.05])
xlabel('\theta_1')
ylabel('\theta_2')
legend('True application set', 'Estimated application set', 'Identification set','True parameters', 'Estimates (scenarios)', 'Estimates (ellispoidal)')
title('Unconstrained MPC')

% ------- Constrained case
figure; hold on;
ellipse(gamma/2*VappHessC,theta0,'r','empty');
ellipse(gamma/2*VappHessInitialC,theta0,'b','empty');
ellipse(Nident/chi2inv(alpha,2)*(T\optInputDesignC.informationMatrix/T),theta0,'k');
plot(thetaHatScenariosC(1,:),thetaHatScenariosC(2,:),'og');
plot(thetaHatC(1,:),thetaHatC(2,:),'xk');
xlim([0.25 0.9])
ylim([0.4 1.4])
xlabel('\theta_1')
ylabel('\theta_2')
legend('True application set', 'Estimated application set', 'Identification set','True parameters', 'Estimates (scenarios)', 'Estimates (ellipsoidal)')
title('Constrained MPC')
