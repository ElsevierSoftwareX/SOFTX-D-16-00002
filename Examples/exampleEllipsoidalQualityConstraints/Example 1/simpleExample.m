% Moose implementation of example based on Example 3.16 in Henrik Jansson, 
% "Experiment design with applications in identification for control", PhD thesis.
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
a          = -0.7;
b          = 0.36;
theta0     = [b a]; 
A0         = 1;
B0         = [0 theta0(1)];
F0         = [1 theta0(2)];
C0         = 1;
D0         = 1;
lambda     = 0.1;             
trueSystem = idpoly(A0,B0,C0,D0,F0,lambda,Ts);
% ----- SYSTEM IDENTIFICATION MODEL AND EXPERIMENT
model      = trueSystem;
G0         = tf(B0,F0,Ts,'Variable','z^-1');       
Nident     = 500;
alpha      = 0.95; 
% ----- ELLIPSOIDAL QUALITY CONSTRAINT
T          = tf([0 (1-0.1353)^2],[1 -2*0.1353 0.1353^2],Ts,'Variable','z^-1');
Wn         = 1; 
Wd         = 1; 
Xn         = -G0; 
Yn         = T;
Yd         = 1; 
Xd         = 0; 
Kn         = 0; 
Kd         = 0; 
gamma      = 0.1;          
wSamp      = linspace(0,pi,30);
% ------- Numerical evaluation of ellipsoidal quality constraint for contours
Nnum       = 150;
Fnum       = zeros(Nnum);
delta      = linspace(-0.1,0.1,Nnum);
Tfr        = squeeze(freqresp(T,wSamp)); 
Tmag       = Tfr'.*Tfr.';
% h          = waitbar(0,'Running contour evaluation of ellipsoidal quality constraint...');
% for ik = 1:Nnum
%    for jk = 1:Nnum
%       Fnum(ik,jk) = F(theta0+[delta(ik) delta(jk)],G0,Tmag,wSamp,Ts);
%    end
%    waitbar(ik/Nnum,h);
% end
% close(h)
% save Fnum
load Fnum
%% INPUT DESIGN USING MOOSE (with frequency grid)
% The optimal input design problem is specified. A 20 parameter FIR
% spectrum is used.
optInputDesignFreq                = oidProblem(model,Nident,'FIR',20);

% The application contraint is added to the input design problem 
optInputDesignFreq.constraints{1} = oidEllipsoidalQualityConstraint(Kd,Kn,Wd,Wn,Xd,Xn,Yd,Yn,gamma^2,alpha,wSamp);

% The input design problem is solved for the minimum input power and the
% corresponding optimal spectral factor is obtained.
[optHFreq,~,~,optPowerFreq]       = solve(optInputDesignFreq, [1 0 0]);
%% INPUT DESIGN USING MOOSE (with KYP)
% The optimal input design problem is specified. A 20 parameter FIR
% spectrum is used.
optInputDesignKYP                = oidProblem(model,Nident,'FIR',20);

% The application contraint is added to the input design problem 
optInputDesignKYP.constraints{1} = oidEllipsoidalQualityConstraint(Kd,Kn,Wd,Wn,Xd,Xn,Yd,Yn,gamma^2,alpha,'KYP');

% The input design problem is solved for the minimum input power and the
% corresponding optimal spectral factor is obtained.
[optHKYP,~,~,optPowerKYP]        = solve(optInputDesignKYP, [1 0 0]);

%% MONTE-CARLO STUDY
NMC            = 100;
thetaHatMCFreq = zeros(2,NMC);
h              = waitbar(0,'Running Monte-Carlo study for solution using frequency grid...');
for tk = 1:NMC
   e                    = sqrt(lambda)*randn(Nident,1);
   u                    = lsim(optHFreq,randn(Nident,1));
   y                    = lsim(trueSystem,u) + e;
   theta                = getpvec(oe(iddata(y,u,Ts),model,...
                          'MaxIter',100,'Tolerance',1e-10));
   thetaHatMCFreq(:,tk) = [theta(1);theta(2)]; 
   waitbar(tk/NMC,h);
end
close(h)
thetaHatMCKYP = zeros(2,NMC);
h             = waitbar(0,'Running Monte-Carlo study for solution using KYP...');
for tk = 1:NMC
   e                   = sqrt(lambda)*randn(Nident,1);
   u                   = lsim(optHKYP,randn(Nident,1));
   y                   = lsim(trueSystem,u) + e;
   theta               = getpvec(oe(iddata(y,u,Ts),model,...
                         'MaxIter',100,'Tolerance',1e-10));
   thetaHatMCKYP(:,tk) = [theta(1);theta(2)]; 
   waitbar(tk/NMC,h);
end
close(h)
%% PLOT RESULTS
figure;
hold on
ellipse(Nident/chi2inv(alpha,2)*(optInputDesignFreq.informationMatrix),theta0,'k');
ellipse(Nident/chi2inv(alpha,2)*(optInputDesignKYP.informationMatrix),theta0,'r','empty');
plot(thetaHatMCFreq(1,:),thetaHatMCFreq(2,:),'xk');
plot(thetaHatMCKYP(1,:),thetaHatMCKYP(2,:),'xr');
contour(theta0(1)+delta,theta0(2)+delta,Fnum',gamma^2,'ShowText','on');
xlabel('\theta_1')
ylabel('\theta_2')
legend('Identification set (freq)', 'True parameters','Identification set (KYP)','Estimates (freq)','Estimates (KYP)','Countours of |\Delta|_{\infty}\leq\gamma')

