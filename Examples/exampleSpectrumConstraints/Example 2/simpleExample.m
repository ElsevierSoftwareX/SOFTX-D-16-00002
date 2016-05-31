% Moose implementation of input design.
%
% The example illustrates the syntax used to set signal spectrum constraints (MIMO case).
%
% The example requires MOOSE2 and yalmip version 20130201 or above.

% Author: Mariette Annergren & Christian A. Larsson
% Copyright (c) 2015 Mariette Annergren & Christian A. Larsson

clear variables
close all
clc
P         = bodeoptions; P.MagUnits = 'abs'; P.FreqScale = 'linear';
wOut      = linspace(1e-5,pi,100);
%% SETUP THE PROBLEM 
% ----- SYSTEM
% All coefficients are to be estimated
Ts         = 1;  
numIn      = 2;
numOut     = 2;
A0{1,1}    = [1];       A0{1,2} = [0 1];          
A0{2,1}    = [0];       A0{2,2} = [1];                    
B0{1,1}    = [1 0.1];   B0{1,2} = [1 0.1];      
B0{2,1}    = [1 0.1];   B0{2,2} = [1 0.1];      
F0         = 1;
C0         = 1;
D0         = 1;
lambda     = eye(numOut);             
trueSystem = idpoly(A0,B0,C0,D0,F0,lambda,Ts);
theta0     = getpvec(trueSystem);
% ----- SYSTEM IDENTIFICATION MODEL AND EXPERIMENT
model      = trueSystem;
Nident     = 500;
alpha      = 0.95;
%% INPUT DESIGN USING MOOSE
% The optimal input design problem is specified. A 20 parameter FIR spectrum is used.
optInputDesign                                                  = oidProblem(model,Nident,'FIR',20);

% Excitation spectrum constraint -- they must be real otherwise they make no sense
%-- constraint is LTI
conA{1,1}                                                       = [1 0.9];   conA{1,2} = [1 0.9];           
conA{2,1}                                                       = [1 0.9];   conA{2,2} = [1 0.9];            
conB{1,1}                                                       = [1];       conB{1,2} = [1];      
conB{2,1}                                                       = [1];       conB{2,2} = [1]; 
conSys                                                          = tf(conA,conB,1,'Variable','z^-1'); 
% Note that the LTI lower / upper bound constraint in MIMO is treated as the linear matrix inequality Phi_u >= / <=  conSys*conSys'. 
%-- constraint is STRUCT
% [magconSys,wSys]                                                = freqresp(conSys*conSys',wOut);
% conSys                                                          = struct('B',real(magconSys),'w',wSys,'type',{'lmi'});
% Note that the STRUCT lower / upper bound constraint in MIMO is treated as the linear matrix or elementwise inequality Phi_u >= / <=  conSys. 
%-- constraint is DOUBLE
% conSys                                                          = [0.3 0.2; 0.2 0.3];
% Note that the DOUBLE lower / upper bound constraint in MIMO is treated as the linear matrix inequality Phi_u >= / <=  conSys. 
optInputDesign.spectrum.excitation.lb                           = conSys; 

% The input design problem is solved for the minimum input power and the
% corresponding optimal spectral factor is obtained.
[optH,yalmipInfo,infoMatrix,objValue,signalPower,spectrumParam] = solve(optInputDesign,[1 0 0]);

%% MONTE-CARLO STUDY
NMC        = 100;
thetaHatMC = zeros(length(theta0),NMC);
distTheta0 = zeros(1,NMC);
h          = waitbar(0,'Running Monte-Carlo study...');
for tk = 1:NMC
   e                = randn(Nident,size(trueSystem.NoiseVariance,2))*chol(lambda);
   u                = lsim(optH,randn(Nident,size(trueSystem,2)));
   y                = sim(trueSystem,[u e]);
   modelEst         = arx(iddata(y,u,Ts),model,...
                        'MaxIter',100,'Tolerance',1e-10);
   theta            = getpvec(modelEst);
   thetaHatMC(:,tk) = theta; 
   distTheta0(1,tk) = norm(thetaHatMC(:,tk)-theta0);
   waitbar(tk/NMC,h);
end
close(h)
%% PLOT RESULTS
plotResultsMIMO