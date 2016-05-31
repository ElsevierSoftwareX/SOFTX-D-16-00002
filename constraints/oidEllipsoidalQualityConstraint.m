% Representation of the ellipsoidal quality constraint
%
% F(w,theta) <= gamma_e for all w and theta in 
% 
% E = {theta : (theta-theta0)'R(theta-theta0) <= 1}, where
%
%               [W_nG(theta) + X_n]^* Y_n^*Y_n [W_nG(theta) + X_n] + K_n^*K_n
% F(w,theta) = -----------------------------------------------------------------, and R >= 0,
%               [W_dG(theta) + X_d]^* Y_d^*Y_d [W_dG(theta) + X_d] + K_d^*K_d
% where theta only paramterize the input/output model NOT the noise model.
%
% Construction:
% con = oidEllipsoidalQualityConstraint(Kd,Kn,Wd,Wn,Xd,Xn,Yd,Yn,gamma,conf,wSamp)
%
% Creates an ellipsoidal quality constraint. The arguments are
%   W_n,W_d,X_n,X_d:    Finite dimensional stable transfer functions. 
%   Y_n,Y_d,K_n,K_d:    Finite dimensional stable transfer functions.
%   gamma:              Defines the upper bound on F(w,theta). 
%   conf:               Confidence level of the estimates. conf is equal to the
%                       probability used in designing the system identification 
%                       ellipsoid.
%   wSamp:              A vector of frequency points where the constraint is evaluated
%                       or a string 'kyp' which uses the KYP lemma for the constraint.
% 
%  The ellipsoid E is the system identification set. The matrix R is set 
%  automatically to the properly scaled information matrix. 
% 
%  If wSamp is not specified, by default, a uniform frequency grid on [0,pi[ 
%  with 20 frequency points is used to evaluate the ellipsoidal quality constraint.
%  
%  The ellipsoidal quality constraint only supports fully parametrized SISO systems with
%  polynomial A and F.
%
%  Application constraints are added to the constraints property of
%  prob = oidProblem(...). The ith constraint is specified by
%  prob.constraints{i} = oidConstraint(...);
% 
% Theoretical reference (for example): Henrik Jansson, "Experiment design with
% applications in identification for control", PhD thesis, Section 3.7.


% Author: Mariette Annergren & Christian A. Larsson
% Copyright (c) 2015 Mariette Annergren & Christian A. Larsson

classdef oidEllipsoidalQualityConstraint < oidConstraint
   properties
      Kd;
      Kn;
      Xd;
      Xn;
      Wd;
      Wn;
      Yd;
      Yn;
      gamma;
      alpha;
      wSamp;
   end
   methods
      function obj = oidEllipsoidalQualityConstraint(varargin)
         % Order of arguments: Kd,Kn,Wd,Wn,Xd,Xn,Yd,Yn,gamma,alpha,wSamp
         if ~isempty(varargin{1}) && (isnumeric(varargin{1}) || isa(varargin{1},'lti'))
            obj.Kd = tf(varargin{1});
         else
            obj.Kd = 0;
         end
         if ~isempty(varargin{2}) && (isnumeric(varargin{2}) || isa(varargin{2},'lti'))
            obj.Kn = tf(varargin{2});
         else
            obj.Kn = 0;
         end
         if ~isempty(varargin{3}) && (isnumeric(varargin{3}) || isa(varargin{3},'lti'))
            obj.Wd = tf(varargin{3});
         else
            obj.Wd = 0;
         end
         if ~isempty(varargin{4}) && (isnumeric(varargin{4}) || isa(varargin{4},'lti'))
            obj.Wn = tf(varargin{4});
         else
            obj.Wn = 0;
         end
         if ~isempty(varargin{5}) && (isnumeric(varargin{5}) || isa(varargin{5},'lti'))
            obj.Xd = tf(varargin{5});
         else
            obj.Xd = 0;
         end
         if ~isempty(varargin{6}) && (isnumeric(varargin{6}) || isa(varargin{6},'lti'))
            obj.Xn = absorbDelay(tf(varargin{6}));
         else
            obj.Xn = 0;
         end
         if ~isempty(varargin{7}) && (isnumeric(varargin{7}) || isa(varargin{7},'lti'))
            obj.Yd = tf(varargin{7});
         else
            obj.Yd = 0;
         end
         if ~isempty(varargin{8}) && (isnumeric(varargin{8}) || isa(varargin{8},'lti'))
            obj.Yn = tf(varargin{8});
         else
            obj.Yn = 0;
         end
         obj.gamma = varargin{9};
         obj.alpha = varargin{10};
         if nargin > 10
             if strcmpi(varargin{11},'kyp')
                 obj.wSamp = [];
             elseif ~isempty(varargin{11})
                 obj.wSamp = varargin{11};
             else
                 obj.wSamp = linspace(0,pi,20);
             end
         else
              obj.wSamp = linspace(0,pi,20);
         end
      end
      function con = getConstraints(obj,prob)
         con = [];
         
         if ~isa(prob.model_.model,'idpoly')
            throw(MException('VerifyEllipsoidalQualityConstraints:wrongModelType','Ellipsoidal constraints only handles idpoly models'));
         end
         
         na = prob.model_.model.na;
         nb = prob.model_.model.nb;
         nc = prob.model_.model.nc;
         nd = prob.model_.model.nd;
         nf = prob.model_.model.nf;
         nk = prob.model_.model.nk;
         
         % Model check
         if (na && nf) || size(prob.model_,1) > 1 || size(prob.model_,2) > 1
            throw(MException('VerifyEllipsoidalQualityConstraints:wrongModelStructure','Ellipsoidal constraints only handles SISO models with A or F polynomial'));
         end
         
         % Delay vectors
         zi =tf([0 1],1,prob.model_.Ts,'variable','z^-1');
         ZD = tf(zeros(1,na+nf+nb));
         for ik = 1:(na+nf)
            ZD(ik)=zi^ik;
         end
         ZN = tf(zeros(1,na+nf+nb));
         for ik = 1:nb
            ZN((na+nf)+ik)=zi^(ik-1+nk);
         end
         
         % Reduced information matrix for G parameters, permute of F polynomial
         T = [eye(na), zeros(na,nb+nc+nd+nf); zeros(nf,na+nb+nc+nd), eye(nf);...
            zeros(nb,na), eye(nb), zeros(nb,nf+nc+nd)];
         th = T*getpvec(prob.model_.model); % eta_0 in Jansson's thesis
         iF = T*prob.N*prob.spectrum.informationMatrix(prob)/chi2inv(obj.alpha,numel(th))*T';
         R = [-iF, iF*th; th'*iF, 1 - th'*iF*th];
         
         % Construct the F0 and F1 functions
         % Least common denominator
         lcd = tf(poly(pole([obj.Xn obj.Xd obj.Yn obj.Yd obj.Kn obj.Kd obj.Wn obj.Wd])),1,prob.model_.Ts);
         xn = [obj.Wn, obj.Xn, obj.Xn]';
         xd = [obj.Wd, obj.Xd, obj.Xd]';
         zz = [0; 1; 1];
         ZV = [ZN, 0; ZD, 0; zeros(size(ZN)), 1];
         
         Mn   = minreal(lcd'*ZV'*(obj.Yn'*obj.Yn)*(xn*xn')*ZV*lcd + lcd'*ZV'*(obj.Kn'*obj.Kn)*(zz*zz')*ZV*lcd);
         F1tf = 0.5*(Mn + (Mn').');
         F1   = (F1tf);
         
         Md   = minreal(lcd'*ZV'*(obj.Yd'*obj.Yd)*(xd*xd')*ZV*lcd + lcd'*ZV'*(obj.Kd'*obj.Kd)*(zz*zz')*ZV*lcd);
         F0tf = 0.5*(Md + (Md').');
         F0   = (F0tf);
         
         if isempty(obj.wSamp)
             % Get coefficients in F0 and F1
             [numF0tf,~]   = tfdata(F0tf);
             nCoeffF0tf    = ceil(cellfun('length',numF0tf)/2);
             [numF1tf,~]   = tfdata(F1tf);
             nCoeffF1tf    = ceil(cellfun('length',numF1tf)/2);
             nCoeffMaxF0F1 = max([max(max(nCoeffF0tf)) max(max(nCoeffF1tf))]);                
             % FIR spectrum formulation of tau
             nCoeffMaxTau  = 20; % fixed value
             tau           = sdpvar(nCoeffMaxTau,1);
             if nCoeffMaxTau < 2
                 Atau = [];
                 Btau = [];
                 Ctau = [];
                 Dtau = 0.5*tau;
             else
                 Atau = [zeros(1,nCoeffMaxTau-2) zeros(1,1); eye(nCoeffMaxTau-2) zeros(nCoeffMaxTau-2,1)];
                 Btau = [1 zeros(1,nCoeffMaxTau-2)]';
                 Ctau = tau(2:end,1)';
                 Dtau = 0.5*tau(1);
                 ntau = size(Atau,1);
                 Qtau = sdpvar(ntau,ntau,'symmetric');
                 con  = [con;[Qtau-Atau'*Qtau*Atau, Ctau'-Atau'*Qtau*Btau; Ctau - Btau'*Qtau*Atau, Dtau+Dtau'-Btau'*Qtau*Btau] >= 0];
             end
             % Get coefficients of F0 and F1
             coeffF0tf = zeros(size(F0tf,1),size(F0tf,2),nCoeffMaxF0F1);
             coeffF1tf = zeros(size(F1tf,1),size(F1tf,2),nCoeffMaxF0F1);     
             for row = 1:size(F0tf,1)
                 for col = 1:size(F0tf,2) % however, the matrix is square and size(F0tf) == size(F1tf)
                     coeffF0tf(row,col,1:nCoeffF0tf(row,col)) = numF0tf{row,col}(nCoeffF0tf(row,col):end);
                     coeffF1tf(row,col,1:nCoeffF1tf(row,col)) = numF1tf{row,col}(nCoeffF1tf(row,col):end);
                 end
             end
             coeffGammaF0F1 = obj.gamma*coeffF0tf-coeffF1tf; % tensor: [dim(eta)+1 x dim(eta)+1 x number of lags]
             % FIR spectrum formulation of Lambda
             % Construct the positive real part of Lambda
             nCoeffMaxLambda = nCoeffMaxTau+nCoeffMaxF0F1-1;
             Lambda = [];
             for index = 0:nCoeffMaxLambda-1
                 switchLag = 1;
                 for coeffTau = -(nCoeffMaxTau-1):nCoeffMaxTau-1
                     for coeffF0F1 = -(nCoeffMaxF0F1-1):nCoeffMaxF0F1-1
                         if (coeffTau+coeffF0F1)==index
                             if switchLag
                                 if isempty(Lambda)
                                     Lambda    = -R;
                                     switchLag = 0;
                                 else
                                     Lambda = cat(3,Lambda,tau(abs(coeffTau)+1)*coeffGammaF0F1(:,:,abs(coeffF0F1)+1));
                                     switchLag = 0;
                                 end
                             else
                                 if index == 0
                                     Lambda = Lambda + tau(abs(coeffTau)+1)*coeffGammaF0F1(:,:,abs(coeffF0F1)+1);
                                 else
                                     Lambda(:,:,index+1) = Lambda(:,:,index+1) + tau(abs(coeffTau)+1)*coeffGammaF0F1(:,:,abs(coeffF0F1)+1);
                                 end
                             end
                         end
                     end
                 end
             end
             if nCoeffMaxLambda < 2
                 ALambda = [];
                 BLambda = [];
                 CLambda = [];
                 DLambda = 0.5*Lambda;
                 con     = [con; Lambda >= 0];
             else
                 ALambda = zeros(size(F0tf,1)*(nCoeffMaxLambda-1));
                 ALambda(size(F0tf,1)+1:size(F0tf,1)*(nCoeffMaxLambda-1),1:size(F0tf,1)*(nCoeffMaxLambda-2)) = eye(size(F0tf,1)*(nCoeffMaxLambda-2));
                 BLambda = [eye(size(F0tf,1)); zeros(size(F0tf,1)*(nCoeffMaxLambda-2),size(F0tf,1))];
                 CLambda = reshape(Lambda(:,:,2:end),size(F0tf,1),size(F0tf,1)*(nCoeffMaxLambda-1));
                 DLambda = 0.5*Lambda(:,:,1);
                 nLambda = size(ALambda,1);
                 QLambda = sdpvar(nLambda,nLambda,'symmetric');
                 con  = [con;[QLambda-ALambda'*QLambda*ALambda, CLambda'-ALambda'*QLambda*BLambda; CLambda - BLambda'*QLambda*ALambda, DLambda+DLambda'-BLambda'*QLambda*BLambda] >= 0];
             end
             
         else
             F0fr = real(freqresp(F0,obj.wSamp/prob.model_.Ts));
             F1fr = real(freqresp(F1,obj.wSamp/prob.model_.Ts));
             tau  = sdpvar(numel(obj.wSamp),1);
             
             for ik = 1:numel(obj.wSamp)
                 eigVec = eig(obj.gamma*F0fr(:,:,ik) - F1fr(:,:,ik)) ;
                 if sum(eigVec<0)==0;
                     warningOutput = warning('on','oidEllipsoidalQualityConstraint:positiveSemiDefinite');
                     warning('oidEllipsoidalQualityConstraint:positiveSemiDefinite',...
                         'Ellipsoidal Quality Constraint:gamma*F0 - F1 should not be positive semidefinite');
                     warning(warningOutput)
                 end
                 con = [con; tau(ik)*(obj.gamma*F0fr(:,:,ik) - F1fr(:,:,ik)) - R > 0];
             end
             con = [con; tau>0];
         end
      end
   end
end