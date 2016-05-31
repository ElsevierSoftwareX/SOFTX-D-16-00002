% Representation of the quality (weighted trace) constraint
%
% (1): trace W(w)P        <= gamma   for all w
%       W(w) = V(w)V(w)^* >= 0       for all w
%
% Two quality constrains that can be represented as the above constraint are
% (2): norm(|T(e^iw)/G_0(e^iw)|^2 dG^*(e^iw,theta_0) P dG(e^iw,theta_0),2)<=gamma
% (3): norm(|T(e^iw)/G_0(e^iw)|^2 dG^*(e^iw,theta_0) P dG(e^iw,theta_0),inf)<=gamma
%
% Construction
% con = oidQualityConstraint(V) defines constraint (1) with gamma = 1.
%
% con = oidQualityConstraint(V,gamma) defines constraint (1) for arbitrary gamma >= 0.
%
% con = oidQualityConstraint(V,gamma,wSamp) gives the option to change how the
% infinite constraint (1) is transformed to a finite dimensional constraint.
% wSamp can be either: 1. A vector of frequency points where (1) is evaluated.
%                      2. The string 'kyp', which uses the KYP lemma for (1).
% By default, a uniform frequency grid on [0,pi[ is used with 20 frequency
% points. This is found to be much more numerically robust than using KYP.
%
% For SISO systems
% con = oidQualityConstraint(V,gamma,wSamp,norm) defines constraint (2) if 
% norm = 2 and constraint (3) if norm = inf, for arbitrary gamma => 0. Note
% that the argument V in these cases defines T. Also note that when
% evaluating constraint (2), neither KYP nor a frequency grid is used.
%
% Theoretical reference (for example): Henrik Jansson, "Experiment design with
% applications in identification for control", PhD thesis, Section 3.6.


% Author: Mariette Annergren & Christian A. Larsson
% Copyright (c) 2015 Mariette Annergren & Christian A. Larsson

classdef oidQualityConstraint < oidConstraint
   properties
      V        % Spectral factor
      wSamp    % Sampling points on frequency axis or string kyp
      gamma    % Scaling
      norm     % For easy handling of constraints (2) and (3)
      scaling  % Scaling factor for numerical conditioning, default = 1 (not fully implemented yet)
   end
   methods
      function obj = oidQualityConstraint(varargin)
         if nargin > 0
            obj.V = varargin{1};
         end
         if nargin > 1 &&  ~isempty(varargin{2})
            obj.gamma = varargin{2};
         else
            obj.gamma = 1;
         end
         if nargin > 2
            if strcmpi(varargin{3},'kyp')
               obj.wSamp = [];
            elseif ~isempty(varargin{3})
               obj.wSamp = varargin{3};
            else
               obj.wSamp = linspace(0,pi,20);
            end
         end
         if nargin > 3
            if varargin{4} == 2
               obj.norm = 2;
            elseif isinf(varargin{4})
               obj.norm = inf;
            else
               obj.norm = [];
            end
         end
         obj.scaling = 1;
      end
            
      function con = getConstraints(obj,prob)
         con = [];
         iF = prob.N*prob.spectrum.informationMatrix(prob);
         n = prob.model_.nparams();
         
         if obj.norm == 2
            if size(prob.model_,1) == 1 && size(prob.model_,2) == 1
               [G,~] = prob.model_.tf();
               [dGcell,~] = prob.model_.gradient();
               dG = tf(zeros(n,1));
               for ik = 1:n
                  dG(ik) = dGcell{ik};
               end
               W = covar(minreal(obj.V/G)*dG,1);
               if W > 1e3;
                   warningOutput = warning('on','oidQualityConstraints:numericalProblems');
                   warning('oidQualityConstraints:numericalProblems','Quality Constraint:Numerical issues encountered when evaluating quality constraint. Is V/G*dG stable?');
                   warning(warningOutput)
               end
               V = chol(W)';
            else
               throw(MException('oidQualityConstraints:wrongConstraint','Quality Constraint:Norm constraint only for SISO'));
            end
            Z = sdpvar(size(V,2));
            con = [[Z, V'; V, iF] >= 0];
         else
            if isinf(obj.norm)
               if size(prob.model_,1) == 1 && size(prob.model_,2) == 1
                  [G,~] = prob.model_.tf();
                  [dGcell,~] = prob.model_.gradient();
                  dG = tf(zeros(n,1));
                  for ik = 1:n
                     dG(ik) = dGcell{ik};
                  end
                  V = minreal(ss(minreal(obj.V/G*dG)),[],0);
               else
               throw(MException('oidQualityConstraints:wrongConstraint','Quality Constraint:Norm constraint only for SISO'));
               end
            else
               V = minreal(ss(obj.V),[],0);
            end
            
            Z = sdpvar(size(V,2));
            if isempty(obj.wSamp)     
               [Am,Bm,Cm,Dm] = ssdata(V);
               nM = size(Am,1);
               Bm = [Bm zeros(nM,n)];
               Cm = [zeros(1,nM); Cm];
               DDm = [Z/obj.scaling, Dm'; Dm iF*obj.scaling];
               
               Q = sdpvar(nM,nM);
               con = [con; [Q-Am'*Q*Am, Cm'-Am'*Q*Bm; Cm - Bm'*Q*Am, DDm-Bm'*Q*Bm] >= 0];
            else
               Vfr = freqresp(V,obj.wSamp/V.Ts);
               for ik = 1:numel(obj.wSamp)
                  con = [con; [Z, Vfr(:,:,ik)'; Vfr(:,:,ik), iF] >= 0];
               end
            end
         end
         con = [con; obj.gamma - trace(Z) >= 0];
      end
   end
end