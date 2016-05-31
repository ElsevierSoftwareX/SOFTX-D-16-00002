% Representation of the application constraint
%
% Vapp(theta) <= 1/gamma   w.p. alpha 
% 
% Construction:
% con = oidApplicationConstraint(V,gamma,conf)
%
% Creates an application constraint. The arguments are
%   V:      V is either the Hessian of Vapp (VappHessian) or a collection of
%           scenarios from Vapp (VappScenarios). The class automatically
%           differentiates between the two by a check of dimensions. If V is a
%           square matrix then V = VappHessian. If V is a fat matrix then V =
%           VappScenarios.
%   gamma:  Defines the upper bound on the allowed performance degradation. 
%   conf:   Confidence level of the estimates. conf is equal to the
%           probability used in designing the system identification ellipsoid.
% 
%  Application constraints are added to the constraints property of
%  prob = oidProblem(...). The ith constraint is specified by
%  prob.constraints{i} = oidConstraint(...);
% 
% Theoretical references (for example): 
% 1. X. Bombois, G. Scorletti, M. Gevers, P. M. J. V. D. Hof, and R. Hildebrand. 
%    Least costly identification experiment for control. Automatica, 42:1651–1662, 2006.
% 2. H. Hjalmarsson. System identification of complex and structured systems. European
%    Journal of Control, 15(34):275 – 310, 2009.


% Author: Mariette Annergren & Christian A. Larsson
% Copyright (c) 2015 Mariette Annergren & Christian A. Larsson

classdef oidApplicationConstraint < oidConstraint
   properties
      VappScenarios % Scenarios of application cost 
                    % rows - elements of parameter, last row - value of Vapp, columns - scenarios, 
      VappHessian   % Hessian of application cost
      gamma         % Upper bound on application cost
      conf          % Confidence level
   end
   methods
      function obj = oidApplicationConstraint(varargin)
         obj = obj@oidConstraint;
         if nargin > 0
             if size(varargin{1},1)< size(varargin{1},2)
                 obj.VappScenarios = varargin{1};
                 obj.VappHessian   = [];
             elseif size(varargin{1},1) == size(varargin{1},2)
                 obj.VappHessian   = varargin{1};
                 obj.VappScenarios = [];
             else
                 throw(MException('VerifyApplicationConstraints:wrongDimensions','Application Constraint:VappScenarios or VappHessian cannot be a tall matrix'));
             end
         end
         if nargin == 3
            obj.gamma = varargin{2};
            obj.conf  = varargin{3};
         end
      end
      function con = getConstraints(obj,problem)
          con = [];
          if ~isempty(obj.VappHessian)
              con = 2*problem.spectrum.informationMatrix(problem) >= ...
                  chi2inv(obj.conf,nparams(problem.model,'free'))*...
                  obj.gamma*obj.VappHessian/problem.N;
          else
              allParam   = getpar(problem.model,'value');
              freeIndex  = getpar(problem.model,'free');
              theta0     = allParam(freeIndex);
              infoMatrix = problem.spectrum.informationMatrix(problem);
              for i = 1:size(obj.VappScenarios,2)
                  con = [con; (theta0-obj.VappScenarios(1:end-1,i))'*infoMatrix*...
                            problem.N*(theta0-obj.VappScenarios(1:end-1,i)) >=...
                            chi2inv(obj.conf,nparams(problem.model,'free'))*obj.gamma*obj.VappScenarios(end,i)];                      
              end
          end
      end
   end
end