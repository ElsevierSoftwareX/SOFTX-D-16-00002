% Representation of constraints for input design. 
%
% The oidConstraint class is the base class for all constraints (except 
% spectrum related ones) in the MOOSE2 input design toolbox. The class is
% not user-interfacing and cannot be instantiated. It is only used to 
% define the interface of constraints.
%
% See also oidApplicationConstraint, oidEllipsoidalConstraint,
% oidQualityConstraint


% Author: Mariette Annergren & Christian A. Larsson
% Copyright (c) 2015 Mariette Annergren & Christian A. Larsson

classdef oidConstraint
   properties
   end
   methods
      function obj = oidConstraint()
      end
   end
   methods (Abstract)
      con = getConstraints(obj,problem);
   end
end