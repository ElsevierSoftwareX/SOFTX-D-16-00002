% Representation of models for input design
%
%  The oidModel class is the base class for all parametric models in the
%  MOOSE2 input design toolbox. The class is not user-interfacing and cannot be
%  instantiated. It is only used to define the interface of parametric models.
%
% See also oidModelPoly


% Author: Mariette Annergren & Christian A. Larsson
% Copyright (c) 2015 Mariette Annergren & Christian A. Larsson

classdef (Abstract=true) oidModel  
   properties (Abstract=true) 
      model;
      NoiseVariance;
      Ts;
   end
   methods (Abstract=true) 
      n       = nparams(obj);
      sz      = size(obj);
      [G,H]   = tf(obj);
      [dG,dH] = gradient(obj);
   end
end