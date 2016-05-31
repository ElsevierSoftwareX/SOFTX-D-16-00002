% Representation of signal spectrum for input design
%
% The oidSpectrum class is the base class for all parametric spectra in the
% MOOSE2 input design toolbox. The class is not user-interfacing and cannot be
% instantiated. It is only used to define the interface of parametric spectra.
%
% See also oidSpectrumFDP, oidSpectrumMA


% Author: Mariette Annergren & Christian A. Larsson
% Copyright (c) 2015 Mariette Annergren & Christian A. Larsson

classdef (Abstract=true) oidSpectrum  
   properties
      type;
   end
   properties (Dependent = true, Abstract=true)  
      M;
   end
   methods
      function obj = oidSpectrum(type)
         obj.type = type;
      end
   end
   methods (Abstract=true)
      con = getConstraints(obj,cons);
      iF  = informationMatrix(obj,model);
      Pr  = excitationPower(obj);
      Pu  = inputPower(obj);
      Py  = outputPower(obj,model);
   end
end