% Representation of transfer function models


% Author: Mariette Annergren & Christian A. Larsson
% Copyright (c) 2015 Mariette Annergren & Christian A. Larsson

classdef oidModelPoly < oidModel
   properties
      model;
   end
   properties (Dependent=true)
      NoiseVariance;
      Ts;
   end
   methods
      function obj = oidModelPoly(model)
         if isa(model,'idpoly')
            obj.model = model;
         else
            throw(MException('VerifyModel:incompatibleModel',['Model should '...
               'be specified as a idpoly object']));
         end
      end
      function [G,H] = tf(obj) 
         G = absorbDelay(tf(obj.model,'measured'));
         H = absorbDelay(minreal(tf(obj.model,'noise')/chol(obj.model.NoiseVariance)));
      end
      function [dG,dH] = gradient(obj)  
         % Gradients of the transfer functions w.r.t. the parameter vector.
         % The parameter numbering follows the numbering of the System
         % Identification Toolbox, i.e., A,B,C,D,F.
         n = nparams(obj.model,'free');
         [p,m] = size(obj.model);
         Ts = obj.model.Ts;
         
         % Number of parameters
         na = obj.model.na;
         nb = obj.model.nb;
         nc = obj.model.nc;
         nd = obj.model.nd;
         nf = obj.model.nf;
         
         if p > 1
            C = num2cell(zeros(p));
            D = num2cell(ones(p));
            for i = 1:p
               C{i,i} = obj.model.c{i};
               D{i,i} = obj.model.d{i};
            end
         else
            C = obj.model.c;
            D = obj.model.d;
         end
         
         % Transfer functions
         A = tf(obj.model.A,1,Ts,'variable','z^-1');
         Ai = minreal(eye(p)/A);
         BF = tf(obj.model.B,obj.model.F,Ts,'variable','z^-1');
         CD = tf(C,D,Ts,'variable','z^-1');
         D = tf(D,1,Ts,'variable','z^-1');
         F = tf(obj.model.F,1,Ts,'variable','z^-1');
         
         % Allocate memory for derivatives
         dG = cell(n,1);
         dH = cell(n,1);
         
         dA = tf(zeros(p,p,n));
         dB = tf(zeros(p,m,n));
         dC = tf(zeros(p,p,n));
         dD = tf(zeros(p,p,n));
         dF = tf(zeros(p,m,n));
         
         zi =tf([0 1],1,Ts,'variable','z^-1');
         
         % Derivatives of transfer functions
         % Loop over outputs
         paramOffset = 0;
         for i = 1:p
            % dA
            for j = 1:p
               for k = 1:na(i,j)
                  if obj.model.Structure.a(i,j).Free(k+1)
                     dA(i,j,k+paramOffset) = zi^k;
                  end
               end
               if ~isempty(k)
                  paramOffset = paramOffset + k;
               end
            end
            
            % dB
            for j = 1:m
               for k = 1:nb(i,j)
                  if obj.model.Structure.b(i,j).Free(k+obj.model.nk(i,j))
                     dB(i,j,k+paramOffset) =...
                        minreal(zi^(k+obj.model.nk(i,j)-1)/F(i,j));
                  end
               end
               if ~isempty(k)
                  paramOffset = paramOffset + k;
               end
            end
            
            % dF
            for j = 1:m
               for k = 1:nf(i,j)
                  if obj.model.Structure.f(i,j).Free(k+1)
                     dF(i,j,k+paramOffset) = minreal(-zi^k*BF(i,j)/F(i,j));
                  end
               end
               if ~isempty(k)
                  paramOffset = paramOffset + k;
               end
            end
            
            % dC
            for k = 1:nc(i)
               if obj.model.Structure.c(i,1).Free(k+1);
                  dC(i,i,k+paramOffset) = minreal(zi^k/D(i,i));
               end
            end
            if ~isempty(k)
               paramOffset = paramOffset + k;
            end
            
            % dD
            for k = 1:nd(i)
               if obj.model.Structure.d(i,1).Free(k+1);
                  dD(i,i,k+paramOffset) = minreal(zi^k*CD(i,i)/D(i,i));
               end
            end
            if ~isempty(k)
               paramOffset = paramOffset + k;
            end
         end
         
         for k = 1:n;
            dAi = minreal(-Ai*dA(:,:,k)*Ai);
            dG{k} = minreal(dAi*BF + Ai*(dB(:,:,k) + dF(:,:,k)));
            dH{k} = minreal(dAi*CD + Ai*(dC(:,:,k) + dD(:,:,k)));
         end
      end
      function varargout = size(obj,varargin)  
         varargout = cell(nargout,1);
         if nargout > 0
            if nargin > 1
               varargout{1} = size(obj.model,varargin{1});
            else
               for ik = 1:nargout
                  varargout{ik} = size(obj.model,ik);
               end
            end
         else
            if nargin > 1
               size(obj.model,varargin{1})
            else
               size(obj.model)
            end
         end
      end
      function n = nparams(obj)  
         n = nparams(obj.model,'free');
      end
      function noiseVar = get.NoiseVariance(obj)
         noiseVar = obj.model.NoiseVariance;
      end
      function Ts = get.Ts(obj)
         Ts = obj.model.Ts;
      end
      function obj = set.NoiseVariance(obj,noiseVar)
         obj.model.NoiseVariance = noiseVar;
      end
      function obj = set.Ts(obj,Ts)
         obj.model.Ts = Ts;
      end
   end
end