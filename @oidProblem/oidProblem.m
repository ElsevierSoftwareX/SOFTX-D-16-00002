% Representation of the optimal input design problem
%
% minimize          objective function
% with respect to   Phi_r(w)
% subject to        constraints
% 
% The objective function can be:
%  1. a_r*excitationPower + a_u*inputPower + a_y*outputPower
%
%  2. a_r*excitationPower + a_u_r*inputPower_r + a_u_e*inputPower_e
%                           + a_y_r*outputPower_r + a_y_e*outputPower_e
%
%  3. trace(P), where P is the covariance matrix of the estimated parameters
%
%  4. max(eig((P))), where P is the covariance matrix of the estimated parameters
%
%  5. det(P), where P is the covariance matrix of the estimated parameters
%
% The constraints can be any number of combinations of:
%  1. Spectrum constraints:
%     r_lb(w) <= Phi_r(w) <= r_ub(w)
%     u_lb(w) <= Phi_u(w) <= u_ub(w)
%     y_lb(w) <= Phi_y(w) <= y_ub(w)
%
%  2. Power constraints:
%     rP_lb(w) <= excitationPower  <= rP_ub
%     uP_lb(w) <= inputPower       <= uP_ub
%     yP_lb(w) <= outputPower      <= yP_ub
%
%  3. Application constraints
%     Vapp(theta) <= 1/gamma_a   w.p. alpha
%
%  4. Quality (weighted trace) constraints
%     trace W(w)P       <= gamma_q  for all w
%     W(w) = V(w)V(w)^* >= 0        for all w
%
%  5. Ellipsoidal quality constraints
%     F(w,e) <= gamma_e for all w and e in E = {e : (e-e0)'R(e-e0) <= 1}, where
%
%              [W_nG(e) + X_n]^* Y_n [W_nG(e) + X_n] + K_n
%     F(w,e) = --------------------------------------------, and R >= 0.
%              [W_dG(e) + X_d]^* Y_d [W_dG(e) + X_d] + K_d
%
% Construction:
%  prob = oidProblem(MODEL)
%  prob = oidProblem(MODEL,N)
%  prob = oidProblem(MODEL,N,type,M)
%  prob = oidProblem(MODEL,N,type,M,C,cParam)
%
%  Creates an input design problem for the model MODEL. The arguments are
%     MODEL:  Model of system to design input for, idpoly object.
%     N:      Experiment length.
%     type:   Input spectrum type, 'AR'            - Autoregressive (not yet implemented),
%                                  'MA' (or 'FIR') - Moving average,
%                                  'LA'            - Laguerre parameterization (not yet implemented),
%                                  'MS'            - Multisine (not yet implemented),
%     M:      Number of parameters in spectrum.
%     C:      Controller for closed loop design, a transfer function for the
%              fixed controller case and the order of the controller for the
%              free controller case (the latter not yet implemented).
%     cParam: Set to 'fixed' for a fixed controller or 'free' to include the
%               controller parameters in the design (the latter not yet implemented).
%     Default values are N=1, 'type'='MA', M=20, C=0, 'cParam'='fixed'.
%
% Adding constraints:
%  Spectrum and power constraints are specified directly in the spectrum
%  property of prob = oidProblem(...).
%   prob.spectrum.signal.ub:        Frequency-wise upper bound on spectrum
%   prob.spectrum.signal.lb:        Frequency-wise lower bound on spectrum
%   prob.spectrum.signal.power.ub:  Upper bound on total power constraint (scalar)
%   prob.spectrum.signal.power.lb:  Lower bound on total power constraint (scalar)
%  signal is either 'excitation' for excitation spectrum, 'input' for input spectrum and 
%  'output' for output spectrum. The excitation and input spectrum coincide
%  in an open-loop set-up. 
%
%  Quality and application constraints are added to the constraints property of
%  prob = oidProblem(...). The ith constraint is specified by
%   prob.constraints{i} = oidConstraint(...);
%  The specification depends on the specific type of constraint.
%  Type help oidApplicationConstraint 
%       help oidQualityConstraint 
%       help oidEllipsoidalQualityConstraint
%  for details on specifying constraints.
%
% Solving the optimization:
%  The optimal input design problem is solved using
%   optH = solve(prob,objective)
%  where optH is a spectral factor of the optimal input spectrum. objective
%  is dependent on which objective function that is considered. If
%  1. objective = [a_r,a_u,a_y],
%  2. objective = [a_r,a_u_r,a_u_e,a_y_r,a_y_e]
%  3. objective = 'A',
%  4. objective = 'E',
%  5. objective = 'D'.
%   optH = solve(prob,objective,opts) can be used to set options for the
%   optimization. opts is an options structure created using SDPSETTINGS.
%
% See also oidSpectrum, oidSpectrumMA, oidConstraint, oidApplicationConstraint, 
% oidQualityConstraint oidEllipsoidalQualityConstraint, sdpsettings

% Author: Mariette Annergren & Christian A. Larsson
% Copyright (c) 2015 Mariette Annergren & Christian A. Larsson

classdef oidProblem
   properties (Dependent = true)        
      model
   end
   properties
      controller;       % Type: struct
      spectrum;         % Type: oidSpectrum
      constraints;      % Type: cell array of oidConstraint
      N;                % Type: experiment length (scalar)
   end
   properties (Dependent = true)        
      opts;            
   end
   properties (Hidden = true)          
      model_;
      opts_;
   end
   methods
      function obj = oidProblem(varargin)
        if nargin > 0
            % Model (argument 1)
            if nargin >= 1
               if isa(varargin{1},'idpoly')
                  obj.model_ = oidModelPoly(varargin{1});
               else
                  error('The model should be specified as an idpoly object.')
               end
            end
            
            % Experiment length (argument 2)
            if nargin >= 2
               if ~isnumeric(varargin{2}) || ~numel(varargin{2})== 1 || ~varargin{2}>0
                  error('Experiment length (samples) should be a positive scalar.');
               end
               obj.N = varargin{2};
            else
               obj.N = 1;
            end
            
            % Spectrum type and number of spectrum parameters (argument 3 and 4)
            if nargin >= 4
               if ~isnumeric(varargin{4}) || ~numel(varargin{4})==1 || ~varargin{4}>0
                  error('Number of spectrum coefficient should be a positive scalar.');
               end
               if strcmpi(varargin{3},'MA') || strcmpi(varargin{3},'FIR')
                  obj.spectrum = oidSpectrumMA(varargin{4},size(obj.model,2));
               elseif strcmpi(varargin{3},'AR')
                  error('AR spectrum not implemented yet. Try MA (FIR) spectrum.');
               elseif strcmpi(varargin{3},'LA')
                  error('LA spectrum not implemented yet. Try MA (FIR) spectrum.');
               elseif strcmpi(varargin{3},'MS')
                  error('Multisine spectrum not implemented yet. Try MA (FIR) spectrum.');
               else
                  error('Spectrum type should be MA (FIR), AR, LA or MS.');
               end
            else
               obj.spectrum = oidSpectrumMA(20,size(obj.model,2));
            end
            
            % Parse and set controller (argument 5 and 6)
            obj.controller =...
               struct('K',tf(zeros(fliplr(size(obj.model)))),'free',false,'openLoop',true); % controller is a struct with fields 'K' and 'free'.
            if nargin >= 5
               if ~isa(varargin{5},'tf')
                  error('Controller should be an LTI transfer function.');
               else
                  obj.controller.K = varargin{5};
                  obj.controller.openLoop = false;
               end
            end
            if nargin >= 6
               if strcmpi(varargin{6},'free')
                  obj.controller.free = true;
                  error('FREE not implemented yet. Try FIXED.');
               elseif strcmpi(varargin{6},'fixed');
                  obj.controller.free = false;
               else
                  error('Controller can be FIXED or FREE');
               end
            end
            
            % Default optimization options (argument 7)
            obj.opts_              = sdpsettings('solver','sdpt3');
            obj.opts_.sdpt3.gaptol = 1e-8;
            obj.opts_.inftol       = 1e-5;
            obj.opts_.steptop      = 1e-4;
            obj.opts_.maxit        = 200;
            obj.opts_.sdpt3.maxit  = 200;
         end

         if exist('yalmip','file') == 0
            warning('MOOSE could not detect yalmip in the Matlab path')
         end
      end 
      function varargout = solve(obj,coeff,varargin)
         con = [];
         % Objective
         if ischar(coeff)
             if strcmp(coeff,'A'); % A-optimality
                 infoMatrix = obj.spectrum.informationMatrix(obj);
                 nTheta     = length(infoMatrix);
                 u          = sdpvar(nTheta,1);
                 for i = 1:nTheta
                     e   = zeros(nTheta,1); e(i) = 1;
                     con = [con; [infoMatrix e; e' u(i)]>=0];
                 end
                 h = sum(u);    
             end 
             if strcmp(coeff,'E');  % E-optimality
                 infoMatrix = obj.spectrum.informationMatrix(obj);
                 nTheta     = length(infoMatrix);
                 t          = sdpvar(1,1);
                 con        = [con; infoMatrix >= t*eye(nTheta)];
                 h          = -t;
             end 
         if strcmp(coeff,'D'); h = -logdet(obj.spectrum.informationMatrix(obj)); end % D-optimality                        
         elseif length(coeff) == 3
             if coeff(1) == 0; Pr     = 0; else Pr = obj.spectrum.excitationPower(obj); end
             if coeff(2) == 0; Pu.tot = 0; else Pu = obj.spectrum.inputPower(obj);      end
             if coeff(3) == 0; Py.tot = 0; else Py = obj.spectrum.outputPower(obj);     end
             h = coeff(1)*Pr + coeff(2)*Pu.tot + coeff(3)*Py.tot;
         elseif length(coeff) == 5
             if coeff(1) == 0;                  Pr   = 0;           else Pr = obj.spectrum.excitationPower(obj); end
             if coeff(2) == 0 || coeff(3) == 0; Pu.r = 0; Pu.e = 0; else Pu = obj.spectrum.inputPower(obj);      end
             if coeff(4) == 0 || coeff(5) == 0; Py.r = 0; Py.e = 0; else Py = obj.spectrum.outputPower(obj);     end
             h = coeff(1)*Pr + coeff(2)*Pu.r + coeff(3)*Pu.e + coeff(4)*Py.r+ coeff(5)*Py.e;
         else
             error('objectiveFunction::Objective function is declared wrong.')
         end
         
         % Spectral constraints
         try
            con = [con; obj.spectrum.getConstraints(obj)];
         catch ME
            throw(ME)
         end
         
         % Constraints
         for k = 1:numel(obj.constraints)
            try
               con = [con; obj.constraints{k}.getConstraints(obj)];
            catch ME
               throw(ME)
            end
         end
         
         % Yalmip options
         if nargin > 2
            opts = varargin{1};
         else
            opts = obj.opts;
         end
         
         % Solution
         try
            d = solvesdp(con,h,opts);
         catch ME
            d = [];
            throw(ME);
         end
         if ~isempty(d)
             if d.problem == 0
                 disp('Input design problem is FEASIBLE')
             elseif d.problem == 1
                 throw(MException('oidProblem:solution','oidProblem:Input design problem is INFEASIBLE'))
             elseif d.problem == 4
                 warningOutput = warning('on','oidProblem:solution');
                 warning('oidProblem:solution','oidProblem:Numerical problems when solving optimization problem');
                 warning(warningOutput)
             else
                 fprintf('oidProblem:Something happened when solving optimization problem, check yalmiperror(%i)\n',d.problem)
             end
         end
         try
             varargout{1} = obj.spectrum.spectralFactor(obj.model.Ts);
         catch ME
             varargout{1} = [];
             warningOutput = warning('on','oidProblem:numericalProblems');
             warning('oidProblem:numericalProblems',...
                 'oidProblem:Numerical problems when calculating spectral factor.\n Error: %s',ME.message);
             warning(warningOutput)
         end

         % Output
         if nargout > 1
            varargout{2} = d;
         end
         if nargout > 2
            varargout{3} = obj.informationMatrix;
         end
         if nargout > 3
            varargout{4} = value(h);
         end
         if nargout > 4
             power = struct('excitation',value(obj.spectrum.excitationPower),'input',value(obj.spectrum.inputPower(obj)),...
                                    'output',value(obj.spectrum.outputPower(obj)));
             varargout{5} = power;
         end
         if nargout > 5
             varargout{6} = obj.spectrum.parameters;
         end
         if nargout > 6
             varargout{7} = con;
         end
      end
      
      % Additional functions
      function iF = informationMatrix(obj)
         try
            iF = double(obj.spectrum.informationMatrix(obj));
         catch
            iF = [];
         end
      end
      
      function [dG,dH] = modelGradient(obj)
         [dG,dH] = obj.model_.gradient();
      end
      
      function model = get.model(obj)
         model = obj.model_.model;
      end
      function opts = get.opts(obj)
         opts = obj.opts_;
      end
      
      function obj = set.constraints(obj,con)
         for k = 1:numel(con)
            if ~isa(con{k},'oidConstraint')
               throw(MException('VerifyConstraint:IncompatibleConstraint',...
                  ['Constraints should be implemented as subclasses '...
                  'of the @oidConstraint superclass']));
            end
         end
         obj.constraints = con;
      end
      
      function obj = set.controller(obj,con)
         if isa(con,'tf')
            obj.controller.K = con;
         elseif isstruct(con) && isfield(con,'K') && isfield(con,'free')
            if isa(con.K,'tf') || isnumeric(con.K)
               obj.controller = con;
            end
         end
      end
      
      function obj = set.model(obj,model)
         if isa(model,'idpoly') && ~model.Ts == 0
            obj.model_ = oidModelPoly(model);
         elseif isa(model,'tf') && ~model.Ts == 0
            obj.model_ = oidModelPoly(idpoly(model));
            warning('Model converted to output error idpoly model');
         else
            throw(MException('VerifyModel:IncompatibleModel',...
               'Only discrete time idpoly models can be used'));
         end
      end
      
      function obj = set.opts(obj,opts)
         obj.opts_ = opts;
      end
      
      function obj = set.spectrum(obj,spectrum)
         if ~isa(spectrum,'oidSpectrum')
            throw(MException('VerifySpectrum:IncompatibleSpectrum',...
               ['Spectra should be implemented as subclasses of the '...
               '@oidSpectrum superclass']));
         end
         obj.spectrum = spectrum;
      end
   end
   
   methods (Static)         
      function R = acf(G,M)
         [A,B,C,D] = ssdata(minreal(G));
         R1 = B*B';
         R2 = D*D';
         R12 = B*D';
         P = dlyap(A,R1);
         if norm(A*P*A'-P+R1) > 1e2
             throw(MException('Informationmatrix:numericalProblems',...
               ['dlyap fails. Probably the eigenvalues are too close to the boundary of the unit disc. Try modifying the oidProblem slightly.']));
         end
         R(:,:,1) = C*P*C'+R2;
         for i = 1:M-1
            R(:,:,i+1) = C*A^(i-1)*(A*P*C'+R12);
         end
      end
      
      function Gvec = vecTF(G)
         [m,n] = size(G);
         Gvec = tf(zeros(m*n,1));
         for i = 1:n 
            Gvec(m*(i-1)+1:i*m) = G(:,i);
         end
      end
   end
end