% Representation of input spectrum with finite dimensional parameterization
%
%  The oidSpectrumFDP class is the base class for all spectra with finite
%  dimensional parameterization used in the MOOSE2 input design toolbox. The 
%  class is not user-interfacing and cannot be instantiated. 
%
%  The oidSpectrumFDP class introduces frequency-wise constraints for the input
%  and output spectrum, and power constraints for input and output signals.
%
% See also oidSpectrumMA


% Author: Mariette Annergren & Christian A. Larsson
% Copyright (c) 2015 Mariette Annergren & Christian A. Larsson
 
classdef (Abstract) oidSpectrumFDP < oidSpectrum
   % Class for Finite Dimensional Spectrum Parameterizations
   properties
      excitation % r(t) if closed-loop, r(t)=u(t) if open-loop
      input      % u(t), u(t)=r(t) is open-loop
      output     % y(t), the same in both open and closed-loop
   end
   properties (Dependent = true)
      M;          % Number of spectrum parameters
   end
   properties (Hidden = true)
      parameters; % Spectrum parameters
   end
   properties (Dependent = true, Hidden = true)
      n;          % Input dimensions
   end
   methods (Abstract)
      [A,B,C,D] = ss(obj);
      [A,B,C,D] = ssBound(R)
   end
   methods
      function obj = oidSpectrumFDP(type,M,n)
         obj            = obj@oidSpectrum(type);
         obj.parameters = sdpvar(n,n,M,'full');
         obj.output     = struct('lb', [],...
                                 'ub', [],...
                                 'power', struct('lb',   [], 'ub',   []));
         obj.excitation = struct('lb', [],...
                                 'ub', [],...
                                 'power', struct('lb',   [], 'ub',   []));
         obj.input      = struct('lb', [],...
                                 'ub', [],...
                                 'power', struct('lb',   [], 'ub',   []));
      end

      function con = getConstraints(obj,problem)
         con = [];
         % KYP for excitation spectrum
         [A,B,C,D] = obj.ss();
         
         %% Excitation spectrum
         % Lower bound excitation spectrum   
         QexLbCompulsory = sdpvar(size(A,2),size(A,2),'symmetric');
         con             = [con; [QexLbCompulsory-A'*QexLbCompulsory*A  C'-A'*QexLbCompulsory*B;
                                C-B'*QexLbCompulsory*A  D+D'-B'*QexLbCompulsory*B] >= 0]; % LMI CONSTRAINT TO ENSURE EXCITATION SPECTRUM IS POSITIVE ENTITY
         if ~isempty(obj.excitation.lb)
             if isa(obj.excitation.lb,'lti')
                 coeffMax = obj.M + 50;
                 for Coeff = obj.M+1:coeffMax
                     RexLb   = oidProblem.acf(tf(obj.excitation.lb),Coeff);
                     if max(max(abs(RexLb(:,:,end)))) < 1e-6
                         break
                     end
                 end
                 RexTrue                   = cat(3,obj.parameters,zeros(size(obj.parameters,1),size(obj.parameters,1),Coeff-obj.M));
                 [AexLb,BexLb,CexLb,DexLb] = obj.ssBound(RexTrue-RexLb);
                 QexLb                     = sdpvar(size(AexLb,2),size(AexLb,2),'symmetric');
                 con                       = [con; [QexLb-AexLb'*QexLb*AexLb  CexLb'-AexLb'*QexLb*BexLb;
                                              CexLb-BexLb'*QexLb*AexLb  DexLb+DexLb'-BexLb'*QexLb*BexLb] >= 0]; 
             elseif isa(obj.excitation.lb,'struct')
                 BexLb = obj.excitation.lb.B;
                 wexLb = obj.excitation.lb.w;
                 for k = 1:length(wexLb);
                     PexLb = C*((exp(1i*wexLb(k)*problem.model.Ts)*...
                         eye(size(A)) - A)\eye(size(A)))*B + D;
                     PexLb = PexLb+PexLb';
                     if strcmp(obj.excitation.lb.type,'lmi')
                         if ~ishermitian(BexLb(:,:,k))
                             if max(max(abs(BexLb(:,:,k)-BexLb(:,:,k)'))) < 1e-8
                                 warning('Lower bound on excitation is numerically non-Hermitian. The bound is made Hermitian.')
                                 BexLb(:,:,k) = 0.5*(BexLb(:,:,k)+BexLb(:,:,k)');
                             else
                                 warning('Lower bound on excitation is non-Hermitian. Linear matrix constraint is handled as an element-wise constraint.')
                             end
                         end
                         con = [con; -BexLb(:,:,k)+PexLb >= 0];
                     else
                         for i = 1:size(PexLb,1)
                             for j = 1:size(PexLb,2)
                                 con = [con; -squeeze(BexLb(i,j,k))+PexLb(i,j) >= 0];
                             end
                         end
                     end
                 end
             elseif size(obj.excitation.lb) == size(D)
                 QexLb = sdpvar(size(A,2),size(A,2),'symmetric');
                 if ~ishermitian(obj.excitation.lb)
                     if max(max(abs(obj.excitation.lb-obj.excitation.lb'))) < 1e-8
                         warning('Lower bound on excitation is numerically non-Hermitian. The bound is made Hermitian.')
                         obj.excitation.lb = 0.5*(obj.excitation.lb+obj.excitation.lb');
                         con               = [con; [QexLb-A'*QexLb*A  C'-A'*QexLb*B;
                                                 C-B'*QexLb*A  D+D'-obj.excitation.lb-B'*QexLb*B] >= 0];
                     else
                         warning('Lower bound on excitation is non-Hermitian. A frequency-sampled element-wise constraint is used instead of the KYP lemma.')
                         wexLb = linspace(1e-3,pi);
                         for k = 1:length(wexLb);
                             PexLb = C*((exp(1i*wexLb(k)*problem.model.Ts)*...
                                 eye(size(A)) - A)\eye(size(A)))*B + D;
                             PexLb = PexLb+PexLb';
                             for i = 1:size(PexLb,1)
                                 for j = 1:size(PexLb,2)
                                     con = [con; -obj.excitation.lb(i,j)+PexLb(i,j) >= 0];
                                 end
                             end
                         end
                     end
                 else
                 con   = [con; [QexLb-A'*QexLb*A  C'-A'*QexLb*B;
                           C-B'*QexLb*A  D+D'-obj.excitation.lb-B'*QexLb*B] >= 0];
                 end
             else
                 error('Lower bound on excitation spectrum is specified wrong.')
             end
         end
         % Upper bound excitation spectrum
         if ~isempty(obj.excitation.ub)
             % Sample LTI constraint
             if isa(obj.excitation.ub,'lti')
                 [BexUb,wexUb] = freqresp(obj.excitation.ub*obj.excitation.ub',...
                     linspace(1e-3,pi)./problem.model.Ts);
                 if max(max(max(abs(imag(BexUb)))))<1e-6
                     BexUb = real(BexUb);
                 end
             elseif isa(obj.excitation.ub,'struct')
                 BexUb = obj.excitation.ub.B;
                 wexUb = obj.excitation.ub.w;
             elseif size(obj.excitation.ub) == size(D)
                 wexUb = linspace(1e-3,pi)./problem.model.Ts;
                 BexUb = repmat(obj.excitation.ub,1,1,length(wexUb));
             end
             for k = 1:length(wexUb);
                 PexUb = C*((exp(1i*wexUb(k)*problem.model.Ts)*...
                     eye(size(A)) - A)\eye(size(A)))*B + D;
                 PexUb = PexUb+PexUb';
                if isa(obj.excitation.ub,'struct') && ~strcmp(obj.excitation.ub.type,'lmi')
                    for i = 1:size(PexUb,1)
                        for j = 1:size(PexUb,2)
                            con = [con; BexUb(i,j,k)-PexUb(i,j) >= 0];
                        end
                    end
                else
                    if ~ishermitian(BexUb(:,:,k))
                        if max(max(abs(BexUb(:,:,k)-BexUb(:,:,k)'))) < 1e-8
                            warning('Upper bound on excitation is numerically non-Hermitian. The bound is made Hermitian.')
                            BexUb(:,:,k) = 0.5*(BexUb(:,:,k)+BexUb(:,:,k)');
                        else
                            warning('Upper bound on excitation is non-Hermitian. Linear matrix constraint is handled as an element-wise constraint.')
                        end
                    end
                    con = [con; BexUb(:,:,k)-PexUb >= 0];
                end
             end
         end
         %% Input spectrum
         % Lower bound input spectrum
         if ~isempty(obj.input.lb)
             if isa(obj.input.lb,'lti')
                 [~,m] = size(problem.model_);
                 [G,H] = problem.model_.tf();
                 Su_r  = minreal(eye(m)/(eye(m) - problem.controller.K*G));
                 Su_e  = minreal(eye(m)/(eye(m) - problem.controller.K*G)*problem.controller.K*H);
                 if ~(problem.controller.openLoop)
                     % Contribution from excitation signal (r)
                     [outU_r, ~]                = size(Su_r);
                     [GammaUPos_r, GammaUNeg_r] = obj.gammaR(problem,Su_r);
                     inLbCoeff_r                = size(GammaUPos_r,3);
                     Rtemp_r                    = [];
                     for m = 0:(obj.M-1) + (inLbCoeff_r-1)
                         switchLag = 1;
                         for m1 = -(inLbCoeff_r-1):(inLbCoeff_r-1)
                             for m2 = -(obj.M-1):(obj.M-1)
                                 if (m1+m2) == m
                                     if switchLag
                                         if isempty(Rtemp_r)
                                             if m1 >= 0
                                                 Rtemp_r     = reshape(GammaUPos_r(:,:,abs(m1)+1)*vec(obj.parameters(:,:,abs(m2)+1)),outU_r,outU_r);
                                                 switchLag = 0;
                                             else
                                                 Rtemp_r     = reshape(GammaUNeg_r(:,:,abs(m1)+1)*vec(obj.parameters(:,:,abs(m2)+1)),outU_r,outU_r);
                                                 switchLag = 0;
                                             end
                                         else
                                             if m1 >= 0
                                                 Rtemp_r     = cat(3,Rtemp_r,reshape(GammaUPos_r(:,:,abs(m1)+1)*vec(obj.parameters(:,:,abs(m2)+1)),outU_r,outU_r));
                                                 switchLag = 0;
                                             else
                                                 Rtemp_r     = cat(3,Rtemp_r,reshape(GammaUNeg_r(:,:,abs(m1)+1)*vec(obj.parameters(:,:,abs(m2)+1)),outU_r,outU_r));
                                                 switchLag = 0;
                                             end
                                         end
                                     else
                                         if m == 0
                                             if m1 >= 0
                                                 Rtemp_r = Rtemp_r + reshape(GammaUPos_r(:,:,abs(m1)+1)*vec(obj.parameters(:,:,abs(m2)+1)),outU_r,outU_r);
                                             else
                                                 Rtemp_r = Rtemp_r + reshape(GammaUNeg_r(:,:,abs(m1)+1)*vec(obj.parameters(:,:,abs(m2)+1)),outU_r,outU_r);
                                             end
                                         else
                                             if m1 >= 0
                                                 Rtemp_r(:,:,m+1) = Rtemp_r(:,:,m+1) + reshape(GammaUPos_r(:,:,abs(m1)+1)*vec(obj.parameters(:,:,abs(m2)+1)),outU_r,outU_r);
                                             else
                                                 Rtemp_r(:,:,m+1) = Rtemp_r(:,:,m+1) + reshape(GammaUNeg_r(:,:,abs(m1)+1)*vec(obj.parameters(:,:,abs(m2)+1)),outU_r,outU_r);
                                             end
                                         end
                                     end
                                 end
                             end
                         end
                     end
                     % Contribution from noise signal (e)
                     [outU_e, ~]      = size(Su_e);
                     [GammaUPos_e, ~] = obj.gammaR(problem,Su_e);
                     inLbCoeff_e      = size(GammaUPos_e,3);
                     Lambda           = problem.model_.NoiseVariance;
                     Rtemp_e          = [];
                     for m = 0:(inLbCoeff_e-1)
                         if isempty(Rtemp_e)
                             Rtemp_e = reshape(GammaUPos_e(:,:,abs(m)+1)*vec(Lambda),outU_e,outU_e);
                         else
                             Rtemp_e = cat(3,Rtemp_e,reshape(GammaUPos_e(:,:,abs(m)+1)*vec(Lambda),outU_e,outU_e));
                         end
                     end
                     if size(Rtemp_r,3) > size(Rtemp_e,3)
                         Rtemp_e = cat(3,Rtemp_e,zeros(size(Rtemp_e,1),size(Rtemp_e,2),size(Rtemp_r,3)-size(Rtemp_e,3)));
                     end
                     if size(Rtemp_e,3) > size(Rtemp_r,3)
                         Rtemp_r = cat(3,Rtemp_r,zeros(size(Rtemp_r,1),size(Rtemp_r,2),size(Rtemp_e,3)-size(Rtemp_r,3)));
                     end
                     RinLb                     = oidProblem.acf(tf(obj.input.lb),size(Rtemp_r,3));
                     [AinLb,BinLb,CinLb,DinLb] = obj.ssBound(Rtemp_r+Rtemp_e-RinLb);
                     QinLb                     = sdpvar(size(AinLb,2),size(AinLb,2),'symmetric');
                     con                       = [con; [QinLb-AinLb'*QinLb*AinLb  CinLb'-AinLb'*QinLb*BinLb;
                         CinLb-BinLb'*QinLb*AinLb  DinLb+DinLb'-BinLb'*QinLb*BinLb] >= 0];
                 else % open-loop
                     coeffMax = obj.M + 50;
                     for Coeff = obj.M+1:coeffMax
                         RinLb   = oidProblem.acf(tf(obj.input.lb),Coeff);
                         if max(max(abs(RinLb(:,:,end)))) < 1e-6
                             break
                         end
                     end
                     RinTrue                   = cat(3,obj.parameters,zeros(size(obj.parameters,1),size(obj.parameters,1),Coeff-obj.M));
                     [AinLb,BinLb,CinLb,DinLb] = obj.ssBound(RinTrue-RinLb);
                     QinLb                     = sdpvar(size(AinLb,2),size(AinLb,2),'symmetric');
                     con                       = [con; [QinLb-AinLb'*QinLb*AinLb  CinLb'-AinLb'*QinLb*BinLb;
                         CinLb-BinLb'*QinLb*AinLb  DinLb+DinLb'-BinLb'*QinLb*BinLb] >= 0];
                 end
             elseif isa(obj.input.lb,'struct')
                 [~,m] = size(problem.model_);
                 [G,H] = problem.model_.tf();
                 Su_r  = minreal(eye(m)/(eye(m) - problem.controller.K*G));
                 Su_e  = minreal(eye(m)/(eye(m) - problem.controller.K*G)*problem.controller.K*H);
                 if ~(problem.controller.openLoop)
                     % Contribution from excitation signal (r)
                     [outU_r, ~]                = size(Su_r);
                     [GammaUPos_r, GammaUNeg_r] = obj.gammaR(problem,Su_r);
                     inLbCoeff_r                = size(GammaUPos_r,3);
                     Rtemp_r                    = [];
                     for m = 0:(obj.M-1) + (inLbCoeff_r-1)
                         switchLag = 1;
                         for m1 = -(inLbCoeff_r-1):(inLbCoeff_r-1)
                             for m2 = -(obj.M-1):(obj.M-1)
                                 if (m1+m2) == m
                                     if switchLag
                                         if isempty(Rtemp_r)
                                             if m1 >= 0
                                                 Rtemp_r     = reshape(GammaUPos_r(:,:,abs(m1)+1)*vec(obj.parameters(:,:,abs(m2)+1)),outU_r,outU_r);
                                                 switchLag = 0;
                                             else
                                                 Rtemp_r     = reshape(GammaUNeg_r(:,:,abs(m1)+1)*vec(obj.parameters(:,:,abs(m2)+1)),outU_r,outU_r);
                                                 switchLag = 0;
                                             end
                                         else
                                             if m1 >= 0
                                                 Rtemp_r     = cat(3,Rtemp_r,reshape(GammaUPos_r(:,:,abs(m1)+1)*vec(obj.parameters(:,:,abs(m2)+1)),outU_r,outU_r));
                                                 switchLag = 0;
                                             else
                                                 Rtemp_r     = cat(3,Rtemp_r,reshape(GammaUNeg_r(:,:,abs(m1)+1)*vec(obj.parameters(:,:,abs(m2)+1)),outU_r,outU_r));
                                                 switchLag = 0;
                                             end
                                         end
                                     else
                                         if m == 0
                                             if m1 >= 0
                                                 Rtemp_r = Rtemp_r + reshape(GammaUPos_r(:,:,abs(m1)+1)*vec(obj.parameters(:,:,abs(m2)+1)),outU_r,outU_r);
                                             else
                                                 Rtemp_r = Rtemp_r + reshape(GammaUNeg_r(:,:,abs(m1)+1)*vec(obj.parameters(:,:,abs(m2)+1)),outU_r,outU_r);
                                             end
                                         else
                                             if m1 >= 0
                                                 Rtemp_r(:,:,m+1) = Rtemp_r(:,:,m+1) + reshape(GammaUPos_r(:,:,abs(m1)+1)*vec(obj.parameters(:,:,abs(m2)+1)),outU_r,outU_r);
                                             else
                                                 Rtemp_r(:,:,m+1) = Rtemp_r(:,:,m+1) + reshape(GammaUNeg_r(:,:,abs(m1)+1)*vec(obj.parameters(:,:,abs(m2)+1)),outU_r,outU_r);
                                             end
                                         end
                                     end
                                 end
                             end
                         end
                     end
                     % Contribution from noise signal (e)
                     [outU_e, ~]      = size(Su_e);
                     [GammaUPos_e, ~] = obj.gammaR(problem,Su_e);
                     inLbCoeff_e      = size(GammaUPos_e,3);
                     Lambda           = problem.model_.NoiseVariance;
                     Rtemp_e          = [];
                     for m = 0:(inLbCoeff_e-1)
                         if isempty(Rtemp_e)
                             Rtemp_e = reshape(GammaUPos_e(:,:,abs(m)+1)*vec(Lambda),outU_e,outU_e);
                         else
                             Rtemp_e = cat(3,Rtemp_e,reshape(GammaUPos_e(:,:,abs(m)+1)*vec(Lambda),outU_e,outU_e));
                         end
                     end
                     if size(Rtemp_r,3) > size(Rtemp_e,3)
                         Rtemp_e = cat(3,Rtemp_e,zeros(size(Rtemp_e,1),size(Rtemp_e,2),size(Rtemp_r,3)-size(Rtemp_e,3)));
                     end
                     if size(Rtemp_e,3) > size(Rtemp_r,3)
                         Rtemp_r = cat(3,Rtemp_r,zeros(size(Rtemp_r,1),size(Rtemp_r,2),size(Rtemp_e,3)-size(Rtemp_r,3)));
                     end
                     [AinLb,BBinLb,CinLb,DinLb] = obj.ssBound(Rtemp_r+Rtemp_e);
                     BinLb = obj.input.lb.B;
                     winLb = obj.input.lb.w;
                     for k = 1:length(winLb);
                         PinLb = CinLb*((exp(1i*winLb(k)*problem.model.Ts)*...
                             eye(size(AinLb)) - AinLb)\eye(size(AinLb)))*BBinLb + DinLb;
                         PinLb = PinLb+PinLb';
                         if strcmp(obj.input.lb.type,'lmi')
                             if ~ishermitian(BinLb(:,:,k))
                                 if max(max(abs(BinLb(:,:,k)-BinLb(:,:,k)'))) < 1e-8
                                     warning('Lower bound on input is numerically non-Hermitian. The bound is made Hermitian.')
                                     BinLb(:,:,k) = 0.5*(BinLb(:,:,k)+BinLb(:,:,k)');
                                 else
                                     warning('Lower bound on input is non-Hermitian. Linear matrix constraint is handled as an element-wise constraint.')
                                 end
                             end
                             con = [con; -BinLb(:,:,k)+PinLb >= 0];
                         else
                             for i = 1:size(PinLb,1)
                                 for j = 1:size(PinLb,2)
                                     con = [con; -BinLb(i,j,k)+PinLb(i,j) >= 0];
                                 end
                             end
                         end
                     end
                 else % open-loop
                     BinLb = obj.input.lb.B;
                     winLb = obj.input.lb.w;
                     for k = 1:length(winLb);
                         PinLb = C*((exp(1i*winLb(k)*problem.model.Ts)*...
                             eye(size(A)) - A)\eye(size(A)))*B + D;
                         PinLb = PinLb+PinLb';
                         if strcmp(obj.input.lb.type,'lmi')
                             if ~issymmetric(BinLb(:,:,k))
                                 if max(max(abs(BinLb(:,:,k)-BinLb(:,:,k)'))) < 1e-8
                                     warning('Lower bound on input is numerically non-Hermitian. The bound is made Hermitian.')
                                     BinLb(:,:,k) = 0.5*(BinLb(:,:,k)+BinLb(:,:,k)');
                                 else
                                     warning('Lower bound on input is non-Hermitian. Linear matrix constraint is handled as an element-wise constraint.')
                                 end
                             end
                             con = [con; -BinLb(:,:,k)+PinLb >= 0];
                         else
                             for i = 1:size(PinLb,1)
                                 for j = 1:size(PinLb,2)
                                     con = [con; -BinLb(i,j,k)+PinLb(i,j) >= 0];
                                 end
                             end
                         end
                     end
                 end
             elseif size(obj.input.lb) == size(D)
                 [~,m] = size(problem.model_);
                 [G,H] = problem.model_.tf();
                 Su_r  = minreal(eye(m)/(eye(m) - problem.controller.K*G));
                 Su_e  = minreal(eye(m)/(eye(m) - problem.controller.K*G)*problem.controller.K*H);
                 if ~(problem.controller.openLoop)
                     % Contribution from excitation signal (r)
                     [outU_r, ~]                = size(Su_r);
                     [GammaUPos_r, GammaUNeg_r] = obj.gammaR(problem,Su_r);
                     inLbCoeff_r                = size(GammaUPos_r,3);
                     Rtemp_r                    = [];
                     for m = 0:(obj.M-1) + (inLbCoeff_r-1)
                         switchLag = 1;
                         for m1 = -(inLbCoeff_r-1):(inLbCoeff_r-1)
                             for m2 = -(obj.M-1):(obj.M-1)
                                 if (m1+m2) == m
                                     if switchLag
                                         if isempty(Rtemp_r)
                                             if m1 >= 0
                                                 Rtemp_r     = reshape(GammaUPos_r(:,:,abs(m1)+1)*vec(obj.parameters(:,:,abs(m2)+1)),outU_r,outU_r);
                                                 switchLag = 0;
                                             else
                                                 Rtemp_r     = reshape(GammaUNeg_r(:,:,abs(m1)+1)*vec(obj.parameters(:,:,abs(m2)+1)),outU_r,outU_r);
                                                 switchLag = 0;
                                             end
                                         else
                                             if m1 >= 0
                                                 Rtemp_r     = cat(3,Rtemp_r,reshape(GammaUPos_r(:,:,abs(m1)+1)*vec(obj.parameters(:,:,abs(m2)+1)),outU_r,outU_r));
                                                 switchLag = 0;
                                             else
                                                 Rtemp_r     = cat(3,Rtemp_r,reshape(GammaUNeg_r(:,:,abs(m1)+1)*vec(obj.parameters(:,:,abs(m2)+1)),outU_r,outU_r));
                                                 switchLag = 0;
                                             end
                                         end
                                     else
                                         if m == 0
                                             if m1 >= 0
                                                 Rtemp_r = Rtemp_r + reshape(GammaUPos_r(:,:,abs(m1)+1)*vec(obj.parameters(:,:,abs(m2)+1)),outU_r,outU_r);
                                             else
                                                 Rtemp_r = Rtemp_r + reshape(GammaUNeg_r(:,:,abs(m1)+1)*vec(obj.parameters(:,:,abs(m2)+1)),outU_r,outU_r);
                                             end
                                         else
                                             if m1 >= 0
                                                 Rtemp_r(:,:,m+1) = Rtemp_r(:,:,m+1) + reshape(GammaUPos_r(:,:,abs(m1)+1)*vec(obj.parameters(:,:,abs(m2)+1)),outU_r,outU_r);
                                             else
                                                 Rtemp_r(:,:,m+1) = Rtemp_r(:,:,m+1) + reshape(GammaUNeg_r(:,:,abs(m1)+1)*vec(obj.parameters(:,:,abs(m2)+1)),outU_r,outU_r);
                                             end
                                         end
                                     end
                                 end
                             end
                         end
                     end
                     % Contribution from noise signal (e)
                     [outU_e, ~]      = size(Su_e);
                     [GammaUPos_e, ~] = obj.gammaR(problem,Su_e);
                     inLbCoeff_e      = size(GammaUPos_e,3);
                     Lambda           = problem.model_.NoiseVariance;
                     Rtemp_e          = [];
                     for m = 0:(inLbCoeff_e-1)
                         if isempty(Rtemp_e)
                             Rtemp_e = reshape(GammaUPos_e(:,:,abs(m)+1)*vec(Lambda),outU_e,outU_e);
                         else
                             Rtemp_e = cat(3,Rtemp_e,reshape(GammaUPos_e(:,:,abs(m)+1)*vec(Lambda),outU_e,outU_e));
                         end
                     end
                     if size(Rtemp_r,3) > size(Rtemp_e,3)
                         Rtemp_e = cat(3,Rtemp_e,zeros(size(Rtemp_e,1),size(Rtemp_e,2),size(Rtemp_r,3)-size(Rtemp_e,3)));
                     end
                     if size(Rtemp_e,3) > size(Rtemp_r,3)
                         Rtemp_r = cat(3,Rtemp_r,zeros(size(Rtemp_r,1),size(Rtemp_r,2),size(Rtemp_e,3)-size(Rtemp_r,3)));
                     end
                     [AinLb,BinLb,CinLb,DinLb] = obj.ssBound(Rtemp_r+Rtemp_e);
                     if ~ishermitian(obj.input.lb)
                         if max(max(abs(obj.input.lb-obj.input.lb'))) < 1e-8
                             warning('Lower bound on input is numerically non-Hermitian. The bound is made Hermitian.')
                             obj.input.lb = 0.5*(obj.input.lb+obj.input.lb');
                             QinLb = sdpvar(size(AinLb,2),size(AinLb,2),'symmetric');
                             con   = [con; [QinLb-AinLb'*QinLb*AinLb  CinLb'-AinLb'*QinLb*BinLb;
                                 CinLb-BinLb'*QinLb*AinLb  DinLb+DinLb'-obj.input.lb-BinLb'*QinLb*BinLb] >= 0];
                         else
                             warning('Lower bound on input is non-Hermitian. A frequency-sampled element-wise constraint is used instead of the KYP lemma.')
                             winLb = linspace(1e-3,pi);
                             for k = 1:length(winLb);
                                 PinLb = CinLb*((exp(1i*winLb(k)*problem.model.Ts)*...
                                     eye(size(AinLb)) - AinLb)\eye(size(AinLb)))*BinLb + DinLb;
                                 PinLb = PinLb+PinLb';
                                 for i = 1:size(PinLb,1)
                                     for j = 1:size(PinLb,2)
                                         con = [con; -obj.input.lb(i,j)+PinLb(i,j) >= 0];
                                     end
                                 end
                             end
                         end
                     else
                         QinLb = sdpvar(size(AinLb,2),size(AinLb,2),'symmetric');
                         con   = [con; [QinLb-AinLb'*QinLb*AinLb  CinLb'-AinLb'*QinLb*BinLb;
                             CinLb-BinLb'*QinLb*AinLb  DinLb+DinLb'-obj.input.lb-BinLb'*QinLb*BinLb] >= 0];
                     end
                 else % open-loop
                     if ~issymmetric(obj.input.lb)
                         if max(max(abs(obj.input.lb-obj.input.lb'))) < 1e-8
                             warning('Lower bound on input is numerically non-Hermitian. The bound is made Hermitian.')
                             obj.input.lb = 0.5*(obj.input.lb+obj.input.lb');
                             QinLb = sdpvar(size(A,2),size(A,2),'symmetric');
                             con   = [con; [QinLb-A'*QinLb*A  C'-A'*QinLb*B;
                                 C-B'*QinLb*A  D+D'-obj.input.lb-B'*QinLb*B] >= 0];
                         else
                             warning('Lower bound on input is non-Hermitian. A frequency-sampled element-wise constraint is used instead of the KYP lemma.')
                             winLb = linspace(1e-3,pi);
                             for k = 1:length(winLb);
                                 PinLb = C*((exp(1i*winLb(k)*problem.model.Ts)*...
                                     eye(size(A)) - A)\eye(size(A)))*B + D;
                                 PinLb = PinLb+PinLb';
                                 for i = 1:size(PinLb,1)
                                     for j = 1:size(PinLb,2)
                                         con = [con; -obj.input.lb(i,j)+PinLb(i,j) >= 0];
                                     end
                                 end
                             end
                         end
                     else
                         QinLb = sdpvar(size(A,2),size(A,2),'symmetric');
                         con   = [con; [QinLb-A'*QinLb*A  C'-A'*QinLb*B;
                             C-B'*QinLb*A  D+D'-obj.input.lb-B'*QinLb*B] >= 0];
                     end
                 end
             else
                 error('Lower bound on input spectrum is specified wrong.')
             end
         end
         % Upper bound input spectrum
         if ~isempty(obj.input.ub)
             % Sample LTI constraint
             if isa(obj.input.ub,'lti')
                 [BinUb,winUb] = freqresp(obj.input.ub*obj.input.ub',linspace(1e-3,pi)./problem.model.Ts);
                 if max(max(max(abs(imag(BinUb)))))<1e-6
                     BinUb = real(BinUb);
                 end
             elseif isa(obj.input.ub,'struct')
                 BinUb = obj.input.ub.B;
                 winUb = obj.input.ub.w;
             elseif size(obj.input.ub) == size(D)
                 winUb = linspace(1e-3,pi)./problem.model.Ts;
                 BinUb = repmat(obj.input.ub,1,1,length(winUb));
             end
             [~,m] = size(problem.model_);
             [G,H] = problem.model_.tf();
             Su_r  = minreal(eye(m)/(eye(m) - problem.controller.K*G));
             Su_e  = minreal(eye(m)/(eye(m) - problem.controller.K*G)*problem.controller.K*H);
             if ~(problem.controller.openLoop)
                     % Contribution from excitation signal (r)
                     [outU_r, ~]                = size(Su_r);
                     [GammaUPos_r, GammaUNeg_r] = obj.gammaR(problem,Su_r);
                     inUbCoeff_r                = size(GammaUPos_r,3);
                     Rtemp_r                    = [];
                     for m = 0:(obj.M-1) + (inUbCoeff_r-1)
                         switchLag = 1;
                         for m1 = -(inUbCoeff_r-1):(inUbCoeff_r-1)
                             for m2 = -(obj.M-1):(obj.M-1)
                                 if (m1+m2) == m
                                     if switchLag
                                         if isempty(Rtemp_r)
                                             if m1 >= 0
                                                 Rtemp_r     = reshape(GammaUPos_r(:,:,abs(m1)+1)*vec(obj.parameters(:,:,abs(m2)+1)),outU_r,outU_r);
                                                 switchLag = 0;
                                             else
                                                 Rtemp_r     = reshape(GammaUNeg_r(:,:,abs(m1)+1)*vec(obj.parameters(:,:,abs(m2)+1)),outU_r,outU_r);
                                                 switchLag = 0;
                                             end
                                         else
                                             if m1 >= 0
                                                 Rtemp_r     = cat(3,Rtemp_r,reshape(GammaUPos_r(:,:,abs(m1)+1)*vec(obj.parameters(:,:,abs(m2)+1)),outU_r,outU_r));
                                                 switchLag = 0;
                                             else
                                                 Rtemp_r     = cat(3,Rtemp_r,reshape(GammaUNeg_r(:,:,abs(m1)+1)*vec(obj.parameters(:,:,abs(m2)+1)),outU_r,outU_r));
                                                 switchLag = 0;
                                             end
                                         end
                                     else
                                         if m == 0
                                             if m1 >= 0
                                                 Rtemp_r = Rtemp_r + reshape(GammaUPos_r(:,:,abs(m1)+1)*vec(obj.parameters(:,:,abs(m2)+1)),outU_r,outU_r);
                                             else
                                                 Rtemp_r = Rtemp_r + reshape(GammaUNeg_r(:,:,abs(m1)+1)*vec(obj.parameters(:,:,abs(m2)+1)),outU_r,outU_r);
                                             end
                                         else
                                             if m1 >= 0
                                                 Rtemp_r(:,:,m+1) = Rtemp_r(:,:,m+1) + reshape(GammaUPos_r(:,:,abs(m1)+1)*vec(obj.parameters(:,:,abs(m2)+1)),outU_r,outU_r);
                                             else
                                                 Rtemp_r(:,:,m+1) = Rtemp_r(:,:,m+1) + reshape(GammaUNeg_r(:,:,abs(m1)+1)*vec(obj.parameters(:,:,abs(m2)+1)),outU_r,outU_r);
                                             end
                                         end
                                     end
                                 end
                             end
                         end
                     end
                     % Contribution from noise signal (e)
                     [outU_e, ~]      = size(Su_e);
                     [GammaUPos_e, ~] = obj.gammaR(problem,Su_e);
                     inUbCoeff_e      = size(GammaUPos_e,3);
                     Lambda           = problem.model_.NoiseVariance;
                     Rtemp_e          = [];
                     for m = 0:(inUbCoeff_e-1)
                         if isempty(Rtemp_e)
                             Rtemp_e = reshape(GammaUPos_e(:,:,abs(m)+1)*vec(Lambda),outU_e,outU_e);
                         else
                             Rtemp_e = cat(3,Rtemp_e,reshape(GammaUPos_e(:,:,abs(m)+1)*vec(Lambda),outU_e,outU_e));
                         end
                     end
                     if size(Rtemp_r,3) > size(Rtemp_e,3)
                         Rtemp_e = cat(3,Rtemp_e,zeros(size(Rtemp_e,1),size(Rtemp_e,2),size(Rtemp_r,3)-size(Rtemp_e,3)));
                     end
                     if size(Rtemp_e,3) > size(Rtemp_r,3)
                         Rtemp_r = cat(3,Rtemp_r,zeros(size(Rtemp_r,1),size(Rtemp_r,2),size(Rtemp_e,3)-size(Rtemp_r,3)));
                     end
                 [AinUb,BBinUb,CinUb,DinUb] = obj.ssBound(Rtemp_r+Rtemp_e);
                 for k = 1:length(winUb);
                     PinUb = CinUb*((exp(1i*winUb(k)*problem.model.Ts)*...
                         eye(size(AinUb)) - AinUb)\eye(size(AinUb)))*BBinUb + DinUb;
                     PinUb = PinUb+PinUb';
                     if isa(obj.input.ub,'struct') && ~strcmp(obj.input.ub.type,'lmi')
                         for i = 1:size(PinUb,1)
                             for j = 1:size(PinUb,2)
                                 con = [con; BinUb(i,j,k)-PinUb(i,j) >= 0];
                             end
                         end
                     else
                         if ~issymmetric(BinUb(:,:,k))
                             if max(max(abs(BinUb(:,:,k)-BinUb(:,:,k)'))) < 1e-8
                                 warning('Upper bound on input is numerically non-Hermitian. The bound is made Hermitian.')
                                 BinUb(:,:,k) = 0.5*(BinUb(:,:,k)+BinUb(:,:,k)');
                             else
                                 warning('Upper bound on input is non-Hermitian. Linear matrix constraint is handled as an element-wise constraint.')
                             end
                         end
                         con = [con; BinUb(:,:,k)-PinUb >= 0];
                     end
                 end
             else % open-loop
                 for k = 1:length(winUb);
                     PinUb = C*((exp(1i*winUb(k)*problem.model.Ts)*...
                         eye(size(A)) - A)\eye(size(A)))*B + D;
                     PinUb = PinUb+PinUb';
                     if isa(obj.input.ub,'struct') && ~strcmp(obj.input.ub.type,'lmi')
                         for i = 1:size(PinUb,1)
                             for j = 1:size(PinUb,2)
                                 con = [con; BinUb(i,j,k)-PinUb(i,j) >= 0];
                             end
                         end
                     else
                         if ~ishermitian(BinUb(:,:,k))
                             if max(max(abs(BinUb(:,:,k)-BinUb(:,:,k)'))) < 1e-8
                                 warning('Upper bound on input is numerically non-Hermitian. The bound is made Hermitian.')
                                 BinUb(:,:,k) = 0.5*(BinUb(:,:,k)+BinUb(:,:,k)');
                             else
                                 warning('Upper bound on input is non-Hermitian. Linear matrix constraint is handled as an element-wise constraint.')
                             end
                         end
                         con = [con; BinUb(:,:,k)-PinUb >= 0];
                     end
                 end
             end
         end
         %% Output spectrum
         % Lower bound output spectrum
         if ~isempty(obj.output.lb)
             if isa(obj.output.lb,'lti')
                 [~,m] = size(problem.model_);
                 [G,H] = problem.model_.tf();
                 Sy_r  = minreal(eye(m)/(eye(m) - G*problem.controller.K)*G);
                 Sy_e  = minreal(eye(m)/(eye(m) - G*problem.controller.K)*H);
                 % Contribution from excitation signal (r)
                 [outY_r, ~]                = size(Sy_r);
                 [GammaYPos_r, GammaYNeg_r] = obj.gammaR(problem,Sy_r);
                 ouLbCoeff_r                = size(GammaYPos_r,3);
                 Rtemp_r                    = [];
                 for m = 0:(obj.M-1) + (ouLbCoeff_r-1)
                     switchLag = 1;
                     for m1 = -(ouLbCoeff_r-1):(ouLbCoeff_r-1)
                         for m2 = -(obj.M-1):(obj.M-1)
                             if (m1+m2) == m
                                 if switchLag
                                     if isempty(Rtemp_r)
                                         if m1 >= 0
                                             Rtemp_r     = reshape(GammaYPos_r(:,:,abs(m1)+1)*vec(obj.parameters(:,:,abs(m2)+1)),outY_r,outY_r);
                                             switchLag = 0;
                                         else
                                             Rtemp_r     = reshape(GammaYNeg_r(:,:,abs(m1)+1)*vec(obj.parameters(:,:,abs(m2)+1)),outY_r,outY_r);
                                             switchLag = 0;
                                         end
                                     else
                                         if m1 >= 0
                                             Rtemp_r     = cat(3,Rtemp_r,reshape(GammaYPos_r(:,:,abs(m1)+1)*vec(obj.parameters(:,:,abs(m2)+1)),outY_r,outY_r));
                                             switchLag = 0;
                                         else
                                             Rtemp_r     = cat(3,Rtemp_r,reshape(GammaYNeg_r(:,:,abs(m1)+1)*vec(obj.parameters(:,:,abs(m2)+1)),outY_r,outY_r));
                                             switchLag = 0;
                                         end
                                     end
                                 else
                                     if m == 0
                                         if m1 >= 0
                                             Rtemp_r = Rtemp_r + reshape(GammaYPos_r(:,:,abs(m1)+1)*vec(obj.parameters(:,:,abs(m2)+1)),outY_r,outY_r);
                                         else
                                             Rtemp_r = Rtemp_r + reshape(GammaYNeg_r(:,:,abs(m1)+1)*vec(obj.parameters(:,:,abs(m2)+1)),outY_r,outY_r);
                                         end
                                     else
                                         if m1 >= 0
                                             Rtemp_r(:,:,m+1) = Rtemp_r(:,:,m+1) + reshape(GammaYPos_r(:,:,abs(m1)+1)*vec(obj.parameters(:,:,abs(m2)+1)),outY_r,outY_r);
                                         else
                                             Rtemp_r(:,:,m+1) = Rtemp_r(:,:,m+1) + reshape(GammaYNeg_r(:,:,abs(m1)+1)*vec(obj.parameters(:,:,abs(m2)+1)),outY_r,outY_r);
                                         end
                                     end
                                 end
                             end
                         end
                     end
                 end
                 % Contribution from noise signal (e)
                 [outY_e, ~]      = size(Sy_e);
                 [GammaYPos_e, ~] = obj.gammaR(problem,Sy_e);
                 ouLbCoeff_e      = size(GammaYPos_e,3);
                 Lambda           = problem.model_.NoiseVariance;
                 Rtemp_e          = [];
                 for m = 0:(ouLbCoeff_e-1)
                     if isempty(Rtemp_e)
                         Rtemp_e = reshape(GammaYPos_e(:,:,abs(m)+1)*vec(Lambda),outY_e,outY_e);
                     else
                         Rtemp_e = cat(3,Rtemp_e,reshape(GammaYPos_e(:,:,abs(m)+1)*vec(Lambda),outY_e,outY_e));
                     end
                 end
                 if size(Rtemp_r,3) > size(Rtemp_e,3)
                     Rtemp_e = cat(3,Rtemp_e,zeros(size(Rtemp_e,1),size(Rtemp_e,2),size(Rtemp_r,3)-size(Rtemp_e,3)));
                 end
                 if size(Rtemp_e,3) > size(Rtemp_r,3)
                     Rtemp_r = cat(3,Rtemp_r,zeros(size(Rtemp_r,1),size(Rtemp_r,2),size(Rtemp_e,3)-size(Rtemp_r,3)));
                 end
                 RouLb                     = oidProblem.acf(tf(obj.output.lb),size(Rtemp_r,3));
                 [AouLb,BouLb,CouLb,DouLb] = obj.ssBound(Rtemp_r+Rtemp_e-RouLb);
                 QouLb                     = sdpvar(size(AouLb,2),size(AouLb,2),'symmetric');
                 con                       = [con; [QouLb-AouLb'*QouLb*AouLb  CouLb'-AouLb'*QouLb*BouLb;
                     CouLb-BouLb'*QouLb*AouLb  DouLb+DouLb'-BouLb'*QouLb*BouLb] >= 0];
             elseif isa(obj.output.lb,'struct')
                 [~,m] = size(problem.model_);
                 [G,H] = problem.model_.tf();
                 Sy_r  = minreal(eye(m)/(eye(m) - G*problem.controller.K)*G);
                 Sy_e  = minreal(eye(m)/(eye(m) - G*problem.controller.K)*H);
                 % Contribution from excitation signal (r)
                 [outY_r, ~]                = size(Sy_r);
                 [GammaYPos_r, GammaYNeg_r] = obj.gammaR(problem,Sy_r);
                 ouLbCoeff_r                = size(GammaYPos_r,3);
                 Rtemp_r                    = [];
                 for m = 0:(obj.M-1) + (ouLbCoeff_r-1)
                     switchLag = 1;
                     for m1 = -(ouLbCoeff_r-1):(ouLbCoeff_r-1)
                         for m2 = -(obj.M-1):(obj.M-1)
                             if (m1+m2) == m
                                 if switchLag
                                     if isempty(Rtemp_r)
                                         if m1 >= 0
                                             Rtemp_r     = reshape(GammaYPos_r(:,:,abs(m1)+1)*vec(obj.parameters(:,:,abs(m2)+1)),outY_r,outY_r);
                                             switchLag = 0;
                                         else
                                             Rtemp_r     = reshape(GammaYNeg_r(:,:,abs(m1)+1)*vec(obj.parameters(:,:,abs(m2)+1)),outY_r,outY_r);
                                             switchLag = 0;
                                         end
                                     else
                                         if m1 >= 0
                                             Rtemp_r     = cat(3,Rtemp_r,reshape(GammaYPos_r(:,:,abs(m1)+1)*vec(obj.parameters(:,:,abs(m2)+1)),outY_r,outY_r));
                                             switchLag = 0;
                                         else
                                             Rtemp_r     = cat(3,Rtemp_r,reshape(GammaYNeg_r(:,:,abs(m1)+1)*vec(obj.parameters(:,:,abs(m2)+1)),outY_r,outY_r));
                                             switchLag = 0;
                                         end
                                     end
                                 else
                                     if m == 0
                                         if m1 >= 0
                                             Rtemp_r = Rtemp_r + reshape(GammaYPos_r(:,:,abs(m1)+1)*vec(obj.parameters(:,:,abs(m2)+1)),outY_r,outY_r);
                                         else
                                             Rtemp_r = Rtemp_r + reshape(GammaYNeg_r(:,:,abs(m1)+1)*vec(obj.parameters(:,:,abs(m2)+1)),outY_r,outY_r);
                                         end
                                     else
                                         if m1 >= 0
                                             Rtemp_r(:,:,m+1) = Rtemp_r(:,:,m+1) + reshape(GammaYPos_r(:,:,abs(m1)+1)*vec(obj.parameters(:,:,abs(m2)+1)),outY_r,outY_r);
                                         else
                                             Rtemp_r(:,:,m+1) = Rtemp_r(:,:,m+1) + reshape(GammaYNeg_r(:,:,abs(m1)+1)*vec(obj.parameters(:,:,abs(m2)+1)),outY_r,outY_r);
                                         end
                                     end
                                 end
                             end
                         end
                     end
                 end
                 % Contribution from noise signal (e)
                 [outY_e, ~]      = size(Sy_e);
                 [GammaYPos_e, ~] = obj.gammaR(problem,Sy_e);
                 ouLbCoeff_e      = size(GammaYPos_e,3);
                 Lambda           = problem.model_.NoiseVariance;
                 Rtemp_e          = [];
                 for m = 0:(ouLbCoeff_e-1)
                     if isempty(Rtemp_e)
                         Rtemp_e = reshape(GammaYPos_e(:,:,abs(m)+1)*vec(Lambda),outY_e,outY_e);
                     else
                         Rtemp_e = cat(3,Rtemp_e,reshape(GammaYPos_e(:,:,abs(m)+1)*vec(Lambda),outY_e,outY_e));
                     end
                 end
                 if size(Rtemp_r,3) > size(Rtemp_e,3)
                     Rtemp_e = cat(3,Rtemp_e,zeros(size(Rtemp_e,1),size(Rtemp_e,2),size(Rtemp_r,3)-size(Rtemp_e,3)));
                 end
                 if size(Rtemp_e,3) > size(Rtemp_r,3)
                     Rtemp_r = cat(3,Rtemp_r,zeros(size(Rtemp_r,1),size(Rtemp_r,2),size(Rtemp_e,3)-size(Rtemp_r,3)));
                 end
                 [AouLb,BBouLb,CouLb,DouLb] = obj.ssBound(Rtemp_r+Rtemp_e);
                 BouLb = obj.output.lb.B;
                 wouLb = obj.output.lb.w;
                 for k = 1:length(wouLb);
                     PouLb = CouLb*((exp(1i*wouLb(k)*problem.model.Ts)*...
                         eye(size(AouLb)) - AouLb)\eye(size(AouLb)))*BBouLb + DouLb;
                     PouLb = PouLb+PouLb';
                     if strcmp(obj.output.lb.type,'lmi')
                         if ~issymmetric(BouLb(:,:,k))
                             if max(max(abs(BouLb(:,:,k)-BouLb(:,:,k)'))) < 1e-8
                                 warning('Lower bound on output is numerically non-Hermitian. The bound is made Hermitian.')
                                 BouLb(:,:,k) = 0.5*(BouLb(:,:,k)+BouLb(:,:,k)');
                             else
                                 warning('Lower bound on output is non-Hermitian. Linear matrix constraint is handled as an element-wise constraint.')
                             end
                         end
                         con = [con; -BouLb(:,:,k)+PouLb >= 0];
                     else
                         for i = 1:size(PouLb,1)
                             for j = 1:size(PouLb,2)
                                 con = [con; -BouLb(i,j,k)+PouLb(i,j) >= 0];
                             end
                         end
                     end
                 end
             elseif size(obj.output.lb,1) == size(problem.model_,1) && size(obj.output.lb,2) == size(problem.model_,1)
                 BouLb = obj.output.lb;
                 [~,m] = size(problem.model_);
                 [G,H] = problem.model_.tf();
                 Sy_r  = minreal(eye(m)/(eye(m) - G*problem.controller.K)*G);
                 Sy_e  = minreal(eye(m)/(eye(m) - G*problem.controller.K)*H);
                 % Contribution from excitation signal (r)
                 [outY_r, ~]                = size(Sy_r);
                 [GammaYPos_r, GammaYNeg_r] = obj.gammaR(problem,Sy_r);
                 ouLbCoeff_r                = size(GammaYPos_r,3);
                 Rtemp_r                    = [];
                 for m = 0:(obj.M-1) + (ouLbCoeff_r-1)
                     switchLag = 1;
                     for m1 = -(ouLbCoeff_r-1):(ouLbCoeff_r-1)
                         for m2 = -(obj.M-1):(obj.M-1)
                             if (m1+m2) == m
                                 if switchLag
                                     if isempty(Rtemp_r)
                                         if m1 >= 0
                                             Rtemp_r     = reshape(GammaYPos_r(:,:,abs(m1)+1)*vec(obj.parameters(:,:,abs(m2)+1)),outY_r,outY_r);
                                             switchLag = 0;
                                         else
                                             Rtemp_r     = reshape(GammaYNeg_r(:,:,abs(m1)+1)*vec(obj.parameters(:,:,abs(m2)+1)),outY_r,outY_r);
                                             switchLag = 0;
                                         end
                                     else
                                         if m1 >= 0
                                             Rtemp_r     = cat(3,Rtemp_r,reshape(GammaYPos_r(:,:,abs(m1)+1)*vec(obj.parameters(:,:,abs(m2)+1)),outY_r,outY_r));
                                             switchLag = 0;
                                         else
                                             Rtemp_r     = cat(3,Rtemp_r,reshape(GammaYNeg_r(:,:,abs(m1)+1)*vec(obj.parameters(:,:,abs(m2)+1)),outY_r,outY_r));
                                             switchLag = 0;
                                         end
                                     end
                                 else
                                     if m == 0
                                         if m1 >= 0
                                             Rtemp_r = Rtemp_r + reshape(GammaYPos_r(:,:,abs(m1)+1)*vec(obj.parameters(:,:,abs(m2)+1)),outY_r,outY_r);
                                         else
                                             Rtemp_r = Rtemp_r + reshape(GammaYNeg_r(:,:,abs(m1)+1)*vec(obj.parameters(:,:,abs(m2)+1)),outY_r,outY_r);
                                         end
                                     else
                                         if m1 >= 0
                                             Rtemp_r(:,:,m+1) = Rtemp_r(:,:,m+1) + reshape(GammaYPos_r(:,:,abs(m1)+1)*vec(obj.parameters(:,:,abs(m2)+1)),outY_r,outY_r);
                                         else
                                             Rtemp_r(:,:,m+1) = Rtemp_r(:,:,m+1) + reshape(GammaYNeg_r(:,:,abs(m1)+1)*vec(obj.parameters(:,:,abs(m2)+1)),outY_r,outY_r);
                                         end
                                     end
                                 end
                             end
                         end
                     end
                 end
                 % Contribution from noise signal (e)
                 [outY_e, ~]      = size(Sy_e);
                 [GammaYPos_e, ~] = obj.gammaR(problem,Sy_e);
                 ouLbCoeff_e      = size(GammaYPos_e,3);
                 Lambda           = problem.model_.NoiseVariance;
                 Rtemp_e          = [];
                 for m = 0:(ouLbCoeff_e-1)
                     if isempty(Rtemp_e)
                         Rtemp_e = reshape(GammaYPos_e(:,:,abs(m)+1)*vec(Lambda),outY_e,outY_e);
                     else
                         Rtemp_e = cat(3,Rtemp_e,reshape(GammaYPos_e(:,:,abs(m)+1)*vec(Lambda),outY_e,outY_e));
                     end
                 end
                 if size(Rtemp_r,3) > size(Rtemp_e,3)
                     Rtemp_e = cat(3,Rtemp_e,zeros(size(Rtemp_e,1),size(Rtemp_e,2),size(Rtemp_r,3)-size(Rtemp_e,3)));
                 end
                 if size(Rtemp_e,3) > size(Rtemp_r,3)
                     Rtemp_r = cat(3,Rtemp_r,zeros(size(Rtemp_r,1),size(Rtemp_r,2),size(Rtemp_e,3)-size(Rtemp_r,3)));
                 end
                 [AouLb,BBouLb,CouLb,DouLb] = obj.ssBound(Rtemp_r+Rtemp_e);
                 if ~issymmetric(obj.output.lb)
                     if max(max(abs(BouLb(:,:,1)-BouLb(:,:,1)'))) < 1e-8
                         warning('Lower bound on output is numerically non-Hermitian. The bound is made Hermitian.')
                         BouLb = 0.5*(BouLb+BouLb');
                         QouLb = sdpvar(size(AouLb,2),size(AouLb,2),'symmetric');
                         con   = [con; [QouLb-AouLb'*QouLb*AouLb  CouLb'-AouLb'*QouLb*BBouLb;
                             CouLb-BBouLb'*QouLb*AouLb  DouLb+DouLb'-BouLb-BBouLb'*QouLb*BBouLb] >= 0];
                     else
                         warning('Lower bound on output is non-Hermitian. A frequency-sampled element-wise constraint is used instead of the KYP lemma.')
                         wouLb = linspace(1e-3,pi);
                         for k = 1:length(wouLb);
                             PouLb = CouLb*((exp(1i*wouLb(k)*problem.model.Ts)*...
                                 eye(size(AouLb)) - AouLb)\eye(size(AouLb)))*BBouLb + DouLb;
                             PouLb = PouLb+PouLb';
                             for i = 1:size(PouLb,1)
                                 for j = 1:size(PouLb,2)
                                     con = [con; -obj.output.lb(i,j)+PouLb(i,j) >= 0];
                                 end
                             end
                         end
                     end
                 else
                     QouLb = sdpvar(size(AouLb,2),size(AouLb,2),'symmetric');
                     con   = [con; [QouLb-AouLb'*QouLb*AouLb  CouLb'-AouLb'*QouLb*BBouLb;
                         CouLb-BBouLb'*QouLb*AouLb  DouLb+DouLb'-obj.output.lb-BBouLb'*QouLb*BBouLb] >= 0];
                 end
             else
                 error('Lower bound on output spectrum is specified wrong.')
             end
         end
         % Upper bound output spectrum
         if ~isempty(obj.output.ub)
             % Sample LTI constraint
             if isa(obj.output.ub,'lti')
                 [BouUb,wouUb] = freqresp(obj.output.ub*obj.output.ub',...
                     linspace(1e-3,pi)./problem.model.Ts);
                 if max(max(max(abs(imag(BouUb)))))<1e-6
                     BouUb = real(BouUb);
                 end
             elseif isa(obj.output.ub,'struct')
                 BouUb = obj.output.ub.B;
                 wouUb = obj.output.ub.w;
             elseif size(obj.output.ub,1) == size(problem.model_,1) && size(obj.output.ub,2) == size(problem.model_,1)
                 wouUb = linspace(1e-3,pi)./problem.model.Ts;
                 BouUb = repmat(obj.output.ub,1,1,length(wouUb));
             end
             [~,m] = size(problem.model_);
             [G,H] = problem.model_.tf();
             Sy_r  = minreal(eye(m)/(eye(m) - G*problem.controller.K)*G);
             Sy_e  = minreal(eye(m)/(eye(m) - G*problem.controller.K)*H);
             % Contribution from excitation signal (r)
             [outY_r, ~]                = size(Sy_r);
             [GammaYPos_r, GammaYNeg_r] = obj.gammaR(problem,Sy_r);
             ouUbCoeff_r                = size(GammaYPos_r,3);
             Rtemp_r                    = [];
             for m = 0:(obj.M-1) + (ouUbCoeff_r-1)
                 switchLag = 1;
                 for m1 = -(ouUbCoeff_r-1):(ouUbCoeff_r-1)
                     for m2 = -(obj.M-1):(obj.M-1)
                         if (m1+m2) == m
                             if switchLag
                                 if isempty(Rtemp_r)
                                     if m1 >= 0
                                         Rtemp_r     = reshape(GammaYPos_r(:,:,abs(m1)+1)*vec(obj.parameters(:,:,abs(m2)+1)),outY_r,outY_r);
                                         switchLag = 0;
                                     else
                                         Rtemp_r     = reshape(GammaYNeg_r(:,:,abs(m1)+1)*vec(obj.parameters(:,:,abs(m2)+1)),outY_r,outY_r);
                                         switchLag = 0;
                                     end
                                 else
                                     if m1 >= 0
                                         Rtemp_r     = cat(3,Rtemp_r,reshape(GammaYPos_r(:,:,abs(m1)+1)*vec(obj.parameters(:,:,abs(m2)+1)),outY_r,outY_r));
                                         switchLag = 0;
                                     else
                                         Rtemp_r     = cat(3,Rtemp_r,reshape(GammaYNeg_r(:,:,abs(m1)+1)*vec(obj.parameters(:,:,abs(m2)+1)),outY_r,outY_r));
                                         switchLag = 0;
                                     end
                                 end
                             else
                                 if m == 0
                                     if m1 >= 0
                                         Rtemp_r = Rtemp_r + reshape(GammaYPos_r(:,:,abs(m1)+1)*vec(obj.parameters(:,:,abs(m2)+1)),outY_r,outY_r);
                                     else
                                         Rtemp_r = Rtemp_r + reshape(GammaYNeg_r(:,:,abs(m1)+1)*vec(obj.parameters(:,:,abs(m2)+1)),outY_r,outY_r);
                                     end
                                 else
                                     if m1 >= 0
                                         Rtemp_r(:,:,m+1) = Rtemp_r(:,:,m+1) + reshape(GammaYPos_r(:,:,abs(m1)+1)*vec(obj.parameters(:,:,abs(m2)+1)),outY_r,outY_r);
                                     else
                                         Rtemp_r(:,:,m+1) = Rtemp_r(:,:,m+1) + reshape(GammaYNeg_r(:,:,abs(m1)+1)*vec(obj.parameters(:,:,abs(m2)+1)),outY_r,outY_r);
                                     end
                                 end
                             end
                         end
                     end
                 end
             end
             % Contribution from noise signal (e)
             [outY_e, ~]      = size(Sy_e);
             [GammaYPos_e, ~] = obj.gammaR(problem,Sy_e);
             ouUbCoeff_e      = size(GammaYPos_e,3);
             Lambda           = problem.model_.NoiseVariance;
             Rtemp_e          = [];
             for m = 0:(ouUbCoeff_e-1)
                 if isempty(Rtemp_e)
                     Rtemp_e = reshape(GammaYPos_e(:,:,abs(m)+1)*vec(Lambda),outY_e,outY_e);
                 else
                     Rtemp_e = cat(3,Rtemp_e,reshape(GammaYPos_e(:,:,abs(m)+1)*vec(Lambda),outY_e,outY_e));
                 end
             end
             if size(Rtemp_r,3) > size(Rtemp_e,3)
                 Rtemp_e = cat(3,Rtemp_e,zeros(size(Rtemp_e,1),size(Rtemp_e,2),size(Rtemp_r,3)-size(Rtemp_e,3)));
             end
             if size(Rtemp_e,3) > size(Rtemp_r,3)
                 Rtemp_r = cat(3,Rtemp_r,zeros(size(Rtemp_r,1),size(Rtemp_r,2),size(Rtemp_e,3)-size(Rtemp_r,3)));
             end
             [AouUb,BBouUb,CouUb,DouUb] = obj.ssBound(Rtemp_r+Rtemp_e);
             for k = 1:length(wouUb);
                 PouUb = CouUb*((exp(1i*wouUb(k)*problem.model.Ts)*...
                     eye(size(AouUb)) - AouUb)\eye(size(AouUb)))*BBouUb + DouUb;
                 PouUb = PouUb+PouUb';
                 if isa(obj.output.ub,'struct') && ~strcmp(obj.output.ub.type,'lmi')
                     for i = 1:size(PouUb,1)
                         for j = 1:size(PouUb,2)
                             con = [con; BouUb(i,j,k)-PouUb(i,j) >= 0];
                         end
                     end
                 else
                     if ~issymmetric(BouUb(:,:,k))
                         if max(max(abs(BouUb(:,:,k)-BouUb(:,:,k)'))) < 1e-8
                             warning('Upper bound on output is numerically non-Hermitian. The bound is Hermitian.')
                             BouUb(:,:,k) = 0.5*(BouUb(:,:,k)+BouUb(:,:,k)');
                         else
                             warning('Upper bound on output is non-Hermitian. Linear matrix constraint is handled as an element-wise constraint.')
                         end
                     end
                     con = [con; BouUb(:,:,k)-PouUb >= 0];
                 end
             end
         end
         %% Power constraints
         % Bound excitation power
         if ~isempty(obj.excitation.power.ub)
             con = [con; obj.excitationPower <= obj.excitation.power.ub];
         end
         if ~isempty(obj.excitation.power.lb)
             con = [con; obj.excitationPower >= obj.excitation.power.lb];
         end
         
         % Bound input power (contribution from r and e)
         if ~isempty(obj.input.power.ub)
             Pu  = obj.inputPower(problem);
             con = [con; Pu.tot <= obj.input.power.ub];
         end
         if ~isempty(obj.input.power.lb)
             Pu  = obj.inputPower(problem);
             con = [con; Pu.tot >= obj.input.power.lb];
         end
         
         % Bound output power (contribution from r and e)
         if ~isempty(obj.output.power.ub)
             Py  = obj.outputPower(problem);
             con = [con; Py.tot <= obj.output.power.ub];
         end
         if ~isempty(obj.output.power.lb)
             Py  = obj.outputPower(problem);
             con = [con; Py.tot >= obj.output.power.lb];
         end
      end
      %% Additional functions
      function M = get.M(obj)
          sz = size(obj.parameters);
          if numel(sz) == 3
              M = sz(3);
          else
              M = 1;
          end
      end
      function n = get.n(obj)
          n = size(obj.parameters,1);
      end
      function obj = set.excitation(obj,excitation)
          % Check spectrum lower bound.
          if isempty(excitation.lb)
              obj.excitation.lb = [];
          else
              if isa(excitation.lb,'lti') && size(excitation.lb,1) == size(excitation.lb,2) &&...
                      size(excitation.lb,1) == obj.n
                  obj.excitation.lb = excitation.lb;
              elseif isa(excitation.lb,'double') && size(excitation.lb,1) == size(excitation.lb,2) &&...
                      size(excitation.lb,1) == obj.n
                  obj.excitation.lb      = excitation.lb;
              elseif  isfield(excitation.lb,'B') && isfield(excitation.lb,'w')
                  if size(excitation.lb.B,3) == 1
                      obj.excitation.lb.B(1,1,:) = excitation.lb.B;
                  else
                      obj.excitation.lb.B = excitation.lb.B;
                  end
                  obj.excitation.lb.w    = excitation.lb.w;
                  if isfield(excitation.lb,'type')
                      obj.excitation.lb.type = excitation.lb.type;
                  else
                      obj.excitation.lb.type = 'lmi';
                  end
              else
                  throw(MException('','VerifySpectrumConstraint:IncompatibleSpectrumConstraint'));
              end
          end
          
          % Check spectrum upper bound
          if isempty(excitation.ub)
              obj.excitation.ub = [];
          else
              if isa(excitation.ub,'lti') && size(excitation.ub,1) == size(excitation.ub,2) &&...
                      size(excitation.ub,1) == obj.n
                  obj.excitation.ub = excitation.ub;
              elseif isa(excitation.ub,'double') && size(excitation.ub,1) == size(excitation.ub,2) &&...
                      size(excitation.ub,1) == obj.n
                  obj.excitation.ub = excitation.ub;     
              elseif  isfield(excitation.ub,'B') && isfield(excitation.ub,'w')
                  if size(excitation.ub.B,3) == 1
                      obj.excitation.ub.B(1,1,:) = excitation.ub.B;
                  else
                      obj.excitation.ub.B = excitation.ub.B;
                  end
                  obj.excitation.ub.w    = excitation.ub.w;
                  if isfield(excitation.ub,'type')
                      obj.excitation.ub.type = excitation.ub.type;
                  else
                      obj.excitation.ub.type = 'lmi';
                  end
              elseif isempty(excitation.ub) 
                  obj.excitation.ub = [];
              else
                  throw(MException('','VerifySpectrumConstraint:IncompatibleSpectrumConstraint'));
              end
          end
          % Check that upper power constraint on excitation is real, scalar and positive.
          if isempty(excitation.power.ub)
              obj.excitation.power.ub = [];
          elseif isreal(excitation.power.ub) && numel(excitation.power.ub) == 1 && excitation.power.ub >= 0
              % Infinite upper bound treated as empty constraint
              if isinf(excitation.power.ub)
                  obj.excitation.power.ub = [];
              else
                  obj.excitation.power.ub = excitation.power.ub;
              end
          else
              throw(MException('','VerifySpectrumConstraint:IncompatiblePowerConstraint'));
          end
          
          % Check that lower power constraint on excitation is real, scalar and positive.
          if isempty(excitation.power.lb)
              obj.excitation.power.lb = [];
          elseif isreal(excitation.power.lb) && numel(excitation.power.lb) == 1 && excitation.power.lb >= 0
              % Zero lower bound treated as empty constraint
              if isinf(excitation.power.lb)
                  obj.excitation.power.lb = [];
              else
                  obj.excitation.power.lb = excitation.power.lb;
              end
          else
              throw(MException('','VerifySpectrumConstraint:IncompatiblePowerConstraint'));
          end
          
      end
      function obj = set.input(obj,input)
          % Check spectrum lower bound.
          if isempty(input.lb)
              obj.input.lb = [];
          else
              if isa(input.lb,'lti') && size(input.lb,1) == size(input.lb,2) &&...
                      size(input.lb,1) == obj.n
                  obj.input.lb = input.lb;
              elseif isa(input.lb,'double') && size(input.lb,1) == size(input.lb,2) &&...
                      size(input.lb,1) == obj.n
                  obj.input.lb = input.lb;
              elseif  isfield(input.lb,'B') && isfield(input.lb,'w')
                  if size(input.lb.B,3) == 1
                      obj.input.lb.B(1,1,:) = input.lb.B;
                  else
                      obj.input.lb.B = input.lb.B;
                  end
                  obj.input.lb.w    = input.lb.w;
                  if isfield(input.lb,'type')
                      obj.input.lb.type = input.lb.type;
                  else
                      obj.input.lb.type = 'lmi';
                  end
              elseif isempty(input.lb)
                  obj.input.lb = [];
              else
                  throw(MException('','VerifySpectrumConstraint:IncompatibleSpectrumConstraint'));
              end
          end
          
          % Check spectrum upper bound
          if isempty(input.ub)
              obj.input.ub = [];
          else
              if isa(input.ub,'lti') && size(input.ub,1) == size(input.ub,2) &&...
                      size(input.ub,1) == obj.n
                  obj.input.ub = input.ub;
              elseif isa(input.ub,'double') && size(input.ub,1) == size(input.ub,2) &&...
                      size(input.ub,1) == obj.n
                  obj.input.ub            = input.ub;
              elseif  isfield(input.ub,'B') && isfield(input.ub,'w')
                  if size(input.ub.B,3) == 1
                      obj.input.ub.B(1,1,:) = input.ub.B;
                  else
                      obj.input.ub.B = input.ub.B;
                  end
                  obj.input.ub.w    = input.ub.w;
                  if isfield(input.ub,'type')
                      obj.input.ub.type = input.ub.type;
                  else
                      obj.input.ub.type = 'lmi';
                  end
              elseif isempty(input.ub) 
                  obj.input.ub = [];
              else
                  throw(MException('','VerifySpectrumConstraint:IncompatibleSpectrumConstraint'));
              end
          end
          % Check that upper power constraint on input is real, scalar and positive.
          if isempty(input.power.ub)
              obj.input.power.ub = [];
          elseif isreal(input.power.ub) && numel(input.power.ub) == 1 && input.power.ub >= 0
              % Infinite upper bound treated as empty constraint
              if isinf(input.power.ub)
                  obj.input.power.ub = [];
              else
                  obj.input.power.ub = input.power.ub;
              end
          else
              throw(MException('','VerifySpectrumConstraint:IncompatiblePowerConstraint'));
          end
          
          % Check that lower power constraint on input is real, scalar and positive.
          if isempty(input.power.lb)
              obj.input.power.lb = [];
          elseif isreal(input.power.lb) && numel(input.power.lb) == 1 && input.power.lb >= 0
              % Zero lower bound treated as empty constraint
              if isinf(input.power.lb)
                  obj.input.power.lb = [];
              else
                  obj.input.power.lb = input.power.lb;
              end
          else
              throw(MException('','VerifySpectrumConstraint:IncompatiblePowerConstraint'));
          end
          
      end
      function obj = set.output(obj,output)
          % Check spectrum lower bound.
          if isempty(output.lb)
              obj.output.lb = [];
          else
              if isa(output.lb,'lti') && size(output.lb,1) == size(output.lb,2) 
                  obj.output.lb = output.lb;
              elseif isa(output.lb,'double') && size(output.lb,1) == size(output.lb,2)
                  obj.output.lb            = output.lb;
              elseif  isfield(output.lb,'B') && isfield(output.lb,'w')
                  if size(output.lb.B,3) == 1
                      obj.output.lb.B(1,1,:) = output.lb.B;
                  else
                      obj.output.lb.B = output.lb.B;
                  end
                  obj.output.lb.w          = output.lb.w;
                  if isfield(output.lb,'type')
                      obj.output.lb.type = output.lb.type;
                  else
                      obj.output.lb.type = 'lmi';
                  end
              elseif isempty(output.lb) 
                  obj.output.lb = [];
              else
                  throw(MException('','VerifySpectrumConstraint:IncompatibleSpectrumConstraint'));
              end
          end
          
          % Check spectrum upper bound
          if isempty(output.ub)
              obj.output.ub = [];
          else
              if isa(output.ub,'lti') && size(output.ub,1) == size(output.ub,2)
                  obj.output.ub = output.ub;
              elseif isa(output.ub,'double') && size(output.ub,1) == size(output.ub,2) 
                  obj.output.ub            = output.ub;
              elseif  isfield(output.ub,'B') && isfield(output.ub,'w')
                  if size(output.ub.B,3) == 1
                      obj.output.ub.B(1,1,:) = output.ub.B;
                  else
                      obj.output.ub.B = output.ub.B;
                  end
                  obj.output.ub.w    = output.ub.w;
                  if isfield(output.ub,'type')
                      obj.output.ub.type = output.ub.type;
                  else
                      obj.output.ub.type = 'lmi';
                  end
              elseif isempty(output.ub) 
                  obj.output.ub = [];
              else
                  throw(MException('','VerifySpectrumConstraint:IncompatibleSpectrumConstraint'));
              end
          end
          
          % Check that upper power constraint on output is real, scalar and positive.
          if isempty(output.power.ub)
              obj.output.power.ub = [];
          elseif isreal(output.power.ub) && numel(output.power.ub) == 1 && output.power.ub >= 0
              % Infinite upper bound treated as empty constraint
              if isinf(output.power.ub)
                  obj.output.power.ub = [];
              else
                  obj.output.power.ub = output.power.ub;
              end
          else
              throw(MException('','VerifySpectrumConstraint:IncompatiblePowerConstraint'));
          end
          
          % Check that lower power constraint on output is real, scalar and positive.
          if isempty(output.power.lb)
              obj.output.power.lb = [];
          elseif isreal(output.power.lb) && numel(output.power.lb) == 1 && output.power.lb >= 0
              % Zero lower bound treated as empty constraint
              if isinf(output.power.lb)
                  obj.output.power.lb = [];
              else
                  obj.output.power.lb = output.power.lb;
              end
          else
              throw(MException('','VerifySpectrumConstraint:IncompatiblePowerConstraint'));
          end
          
      end
      function obj = set.M(obj,M)
          if isnumeric(M) && M > 0
              obj.parameters = sdpvar(obj.n,obj.n,M);
          end
      end
      function obj = set.n(obj,n)
          if isnumeric(n) && n > 1
              obj.parameters = sdpvar(n,n,obj.M);
          end
      end
   end
end