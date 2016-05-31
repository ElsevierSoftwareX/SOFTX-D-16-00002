% Representation of MA (FIR) spectra


% Author: Mariette Annergren & Christian A. Larsson
% Copyright (c) 2015 Mariette Annergren & Christian A. Larsson

classdef oidSpectrumMA < oidSpectrumFDP
   methods
      function obj = oidSpectrumMA(M,n)
         obj = obj@oidSpectrumFDP('MA',M,n);
      end
 
      function iF = informationMatrix(obj,prob)
         n = nparams(prob.model_);
         [p,m] = size(prob.model_);
         
         % Fixed controller or open loop
         if ~prob.controller.free
            % Construct Gamma_r and Gamma_e
            [G,H] = prob.model_.tf();
            [dG,dH] = prob.model_.gradient();
            Hinv  = eye(p)/H;
            Su    = eye(m)/(eye(m) - prob.controller.K*G);
            Gr = tf(zeros(n,p*m));
            Ge = tf(zeros(n,p*p));
            for i = 1:n
               Gr(i,:) = oidProblem.vecTF(minreal(Hinv*dG{i}*Su).').';
               Ge(i,:) = oidProblem.vecTF(...
                  minreal(Hinv*(dH{i}-dG{i}*Su*prob.controller.K*H)).').';
            end
            % Noise contribution
            Li = eye(p)/prob.model_.NoiseVariance;
            Fe = covar(Ge,chol(kron(Li,prob.model_.NoiseVariance)));
            
            % Reference contribution
            Rr = oidProblem.acf(oidProblem.vecTF(Gr),obj.M);
            CC = sdpvar((p*m)^2*(2*obj.M-1),1);
            RR = zeros(n^2,(m*p)^2*(2*obj.M-1));
            
            % C_0 (DC component)
            CC(1:(p*m)^2,:) = vec(kron(Li,obj.parameters(:,:,1)));
            for j = 1:(p*m)
               for k = 1:n
                  RR(1+(k-1)*n:k*n,1+(j-1)*(m*p):j*m*p) = ...
                     reshape(Rr(:,k+(j-1)*n,1),n,p*m);
               end
            end
            
            % C_k, k = 1,...,M 
            for i = 1:obj.M-1
                RT = Rr(:,:,i+1)';
                CC(1+i*(p*m)^2:(i+1)*(p*m)^2,:) = vec(...
                   kron(Li,obj.parameters(:,:,i+1)));
                CC(1+(obj.M+i-1)*(p*m)^2:(i+obj.M)*(p*m)^2,:) = vec(...
                   kron(Li,obj.parameters(:,:,i+1)'));
                for j = 1:(p*m)
                    for k = 1:n
                        RR(1+(k-1)*n:k*n,1+(j-1)*(m*p)+i*(m*p)^2:j*m*p+i*(m*p)^2)=...
                            reshape(Rr(:,k+(j-1)*n,i+1),n,p*m);
                        RR(1+(k-1)*n:k*n,1+(j-1)*(m*p)+(obj.M+i-1)*(m*p)^2:j*m*p+(i+obj.M-1)*(m*p)^2)=...
                            reshape(RT(:,k+(j-1)*n),n,p*m);
                    end
                end
            end
            Fr = reshape(RR*CC,n,n);
            
            % Total per sample information.
            iF = Fr + Fe;
         else
            error('MA (FIR) spectra currently only handle the fixed controller case');
         end
      end
      function [RRPos, RRNeg] = gammaR(obj,prob,S)
         
         [numOut,numIn] = size(S);
         
         % Fixed controller or open loop
         if ~prob.controller.free
             
            % Construct vectorized Gamma 
            Gr = oidProblem.vecTF(minreal(S)); % vectorized S, size: numOut*numIn x 1
            
            % Construct acf for vectorized Gamma
            coeffMax = 50;
            for Coeff = 1:coeffMax
                RrPos = oidProblem.acf(Gr,Coeff); % size: numOut*numIn x numOut*numIn x M
                if max(max(abs(RrPos(:,:,end)))) < 1e-8
                    break
                end
            end
            RrNeg = permute(RrPos,[2,1,3]);
            RRPos = zeros(numOut*numOut,numIn*numIn,size(RrPos,3)); % allocate memory for reshaped Rr
            RRNeg = zeros(numOut*numOut,numIn*numIn,size(RrNeg,3)); % allocate memory for reshaped Rr
            
            % Properly reshape Rr
            for kM = 1:Coeff
                index = 0;
                for kIn = 1:numIn
                    for kOut = 1:numOut
                        index = index + 1;
                        RRPos(1+(kOut-1)*numOut:kOut*numOut,1+(kIn-1)*numIn:kIn*numIn,kM) = ...
                            reshape(RrPos(:,index,kM),numOut,numIn);
                        RRNeg(1+(kOut-1)*numOut:kOut*numOut,1+(kIn-1)*numIn:kIn*numIn,kM) = ...
                            reshape(RrNeg(:,index,kM),numOut,numIn);
                    end
                end
            end
         else
            error('MA (FIR) spectra currently only handle the fixed controller case');
         end
      end
      function Pr = excitationPower(obj,prob)
         % Excitation power contribution 
         Pr = trace(obj.parameters(:,:,1));
      end
      function Pu = inputPower(obj,prob)
         % Input power contribution from r if closed loop and u if open
         % loop (Pr = Pu if open-loop)
         if ~(prob.controller.openLoop)
             if ~prob.controller.free
                 [~,m]  = size(prob.model_);
                 [G,H]  = prob.model_.tf();
                 Lambda = prob.model_.NoiseVariance;
                 % Contribution from r and e
                 Sr    = eye(m)/(eye(m) - prob.controller.K*G);
                 Se    = Sr*prob.controller.K*H;
                 Rr    = oidProblem.acf(tf(Sr),obj.M);
                 Re    = oidProblem.acf(tf(Se),1); % the noise is assumed to be white
                 Pu_r  = trace(Rr(:,:,1)*obj.parameters(:,:,1));
                 Pu_e  = trace(Re(:,:,1)*Lambda);
                 for m = 2:obj.M
                     Pu_r = Pu_r + 2*trace(Rr(:,:,m)*obj.parameters(:,:,m));
                 end
             else
                 error('MA (FIR) spectra currently only handle the fixed controller case');
             end
         else
             Pu_r = trace(obj.parameters(:,:,1));
             Pu_e = 0;
         end
         Pu.tot = Pu_r+Pu_e;
         Pu.r   = Pu_r;
         Pu.e   = Pu_e;
      end
      function Py = outputPower(obj,prob)
          % Output power contribution from r if closed loop and u if open
          % loop
          if ~prob.controller.free
              [m,~]  = size(prob.model_);
              [G,H]  = prob.model_.tf();
              Lambda = prob.model_.NoiseVariance;
              % Contribution from r and e
              Sr    = eye(m)/(eye(m) - G*prob.controller.K)*G;
              Se    = eye(m)/(eye(m) - G*prob.controller.K)*H;
              Rr    = oidProblem.acf((tf(Sr)).',obj.M);
              Re    = oidProblem.acf((tf(Se)).',1); % the noise is assumed to be white
              Py_r  = trace(Rr(:,:,1)*obj.parameters(:,:,1));
              Py_e  = trace(Re(:,:,1)*Lambda);
              for m = 2:obj.M
                  Py_r = Py_r + 2*trace(Rr(:,:,m)*obj.parameters(:,:,m));
              end
          else
              error('MA (FIR) spectra currently only handle the fixed controller case');
          end
         Py.tot = Py_r+Py_e;
         Py.r   = Py_r;
         Py.e   = Py_e;
      end
      function [A,B,C,D] = ss(obj)
         if obj.M < 2
            A = 0;
            B = zeros(1,size(obj.parameters,1));
            C = zeros(size(obj.parameters,1),1);
            D = 0.5*obj.parameters;
         else
            A = zeros(obj.n*(obj.M-1));
            A(obj.n+1:obj.n*(obj.M-1),1:obj.n*(obj.M-2)) = eye(obj.n*(obj.M-2));
            B = [eye(obj.n); zeros(obj.n*(obj.M-2),obj.n)];
            C = reshape(obj.parameters(:,:,2:end),obj.n,obj.n*(obj.M-1));
            D = 0.5*obj.parameters(:,:,1);
         end
      end
      function [A,B,C,D] = ssBound(obj,R)
          n = size(R,2);
          if length(size(R)) == 3
              M = size(R,3);
          else
              M = 1;
          end
         if M < 2
            A = 0;
            B = zeros(1,n);
            C = zeros(n,1);
            D = 0.5*R;
         else
            A = zeros(n*(M-1));
            A(n+1:n*(M-1),1:n*(M-2)) = eye(n*(M-2));
            B = [eye(n); zeros(n*(M-2),n)];
%             C = reshape(R(:,:,2:end),n,n*(M-1));
            C = R(:,:,2);
            for index = 3:M
                C = [C R(:,:,index)];
            end
            D = 0.5*R(:,:,1);
         end
      end
      function H = spectralFactor(obj,Ts)
         if obj.M < 2
            H = real(sqrtm(double(obj.parameters)))*tf(1,1,Ts);
         else
            [A,B,C,D] = obj.ss();
            C = double(C);
            D = double(D);
            R = -D-D';
            S = -B;
            [P,~,~,rep] = dare(A',C',0*A,R,S);
            if rep == -1
               % "Add some noise" to system to get Symplectic spectrum away
               % from unit circle
               P = dare(A',C',0*A,R-(1e-10)*eye(size(R)),S);
            end
            L = D+D'-C*P*C';
            K = -(A*P*C'-B)*pinv(L);
            
            % This is the theoretical spectral factor, it is not used due to numerical
            % issues. An ad hoc method for cleaning up the numerics is implemented as
            % optH.
            % optH = tf(ss(A,K,C,eye(size(C,1)),Ts))*sqrtm(L);
            H = tf();
            for i = 1:obj.n
               [b,a] = ss2tf(A,K,C,eye(size(C,1)),i);
               b(abs(b)<1e-12) = 0;
               b = b(:,1:obj.M);
               b = mat2cell(b,ones(1,obj.n),obj.M);
               H(:,i) = tf(b,a(1:obj.M),Ts);
            end
            H = minreal(H*real(sqrtm(L)));
         end
      end
   end
end