% Author: Mariette Annergren & Christian A. Larsson
% Copyright (c) 2015 Mariette Annergren & Christian A. Larsson

function [u,exitFlag] = simpleMPC(A,B,C,Nu,Ny,lambda,x0,up,r,Ub,Yb)
   % Construct matrices for evolution of system over prediction horizon as
   % function of input sequence.
   Psi          = zeros(size(C,1)*Ny,size(A,1));
   Upsilon      = zeros(size(C,1)*Ny,size(B,2)*Nu);
   Ga           = eye(size(B,2)*Nu,size(B,2)*Nu)...
                    -diag(ones(size(B,2)*(Nu-1),1),size(B,2));
   Up           = zeros(size(C,1)*Ny,size(B,2));
   for i = 1:Ny
      Psi(size(C,1)*(i-1)+1:size(C,1)*i,:) = C*A^(Ny-i+1);
      Up(size(C,1)*(i-1)+1:size(C,1)*i,:) = C*A^(Ny-i)*B;
   end
   for i = 2:Nu
      Upsilon(1:size(C,1)*(i+Ny-Nu),size(B,2)*(i-1)+1:size(B,2)*i) = ...
         Up(end-size(C,1)*(i+Ny-Nu)+1:end,:);
   end
   % Take care of the case when Ny > Nu
   for i = 1:Ny-Nu+1
      Upsilon(1:i*size(C,1),1:size(B,2)) = ...
         Upsilon(1:i*size(C,1),1:size(B,2)) + Up(end-size(C,1)*i+1:end,:);
   end
   
   % Reference
   R         = repmat(r,Ny,1);
   
   % Constraints
   A         = [];
   b         = [];
   if ~isinf(Ub(:,1))
      A = [A; eye(size(B,2)*Nu)];
      b = [b; repmat(Ub(:,1),Nu,1)];
   end
   if ~isinf(Ub(:,2))
      A = [A; -eye(size(B,2)*Nu)];
      b = [b; -repmat(Ub(:,2),Nu,1)];
   end
   if ~isinf(Yb(:,1))
      A = [A; Upsilon];
      b = [b; -Psi*x0 + repmat(Yb(:,1),Ny,1)];
   end
   if ~isinf(Yb(:,2))
      A = [A; -Upsilon];
      b = [b; Psi*x0 - repmat(Yb(:,2),Ny,1)];
   end
   
   % Setup QP
   f              = 2*(x0'*Psi'-R')*Upsilon...
                    - 2*[zeros(size(B,2)*(Nu-1),1);up]'*lambda(end-size(B,2)+1:end,:)*Ga;
   H              = (Upsilon'*Upsilon) + lambda*(Ga'*Ga);
   H              = (H+H')/2;
   
   % Solve optimization
   opts           = optimset('Display','off');
   [u,~,exitFlag] = quadprog(2*H,f,A,b,[],[],[],[],[],opts);
   if exitFlag == 1
      u = u(end);
   else
      u = 0;
   end
end