% Author: Mariette Annergren & Christian A. Larsson
% Copyright (c) 2015 Mariette Annergren & Christian A. Larsson

if length(theta0) == 2
    figure(1);
    hold on
    ellipse(gamma/2*VappHessian,theta0,'r','empty');
    ellipse(Nident/chi2inv(alpha,2)*(infoMatrix),theta0,'k');
    plot(thetaHatMC(1,:),thetaHatMC(2,:),'xk');
    xlabel('\theta_1')
    ylabel('\theta_2')
    legend('Application set','Identification set','True parameters', 'Estimates')
end

figure(2);
hold on
plot(distTheta0)
xlabel('MC-run')
ylabel('Euclidean distance to \theta^0')

% Check of excitation
figure(3)
suptitle('\Phi_r (real part)')
hold on
constraintLB = 'lmi';
constraintUB = 'lmi';
if ~isempty(optInputDesign.spectrum.excitation.lb)
    if isa(optInputDesign.spectrum.excitation.lb,'lti')
        [magLB] = (freqresp((optInputDesign.spectrum.excitation.lb*optInputDesign.spectrum.excitation.lb'),wOut));
        index   = 0;
        for i = 1:size(trueSystem,2)
            for j = 1:size(trueSystem,2)
                index = index + 1;
                subplot(size(trueSystem,2),size(trueSystem,2),index); hold on;
                plot(wOut',squeeze(magLB(i,j,:)),'r-^') % plots only real part
            end
        end
    elseif isa(optInputDesign.spectrum.excitation.lb,'struct')
        index = 0;
        magLB = (optInputDesign.spectrum.excitation.lb.B);
        for i = 1:size(trueSystem,2)
            for j = 1:size(trueSystem,2)
                index = index + 1;
                subplot(size(trueSystem,2),size(trueSystem,2),index); hold on;
                plot(squeeze(optInputDesign.spectrum.excitation.lb.w),squeeze(magLB(i,j,:)),'r-^') % plots only real part
            end
        end
        if ~strcmp(optInputDesign.spectrum.excitation.lb.type,'lmi')
            constraintLB = 'element-wise';
        end
    else
        index = 0;
        magLB = repmat(optInputDesign.spectrum.excitation.lb,1,1,length(wOut));
        for i = 1:size(trueSystem,2)
            for j = 1:size(trueSystem,2)
                index = index + 1;
                subplot(size(trueSystem,2),size(trueSystem,2),index); hold on;
                plot(wOut,optInputDesign.spectrum.excitation.lb(i,j)*ones(length(wOut),1),'r-^') % plots only real part
            end
        end
    end
else
    [magLB] = (freqresp(tf(zeros(size(trueSystem,2),size(trueSystem,2))),wOut));
    index     = 0;
    for i = 1:size(trueSystem,2)
        for j = 1:size(trueSystem,2)
            index = index + 1;
            subplot(size(trueSystem,2),size(trueSystem,2),index); hold on;
            plot(wOut,squeeze(magLB(i,j,:))','r-^')  % plots only real part
        end
    end
end
if ~isempty(optInputDesign.spectrum.excitation.ub)
    if isa(optInputDesign.spectrum.excitation.ub,'lti')
        [magUB] = (freqresp(((optInputDesign.spectrum.excitation.ub*optInputDesign.spectrum.excitation.ub')),wOut));
        index   = 0;
        for i = 1:size(trueSystem,2)
            for j = 1:size(trueSystem,2)
                index = index + 1;
                subplot(size(trueSystem,2),size(trueSystem,2),index); hold on;
                plot(wOut,squeeze(magUB(i,j,:))','r-V')  % plots only real part
            end
        end
    elseif isa(optInputDesign.spectrum.excitation.ub,'struct')
        index = 0;
        magUB = (optInputDesign.spectrum.excitation.ub.B);
        for i = 1:size(trueSystem,2)
            for j = 1:size(trueSystem,2)
                index = index + 1;
                subplot(size(trueSystem,2),size(trueSystem,2),index); hold on;
                plot(squeeze(optInputDesign.spectrum.excitation.ub.w),squeeze(optInputDesign.spectrum.excitation.ub.B(i,j,:)),'r-V')  % plots only real part
            end
        end
        if ~strcmp(optInputDesign.spectrum.excitation.ub.type,'lmi')
            constraintUB = 'element-wise';
        end
    else
        index = 0;
        magUB = repmat(optInputDesign.spectrum.excitation.ub,1,1,length(wOut));
        for i = 1:size(trueSystem,2)
            for j = 1:size(trueSystem,2)
                index = index + 1;
                subplot(size(trueSystem,2),size(trueSystem,2),index); hold on;
                plot(wOut,optInputDesign.spectrum.excitation.ub(i,j)*ones(length(wOut),1),'r-V')  % plots only real part
            end
        end
    end
    [magTrue] = (freqresp((optH*optH'),wOut)); % is real
    index     = 0;
    for i = 1:size(squeeze(magTrue(:,:,1)),2)
        for j = 1:size(squeeze(magTrue(:,:,1)),2)
            index = index + 1;
            subplot(size(squeeze(magTrue(:,:,1)),2),size(squeeze(magTrue(:,:,1)),2),index); hold on;
            plot(wOut,squeeze(magTrue(i,j,:))','-*b')
%             legendsVec = get(gca,'children');
            if strcmp(constraintLB,'lmi') && strcmp(constraintUB,'lmi')
                legend('lower bound (lmi)', 'upper bound (lmi)','optimal \Phi_r')
            elseif strcmp(constraintLB,'lmi') && ~strcmp(constraintUB,'lmi')
                legend('lower bound (lmi)', 'upper bound (element-wise)','optimal \Phi_r')
            elseif ~strcmp(constraintLB,'lmi') && strcmp(constraintUB,'lmi')
                legend('lower bound (element-wise)', 'upper bound (lmi)','optimal \Phi_r')
            else
                legend('lower bound (element-wise)', 'upper bound (element-wise)','optimal \Phi_r')
            end
        end
    end
else
    [magTrue] = (freqresp((optH*optH'),wOut)); % is real
    index     = 0;
    for i = 1:size(trueSystem,2)
        for j = 1:size(trueSystem,2)
            index = index + 1;
            subplot(size(trueSystem,2),size(trueSystem,2),index); hold on;
            plot(wOut,squeeze(magTrue(i,j,:)),'-*b')
            if strcmp(constraintLB,'lmi')
                legend('lower bound (lmi)','optimal \Phi_r')
            else
                legend('lower bound (element-wise)','optimal \Phi_r')
            end
        end
    end
end
r             = lsim(optH,randn(10000,size(trueSystem,2)));
powerExOpt    = signalPower.excitation
powerExCheck1 = sum(var(r))
powerExCheck2 = 0;
for i = 1:size(trueSystem,2)
    powerExCheck2 = powerExCheck2 + trapz(wOut,real(magTrue(i,i,:)))/pi;
end
powerExCheck2
powerConCheck = 0;
for i = 1:size(trueSystem,2)
    powerConCheck = powerConCheck + trapz(wOut,real(magLB(i,i,:)))/pi;
end
powerConCheck

% Check of eigenvalues related to excitation spectra
if ~isempty(optInputDesign.spectrum.excitation.lb)
    eigenValuesExcitation_LowerBound = zeros(size(trueSystem,2),size(magTrue,3));
    for index = 1:size(magTrue,3)
        eigenValuesExcitation_LowerBound(:,index) = (eig(squeeze((magTrue(:,:,index)))-squeeze((magLB(:,:,index)))));
        figure(4); hold on; plot(real(eigenValuesExcitation_LowerBound(:,index)),imag(eigenValuesExcitation_LowerBound(:,index)),'*')
    end
    title('Eigenvalues of \Phi_r-lower bound for different \omega')
    xlabel('real part (should be non-negative)')
    ylabel('imaginary part (should be zero)')
end
if ~isempty(optInputDesign.spectrum.excitation.ub)
    eigenValuesExcitation_UpperBound = zeros(size(trueSystem,2),size(magTrue,3));
    for index = 1:size(magTrue,3)
        eigenValuesExcitation_UpperBound(:,index) = (eig(squeeze((magTrue(:,:,index)))-squeeze((magUB(:,:,index)))));
        figure(5); hold on; plot(real(eigenValuesExcitation_UpperBound(:,index)),imag(eigenValuesExcitation_UpperBound(:,index)),'*')
    end
    title('Eigenvalues of \Phi_r-upper bound for different \omega')
    xlabel('real part (should be non-negative)')
    ylabel('imaginary part (should be zero)')
end

eigenValuesExcitation = zeros(size(trueSystem,2),size(magTrue,3));
for index = 1:size(magTrue,3)
    eigenValuesExcitation(:,index) = (eig(squeeze((magTrue(:,:,index)))));
    figure(6); hold on; plot(real(eigenValuesExcitation(:,index)),imag(eigenValuesExcitation(:,index)),'*')
end
xlabel('real part (should be non-negative)')
ylabel('imaginary part (should be zero)')
title('Eigenvalues of \Phi_r for different \omega')


% Check of input
figure(7);
suptitle('\Phi_u (real part)')
hold on
G  = minreal(tf(model));
m  = size(G,2);
Su = minreal(eye(m)/(eye(m) - optInputDesign.controller.K*G));
constraintLB = 'lmi';
constraintUB = 'lmi';
if ~isempty(optInputDesign.spectrum.input.lb)
    if isa(optInputDesign.spectrum.input.lb,'lti')
        [magLB] = (freqresp((optInputDesign.spectrum.input.lb*optInputDesign.spectrum.input.lb'),wOut));
        index   = 0;
        for i = 1:size(trueSystem,2)
            for j = 1:size(trueSystem,2)
                index = index + 1;
                subplot(size(trueSystem,2),size(trueSystem,2),index); hold on;
                plot(wOut',squeeze(magLB(i,j,:)),'r-^')
            end
        end
    elseif isa(optInputDesign.spectrum.input.lb,'struct')
        index = 0;
        magLB = (optInputDesign.spectrum.input.lb.B);
        for i = 1:size(trueSystem,2)
            for j = 1:size(trueSystem,2)
                index = index + 1;
                subplot(size(trueSystem,2),size(trueSystem,2),index); hold on;
                plot(squeeze(optInputDesign.spectrum.input.lb.w),squeeze(optInputDesign.spectrum.input.lb.B(i,j,:)),'r-^')
            end
        end
        if ~strcmp(optInputDesign.spectrum.input.lb.type,'lmi')
            constraintLB = 'element-wise';
        end
    else
        index = 0;
        magLB = repmat(optInputDesign.spectrum.input.lb,1,1,length(wOut));
        for i = 1:size(trueSystem,2)
            for j = 1:size(trueSystem,2)
                index = index + 1;
                subplot(size(trueSystem,2),size(trueSystem,2),index); hold on;
                plot(wOut,optInputDesign.spectrum.input.lb(i,j)*ones(length(wOut),1),'r-^')
            end
        end
    end
else
    [magLB] = (freqresp(tf(zeros(size(trueSystem,2),size(trueSystem,2))),wOut));
    index     = 0;
    for i = 1:size(trueSystem,2)
        for j = 1:size(trueSystem,2)
            index = index + 1;
            subplot(size(trueSystem,2),size(trueSystem,2),index); hold on;
            plot(wOut,squeeze(magLB(i,j,:))','r-^')
        end
    end
end
if ~isempty(optInputDesign.spectrum.input.ub)
    if isa(optInputDesign.spectrum.input.ub,'lti')
        [magUB] = (freqresp(((optInputDesign.spectrum.input.ub*optInputDesign.spectrum.input.ub')),wOut));
        index     = 0;
        for i = 1:size(trueSystem,2)
            for j = 1:size(trueSystem,2)
                index = index + 1;
                subplot(size(trueSystem,2),size(trueSystem,2),index); hold on;
                plot(wOut,squeeze(magUB(i,j,:))','r-V')
            end
        end
    elseif isa(optInputDesign.spectrum.input.ub,'struct')
        index = 0;
        magUB = (optInputDesign.spectrum.input.ub.B);
        for i = 1:size(trueSystem,2)
            for j = 1:size(trueSystem,2)
                index = index + 1;
                subplot(size(trueSystem,2),size(trueSystem,2),index); hold on;
                plot(squeeze(optInputDesign.spectrum.input.ub.w),squeeze(optInputDesign.spectrum.input.ub.B(i,j,:)),'r-V')
            end
        end
        if ~strcmp(optInputDesign.spectrum.input.ub.type,'lmi')
            constraintUB = 'element-wise';
        end
    else
        index = 0;
        magUB = repmat(optInputDesign.spectrum.input.ub,1,1,length(wOut));
        for i = 1:size(trueSystem,2)
            for j = 1:size(trueSystem,2)
                index = index + 1;
                subplot(size(trueSystem,2),size(trueSystem,2),index); hold on;
                plot(wOut,optInputDesign.spectrum.input.ub(i,j)*ones(length(wOut),1),'r-V')
            end
        end
    end
    [magTrue] = (freqresp(minreal(Su*(optH*optH')*Su'),wOut));
    index     = 0;
    for i = 1:size(squeeze(magTrue(:,:,1)),2)
        for j = 1:size(squeeze(magTrue(:,:,1)),2)
            index = index + 1;
            subplot(size(squeeze(magTrue(:,:,1)),2),size(squeeze(magTrue(:,:,1)),2),index); hold on;
            plot(wOut,squeeze(magTrue(i,j,:))','-*b')
            if strcmp(constraintLB,'lmi') && strcmp(constraintUB,'lmi')
                legend('lower bound (lmi)', 'upper bound (lmi)','optimal \Phi_u')
            elseif strcmp(constraintLB,'lmi') && ~strcmp(constraintUB,'lmi')
                legend('lower bound (lmi)', 'upper bound (element-wise)','optimal \Phi_u')
            elseif ~strcmp(constraintLB,'lmi') && strcmp(constraintUB,'lmi')
                legend('lower bound (element-wise)', 'upper bound (lmi)','optimal \Phi_u')
            else
                legend('lower bound (element-wise)', 'upper bound (element-wise)','optimal \Phi_u')
            end
        end
    end
else
    [magTrue] = (freqresp(minreal(Su*(optH*optH')*Su'),wOut));
    index     = 0;
    for i = 1:size(trueSystem,2)
        for j = 1:size(trueSystem,2)
            index = index + 1;
            subplot(size(trueSystem,2),size(trueSystem,2),index); hold on;
            plot(wOut,squeeze(magTrue(i,j,:)),'-*b')
            if strcmp(constraintLB,'lmi')
                legend('lower bound (lmi)','optimal \Phi_u')
            else
                legend('lower bound (element-wise)','optimal \Phi_u')
            end
        end
    end
end
    
u_r                 = lsim(Su*optH,randn(10000,size(trueSystem,2)));
powerInputOpt_r     = value(signalPower.input.r)
powerInput_r_Check1 = sum(var(u_r))
powerInput_r_Check2 = 0;
for i = 1:size(trueSystem,2)
    powerInput_r_Check2 = powerInput_r_Check2 + trapz(wOut,real(magTrue(i,i,:)))/pi;
end
powerInput_r_Check2
[~, H]              = optInputDesign.model_.tf();
u_e                 = lsim(Su*optInputDesign.controller.K*H,randn(10000,size(trueSystem,2))*chol(lambda));
powerInputOpt_e     = value(signalPower.input.e)
powerInput_e_Check1 = sum(var(u_e))
powerConCheck = 0;
for i = 1:size(trueSystem,2)
    powerConCheck = powerConCheck + trapz(wOut,real(magLB(i,i,:)))/pi;
end
powerConCheck

% Check of eigenvalues related to input spectra
if ~isempty(optInputDesign.spectrum.input.lb)
    eigenValuesInput_LowerBound = zeros(size(trueSystem,2),size(magTrue,3));
    for index = 1:size(magTrue,3)
        eigenValuesInput_LowerBound(:,index) = (eig(squeeze((magTrue(:,:,index)))-squeeze((magLB(:,:,index)))));
        figure(8); hold on; plot(real(eigenValuesInput_LowerBound(:,index)),imag(eigenValuesInput_LowerBound(:,index)),'*')
    end
    plot(0.9:0.1:size(trueSystem,2)+0.1,zeros(length(0.9:0.1:size(trueSystem,2)+0.1),1),'k')
    title('Eigenvalues of \Phi_u-lower bound for different \omega')
    xlabel('real part (should be non-negative)')
    ylabel('imaginary part (should be zero)')
end
if ~isempty(optInputDesign.spectrum.input.ub)
    eigenValuesInput_UpperBound = zeros(size(trueSystem,2),size(magTrue,3));
    for index = 1:size(magTrue,3)
        eigenValuesInput_UpperBound(:,index) = (eig(squeeze((magTrue(:,:,index)))-squeeze((magUB(:,:,index)))));
        figure(9); hold on; plot(real(eigenValuesInput_UpperBound(:,index)),imag(eigenValuesInput_UpperBound(:,index)),'*')
    end
    title('Eigenvalues of \Phi_u-upper bound for different \omega')
    xlabel('real part (should be non-negative)')
    ylabel('imaginary part (should be zero)')
end

eigenValuesInput = zeros(size(trueSystem,2),size(magTrue,3));
for index = 1:size(magTrue,3)
    eigenValuesInput(:,index) = (eig(squeeze((magTrue(:,:,index)))));
    figure(10); hold on; plot(real(eigenValuesInput(:,index)),imag(eigenValuesInput(:,index)),'*')
end
xlabel('real part (should be non-negative)')
ylabel('imaginary part (should be zero)')
title('Eigenvalues of \Phi_u for different \omega')


% Check of output
figure(11);
suptitle('\Phi_y (real part)')
hold on
m = size(G,1);
Sy = (eye(m)/(eye(m) - G*optInputDesign.controller.K)*G);
constraintLB = 'lmi';
constraintUB = 'lmi';
if ~isempty(optInputDesign.spectrum.output.lb)
    if isa(optInputDesign.spectrum.output.lb,'lti')
        [magLB] = (freqresp((optInputDesign.spectrum.output.lb*optInputDesign.spectrum.output.lb'),wOut));
        index   = 0;
        for i = 1:size(trueSystem,1)
            for j = 1:size(trueSystem,1)
                index = index + 1;
                subplot(size(trueSystem,1),size(trueSystem,1),index); hold on;
                plot(wOut',squeeze(magLB(i,j,:)),'r-^')
            end
        end
    elseif isa(optInputDesign.spectrum.output.lb,'struct')
        index = 0;
        magLB = optInputDesign.spectrum.output.lb.B;
        for k = 1:size(magLB,3)
            if max(max(abs(magLB(:,:,k)-magLB(:,:,k)'))) < 1e-8
                magLB(:,:,k) = 0.5*(magLB(:,:,k)+magLB(:,:,k)');
            end
        end
        for i = 1:size(trueSystem,1)
            for j = 1:size(trueSystem,1)
                index = index + 1;
                subplot(size(trueSystem,1),size(trueSystem,1),index); hold on;
                plot(squeeze(optInputDesign.spectrum.output.lb.w),squeeze(optInputDesign.spectrum.output.lb.B(i,j,:)),'r-^')
            end
        end
        if ~strcmp(optInputDesign.spectrum.output.lb.type,'lmi')
            constraintLB = 'element-wise';
        end
    else
        index = 0;
        magLB = repmat(optInputDesign.spectrum.output.lb,1,1,length(wOut));
        for i = 1:size(trueSystem,1)
            for j = 1:size(trueSystem,1)
                index = index + 1;
                subplot(size(trueSystem,1),size(trueSystem,1),index); hold on;
                plot(wOut,optInputDesign.spectrum.output.lb(i,j)*ones(length(wOut),1),'r-^')
            end
        end
    end
else
    [magLB] = (freqresp(tf(zeros(size(trueSystem,1),size(trueSystem,1))),wOut));
    index     = 0;
    for i = 1:size(trueSystem,1)
        for j = 1:size(trueSystem,1)
            index = index + 1;
            subplot(size(trueSystem,1),size(trueSystem,1),index); hold on;
            plot(wOut,squeeze(magLB(i,j,:))','r-^')
        end
    end
end
if ~isempty(optInputDesign.spectrum.output.ub)
    if isa(optInputDesign.spectrum.output.ub,'lti')
        [magUB] = (freqresp(((optInputDesign.spectrum.output.ub*optInputDesign.spectrum.output.ub')),wOut));
        index   = 0;
        for i = 1:size(trueSystem,1)
            for j = 1:size(trueSystem,1)
                index = index + 1;
                subplot(size(trueSystem,1),size(trueSystem,1),index); hold on;
                plot(wOut,squeeze(magUB(i,j,:))','r-V')
            end
        end
    elseif isa(optInputDesign.spectrum.output.ub,'struct')
        index = 0;
        magUB = (optInputDesign.spectrum.output.ub.B);
        for i = 1:size(trueSystem,1)
            for j = 1:size(trueSystem,1)
                index = index + 1;
                subplot(size(trueSystem,1),size(trueSystem,1),index); hold on;
                plot(squeeze(optInputDesign.spectrum.output.ub.w),squeeze(optInputDesign.spectrum.output.ub.B(i,j,:)),'r-V')
            end
        end
        if ~strcmp(optInputDesign.spectrum.output.ub.type,'lmi')
            constraintUB = 'element-wise';
        end
    else
        index = 0;
        magUB = repmat(optInputDesign.spectrum.output.ub,1,1,length(wOut));
        for i = 1:size(trueSystem,1)
            for j = 1:size(trueSystem,1)
                index = index + 1;
                subplot(size(trueSystem,1),size(trueSystem,1),index); hold on;
                plot(wOut,optInputDesign.spectrum.output.ub(i,j)*ones(length(wOut),1),'r-V')
            end
        end
    end
    [magTrue] = (freqresp(minreal(Sy*(optH*optH')*Sy'),wOut));
    index     = 0;
    for i = 1:size(squeeze(magTrue(:,:,1)),1)
        for j = 1:size(squeeze(magTrue(:,:,1)),1)
            index = index + 1;
            subplot(size(squeeze(magTrue(:,:,1)),1),size(squeeze(magTrue(:,:,1)),1),index); hold on;
            plot(wOut,squeeze(magTrue(i,j,:))','-*b')
%             legendsVec = get(gca,'children');
            if strcmp(constraintLB,'lmi') && strcmp(constraintUB,'lmi')
                legend('lower bound (lmi)', 'upper bound (lmi)','optimal \Phi_y')
            elseif strcmp(constraintLB,'lmi') && ~strcmp(constraintUB,'lmi')
                legend('lower bound (lmi)', 'upper bound (element-wise)','optimal \Phi_y')
            elseif ~strcmp(constraintLB,'lmi') && strcmp(constraintUB,'lmi')
                legend('lower bound (element-wise)', 'upper bound (lmi)','optimal \Phi_y')
            else
                legend('lower bound (element-wise)', 'upper bound (element-wise)','optimal \Phi_y')
            end
        end
    end
else
    [magTrue] = (freqresp((Sy*(optH*optH')*Sy'),wOut));
    index     = 0;
    for i = 1:size(trueSystem,1)
        for j = 1:size(trueSystem,1)
            index = index + 1;
            subplot(size(trueSystem,1),size(trueSystem,1),index); hold on;
            plot(wOut,squeeze(magTrue(i,j,:)),'-*b')
            if strcmp(constraintLB,'lmi')
                legend('lower bound (lmi)','optimal \Phi_y')
            else
                legend('lower bound (element-wise)','optimal \Phi_y')
            end
        end
    end
end
    
u_r                  = lsim(Sy*optH,randn(10000,size(trueSystem,2)));
powerOutputOpt_r     = value(signalPower.output.r)
powerOutput_r_Check1 = sum(var(u_r))
powerOutput_r_Check2 = 0;
for i = 1:size(trueSystem,1)
    powerOutput_r_Check2 = powerOutput_r_Check2 + trapz(wOut,real(magTrue(i,i,:)))/pi;
end
powerOutput_r_Check2
u_e                  = lsim(eye(m)/(eye(m) - G*optInputDesign.controller.K)*H,randn(10000,size(trueSystem,2))*chol(lambda));
powerOutput_e_Opt    = value(signalPower.output.e)
powerOutput_e_Check1 = sum(var(u_e))
powerConCheck        = 0;
for i = 1:size(trueSystem,2)
    powerConCheck = powerConCheck + trapz(wOut,real(magLB(i,i,:)))/pi;
end
powerConCheck

% Check of eigenvalues related to output spectra
if ~isempty(optInputDesign.spectrum.output.lb)
    eigenValuesOutput_LowerBound = zeros(size(trueSystem,2),size(magTrue,3));
    for index = 1:size(magTrue,3)
        eigenValuesOutput_LowerBound(:,index) = (eig(squeeze((magTrue(:,:,index)))-squeeze((magLB(:,:,index)))));
        figure(12); hold on; plot(real(eigenValuesOutput_LowerBound(:,index)),imag(eigenValuesOutput_LowerBound(:,index)),'*')
    end
    title('Eigenvalues of \Phi_y-lower bound for different \omega')
    xlabel('real part (should be non-negative)')
    ylabel('imaginary part (should be zero)')
end
if ~isempty(optInputDesign.spectrum.output.ub)
    eigenValuesOutput_UpperBound = zeros(size(trueSystem,2),size(magTrue,3));
    for index = 1:size(magTrue,3)
        eigenValuesOutput_UpperBound(:,index) = (eig(squeeze((magTrue(:,:,index)))-squeeze((magUB(:,:,index)))));
        figure(13); hold on; plot(real(eigenValuesOutput_UpperBound(:,index)),imag(eigenValuesOutput_UpperBound(:,index)),'*')
    end
    title('Eigenvalues of \Phi_y-upper bound for different \omega')
    xlabel('real part (should be non-negative)')
    ylabel('imaginary part (should be zero)')
end

eigenValuesOutput = zeros(size(trueSystem,2),size(magTrue,3));
for index = 1:size(magTrue,3)
    eigenValuesOutput(:,index) = (eig(squeeze((magTrue(:,:,index)))));
    figure(14); hold on; plot(real(eigenValuesOutput(:,index)),imag(eigenValuesOutput(:,index)),'*')
end
xlabel('real part (should be non-negative)')
ylabel('imaginary part (should be zero)')
title('Eigenvalues of \Phi_y for different \omega')

