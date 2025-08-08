function [triggerDelays, acqDuration] = calcTriggerDelays_general( ...
    grad0, probeType, probeRadius, probeT2star, signalCutoff)
    switch probeType
        case 'H'
          gamma = 2.675e5;   % [rad/T] gyromagnetic ratio for proton.
        case 'F'
          gamma = 2.518e5;   % [rad/T] gyromagnetic ratio for fluorine.
        otherwise
          error('Unsupported probe type: %s. Supported types: H, F.', probeType);
    end

    dt     = grad0.dt;
    t_grad = grad0.t_grad(1,:);
    nGrad  = grad0.nGrad;
    
    acqDuration       = -probeT2star * log(signalCutoff); 
    nPointAcqDuration = round(acqDuration/dt);
    acqDuration       = nPointAcqDuration*dt;
    % the max acq duration without considering gradient induced dephasing.
    
    t2s_model = @(t) exp(-t/probeT2star);
    t2sDecay = t2s_model(t_grad(1:nPointAcqDuration+1));
    % figure; plot(t2sDecay)
    %
    N  = 15;
    x   = linspace(-probeRadius, probeRadius,N)';
    y   = linspace(-probeRadius, probeRadius,N)';
    z   = linspace(-probeRadius, probeRadius,N)';
    [X, Y, Z] = meshgrid(x, y, z);
    D = sqrt(X.^2 + Y.^2 + Z.^2);
    insideMask = D <= probeRadius;
    pos = [X(insideMask), Y(insideMask), Z(insideMask)];
    
    M0 = [1; 0; 0];
    nSpin = size(pos, 1);
    
    triggerDelays = [];
    signals = {}; 
    s = 1;
    e = s+nPointAcqDuration;
    iter = 0;
    stop = false;
    while ~stop
        iter = iter+1;
        triggerDelays(end+1) = (s-1) * dt;
        t0 = linspace(s-1, e-1, nPointAcqDuration+1) * dt;
        if e>nGrad+1
            e=nGrad+1;
        end
        g = zeros(nPointAcqDuration+1, 3);
        g(1:e-s+1, :)  = grad0.shape(:, s:e)' / grad0.gamma * 1e3;
        
        [~, ~, Mxyt, ~, ~, ~] = blochsim_CK(zeros(size(g, 1), 1), ...
            g, pos, ones([nSpin 1]), zeros([nSpin 1]), ...
            'M0', M0, 'dt', dt, 'gamma', gamma);
        signal = abs(sum(Mxyt, 1))/nSpin .* t2sDecay;
        
        signal = interp1(t0, signal, t_grad, 'linear',0);
        signals{iter} = signal;
        mask = cell2mat(signals') > signalCutoff;
        indices = find(sum(mask, 1) == 0);
    
        if isempty(indices)
            clc; fprintf('\rFinished: %d segs\n', iter);
            stop = true;
        else
            clc; fprintf('\rProgress: %d, %d\n', iter, indices(1));
            s = indices(1) - 10;
            e = s+nPointAcqDuration;
        end
    end
    [fig] = plot_probeSignal(signals, t_grad, probeType, probeRadius, signalCutoff);
end