function [fig] = plot_stitch_multishot(grad0, triggerDelays, probeType, probeRadius, signalCutoff, nsegs2measure)
%PLOT_STITCH 
%   
    color_facecolor = "#1F1F1F";
    color_label     = "#CCCCCC";
    font_label      = "Arial";
    fig_width       = 800;
    fig_height      = 800;
    position = [(1920-fig_width)/2, (1080-fig_height)/2, fig_width, fig_height];

    figname = 'Stitch Pattern';
    fig = figure('Name', figname, 'Position', position, 'Color', color_facecolor);
    subplot(1,1,1, 'color', color_facecolor, 'xcolor', color_label, 'ycolor', color_label, 'fontname', font_label), hold on;

    [nShot, nDim, nPoint] = size(grad0.shape_all);

    dt = grad0.dt;
    
    for ishot = 1:nShot
        grad     = squeeze(grad0.shape_all(ishot, 1:2, :));
        ktraj    = calc_ktraj_from_grad(grad / grad0.gamma, probeType, false, dt);
        ktrajAbs = sqrt(ktraj(1,:).^2+ ktraj(2,:).^2);
        
        triggerDelays_idx = round(triggerDelays / dt);
        idx2measure       = diff(triggerDelays_idx)+1;
        triggerIdx        = [1 cumsum(idx2measure-1)];
        
        triggerOnTraj = ktraj(:, triggerIdx);
    
        if length(triggerIdx) > 1
            if all(diff(triggerDelays_idx) == diff(triggerDelays_idx(1:2)))
                StitchMode = 'Constant';
            else
                StitchMode = 'Variable';
            end
        else
            StitchMode = 'Standard';
        end
        
        plot(ktraj(1,:), ktraj(2,:), Color='#0072BD'); hold on;
        plot(triggerOnTraj(1,:), triggerOnTraj(2,:), '*', 'MarkerSize', 5, Color='#D95319');
    end
    axis equal; 
    xlabel('k_x [rad/m]', 'Color', color_label, 'fontname', font_label);
    ylabel('k_y [rad/m]', 'Color', color_label, 'fontname', font_label);
    xmax = max(abs(ktraj(1,:)));
    ymax = max(abs(ktraj(2,:)));
    
    xlim([-1.1 * xmax, 1.1 * xmax]);
    ylim([-1.1 * ymax, 1.1 * ymax]);
    
    titlestr = sprintf('%s: %d | Probe: %s | Radius: %.1f mm | Cutoff: %.2f', ...
            StitchMode, nsegs2measure, probeType, probeRadius * 1e3, signalCutoff);
    title(titlestr, 'Color', color_label);

    set(fig, 'Color', color_facecolor);  
    set(fig, 'InvertHardcopy', 'off');  
end

