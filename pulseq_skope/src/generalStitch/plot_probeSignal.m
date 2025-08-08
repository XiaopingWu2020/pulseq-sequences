function [fig] = plot_probeSignal(probeSignals, t, probeType, probeRadius, signalCutoff)
%PLOT_STITCH 
%   
    color_facecolor = "#1F1F1F";
    color_label     = "#CCCCCC";
    font_label      = "Arial";
    fig_width       = 30;
    fig_height      = 5;
    position = [10, 15, fig_width, fig_height];
    figname = 'Stitch Pattern';
    fig = figure('Name', figname, 'Units', 'centimeters', 'Position', position, 'Color', color_facecolor);
    subplot(1,1,1, 'color', color_facecolor, 'xcolor', color_label, 'ycolor', color_label, 'fontname', font_label), hold on;

    nSeg = size(probeSignals, 2);

    hold on;
    for i = 1:nSeg
        ss = probeSignals{i};
        ss(ss < signalCutoff) = 0;
        % plot(t*1e3, ss)
        nonzeroMask = ss > 0;
        % plot(t(nonzeroMask)*1e3, ss(nonzeroMask))
        scatter(t(nonzeroMask)*1e3, ss(nonzeroMask)*100, 1, 'filled')
    end

    xlabel('Time [ms]', 'Color', color_label, 'fontname', font_label);
    ylabel('Probe signal [%]', 'Color', color_label, 'fontname', font_label);
    
    xlim([-1, t(end)*1e3 + 1]);
    ylim([0, 100]);
    
    titlestr = sprintf('nSeg: %d | Probe: %s | Radius: %.1f mm | Cutoff: %.0f%%', ...
            nSeg, probeType, probeRadius * 1e3, signalCutoff*100);
    title(titlestr, 'Color', color_label);

    set(fig, 'Color', color_facecolor);  
    set(fig, 'InvertHardcopy', 'off');  
    set(gca, 'FontName', 'Arial', 'FontSize', 8);        
    set(get(gca, 'XLabel'), 'FontName', 'Arial', 'FontSize', 9); 
    set(get(gca, 'YLabel'), 'FontName', 'Arial', 'FontSize', 9); 
end

