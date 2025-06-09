function [fig] = plot_kspace(ktraj, ktraj_adc)
    color_facecolor = "#1F1F1F";
    color_label     = "#CCCCCC";
    font_label      = "Arial";

    fig_width       = 800;
    fig_height      = 800;
    position = [(1920-fig_width)/2, (1080-fig_height)/2, fig_width, fig_height];

    figname = '2D k-space';
    fig = figure('Name', figname, 'Position', position, 'Color', color_facecolor);
    subplot(1,1,1, 'color', color_facecolor, 'xcolor', color_label, 'ycolor', color_label, 'fontname', font_label), hold on;

    color_traj      = '#1d6996';
    color_adc       = '#a85144';

    plot(ktraj(1,:),ktraj(2,:), 'color', color_traj); % a 2D plot
    hold on;plot(ktraj_adc(1,:),ktraj_adc(2,:), '.', 'MarkerSize', 5, 'MarkerEdgeColor', color_adc); 
    axis equal;

    ax = gca;
    ax.Position = [0.08 0.08 0.9 0.9];      
    ax.LooseInset = [0 0 0 0]; 
    set(fig, 'Color', color_facecolor);  
    set(fig, 'InvertHardcopy', 'off');  
end

