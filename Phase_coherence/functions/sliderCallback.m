% Nested callback function
function sliderCallback(ch)
    cla; hold on;
    for b = 1:nBands
        phase_vals = nan(nLocs,1);
        for iloc = 1:nLocs
            phase_vals(iloc) = phase_avg_all{iloc}(ch,b) + (b-1)*offset;
        end
        plot(1:nLocs, phase_vals, '-o', 'Color', colors(b,:), 'LineWidth', 2);
    end
    xlabel('Location Index'); ylabel('Phase + offset (rad)');
    xticks(1:nLocs); xticklabels(targ_loc);
    title(sprintf('Channel %d', ch));
    grid on;
    legend(legendLabels,'Location','bestoutside');
end