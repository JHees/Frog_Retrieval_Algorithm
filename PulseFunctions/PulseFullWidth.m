function [width, region] = PulseFullWidth(x, y, level)
    % PulseFullWidth. This function calculates the full width of a pulse at a specified level.
    % 此函数用于计算level上脉冲的全宽。
    % It determines the full width at the given level by inverting the pulse measurement and adjusting from the pulse ends.
    % Usage:
    %     width = PulseFullWidth(x, y, level)
    %
    % Input:
    %     x: x-coordinates of the pulse data
    %     y: y-coordinates of the pulse data
    %     level: The level at which the full width is calculated
    %
    % Output:
    %     width: The full width of the pulse at the specified level
    %
    % Note:
    %     The function uses the PulseCrossWidth function to find crossing widths at an inverted level, then subtracts from the total span.
    %
    % AUTHOR: Huang Xingzhao, June 30, 2024

    [w, r] = PulseCrossWidth(x, max(y) - y, max(y) - level);
    [max_w, ind] = max(w);
    width = x(end) - x(1) - max_w;
    region = [r(ind, 2), r(ind, 1)];
end
