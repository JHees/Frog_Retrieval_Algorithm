function [width, crossPosition] = PulseCrossWidth(x, y, level)
    % PulseCrossWidth. This function calculates the width and cross positions at a specified level for a given pulse.
    % 此功能用于计算给定脉冲在指定电平上的宽度和交叉位置。
    % The function determines where the pulse crosses a specified level and computes the width of the pulse at that level.
    % Usage:
    %     [width, crossPosition] = PulseCrossWidth(x, y, level)
    %
    % Input:
    %     x: x-coordinates of the pulse data
    %     y: y-coordinates of the pulse data
    %     level: The level at which the width and cross positions are to be calculated
    %
    % Output:
    %     width: The width of the pulse at the specified level (distance between first rising and last falling edge)
    %     crossPosition: The positions on the x-axis where the pulse crosses the specified level
    %
    % Note:
    %     The function includes a check to ensure continuity in the pulse data and handles cases where the pulse does not cross the specified level.
    %     Calculations use linear interpolation between points for accuracy.
    %
    % AUTHOR: Huang Xingzhao, June 30, 2024

    x = x(:)';
    y = y(:)';

    x = [x, x(end) + x(2) - x(1)];
    y = [y, y(1)];

    ind = (y - level) > 0;
    ind = diff(ind);

    rising = find(ind == 1);
    falling = find(ind == -1);
    if isempty(rising) || isempty(falling)
        width = [];
        crossPosition = [];
        return;
    end
    if rising(1) > falling(1)
        rising = circshift(rising, 1);
    end

    region = [rising', falling'];
    y_interp = (level - y(region)) ./ (y(region + 1) - y(region)); % It wont be divide to 0

    crossPosition = x(region) + y_interp .* (x(region + 1) - x(region));
    if iscolumn(crossPosition)
        crossPosition = crossPosition';
    end
    width = mod(crossPosition(:, 2) - crossPosition(:, 1), x(end) - x(1));

end
