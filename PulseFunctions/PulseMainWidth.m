function [width, region] = PulseMainWidth(x, y, level, centerX)
    % PulseMainWidth. This function calculates the width of the main peak of a pulse at a specified level.
    % 此函数用于计算level上脉冲主峰的宽度。
    % It finds the width and region of the pulse's main peak relative to a given level, using the peak's maximum as the center if not explicitly provided.
    % 
    % Usage:
    %     [width, region] = PulseMainWidth(x, y, level, centerX)
    %
    % Input:
    %     x: x-coordinates of the pulse data
    %     y: y-coordinates of the pulse data, normalized to the peak's maximum
    %     level: The level at which the main peak's width is calculated
    %     centerX: Optional, specifies the center around which to measure the peak's width (defaults to the x-value at maximum y)
    %
    % Output:
    %     width: The width of the pulse's main peak at the specified level
    %     region: The x-coordinates defining this main peak width
    %
    % Note:
    %     If the specified level does not intersect the pulse, or if the peak's calculation is ambiguous, the function returns empty arrays.
    %     Handles multiple intersections around the center, focusing on the primary peak.
    %
    % AUTHOR: Huang Xingzhao, June 30, 2024

    if nargin == 3
        [~, MaxYIndex] = max(y);
        centerX = x(MaxYIndex);
    end
    x = x(:);
    y = y(:) / max(y);
    [w, crossPosition] = PulseCrossWidth(x, y, level);
    if isempty(w)
        width = [];
        region = [];
        return
    end
    if crossPosition(1, 1) > crossPosition(1, 2)
        if (centerX > crossPosition(1, 1) || centerX < crossPosition(1, 2))
            region = crossPosition(1, :);
            width = w(1);
            return
        else
            crossPosition(1, :) = 0;
        end
    end
    cross_index_of_centerX = xor((crossPosition(:, 1) - centerX) >= 0, (crossPosition(:, 2) - centerX) >= 0);
    region = crossPosition(cross_index_of_centerX, :);
    width = w(cross_index_of_centerX);

end
