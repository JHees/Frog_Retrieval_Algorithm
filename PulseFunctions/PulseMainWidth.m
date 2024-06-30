function [width, region] = PulseMainWidth(x, y, level, centerX)
if nargin == 3
    [~, MaxYIndex] = max(y);
    centerX = x(MaxYIndex);
end
x = x(:);
y = y(:)/max(y);
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
