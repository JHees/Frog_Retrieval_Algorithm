function [width, crossPosition] = PulseCrossWidth(x, y, level)
x=x(:)';
y=y(:)';

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
