function width = PulseWidth(x, y, method)
if nargin == 2
    method = 'fwhm';
end

if isrow(y)
    y = y(:);
end
if isrow(x)
    x = x(:);
    x = repmat(x, [1, size(y, 2)]);
end
[x, ind] = sort(x, 1);
y = y(ind);
y = y - min(y);
y = y / max(y, [], 1);

switch method
    case 'fwhm'
        width = zeros(1, size(y, 2));
        for i = 1:size(y, 2)
            width(i) = full_width_half_maximum(x(:, i)', y(:, i)');
        end
    case 'hw1e'
        width = zeros(1, size(y, 2));
        for i = 1:size(y, 2)
            width(i) = half_width_1e(x(:, i)', y(:, i)');
        end
    case 'rms'
        width = root_mean_squared_pulse_width(x, y);
    case 'equ'
        width = equivalent_pulse_width(x, y);
    otherwise
        fprintf("unknow method.")
        midcross
end
end

function width = full_width_half_maximum(x, y)
width = PulseFullWidth(x, y, 0.5);
end

function width = half_width_1e(x, y)
width = PulseFullWidth(x, y, exp(-1)) ./ 2;
end

function rms = root_mean_squared_pulse_width(x, y)
y = y ./ trapz(x,y,1);
rms = sqrt(trapz(x, x.^2 .* y, 1) - trapz(x, x .* y, 1).^2);
end

function width = equivalent_pulse_width(x, y)
width = trapz(x, y, 1) ./ max(y, [], 1);
end
