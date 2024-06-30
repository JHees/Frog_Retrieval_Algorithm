function res = solveFrog_vanilla(I, initialGuess, iterMax, eps, varargin)
% I should be normalize as maximum equal to 1
p = inputParser;
p.addParameter("over_correction", 1);
p.parse(varargin{:});
b = p.Results.over_correction;

N = size(I, 1);

if isempty(initialGuess)
    x = (-N / 2:N / 2 - 1)';
    x_tou = N / 16;
    P = exp(- x.^2/2 ./ x_tou^2) + 0.5 * exp(- (x + N / 8).^2/2 ./ x_tou^2);
else
    P = initialGuess(:);
end

res.err = zeros(1, iterMax);
i = 1;
while i <= iterMax
    if isnumeric(b)
        [P, res.err(i)] = vanilla(I, P, b);
    else
        b = min((1.1)^(1 + i / 5),3);
        [P, res.err(i)] = vanilla(I, P, b);
    end
    if res.err(i) < eps
        break;
    end
    i = i + 1;
end

P = P ./ max(abs(P));
[res.P, res.P_sp] = removeFirstOrderPhase(P);
res.iter = i-1;
res.err = res.err(1:i-1);
end

function sig = data_constraint(I, sig_sp, b)
% if b ==1: sig_sp = sqrt(I) .* exp(1i * angle(sig_sp));
sig_sp = sig_sp .* (sqrt(I) ./ abs(sig_sp)).^b;
sig_sp(isnan(sig_sp)) = 0;
sig = ifft(ifftshift(sig_sp, 1), [], 1);
end

function [P, err] = vanilla(I, P, b)
[~, sig_sp] = TraceGenerate(P);
sig = data_constraint(I, sig_sp, b);
P = sum(sig, 2);
P =P./max(abs(P));
err = TraceError(I, P);
end
