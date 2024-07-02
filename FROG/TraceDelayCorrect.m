function [dt,Delay_time,Delay_N] = TraceDelayCorrect(I, D, F, iter)
    [~, ND] = size(I);
    dt_res = zeros(1, iter + 1);
    if isscalar(D)
        dt_res(1) = D;
        D = (-ND / 2:ND / 2 - 1)' * D;
    else
        dt_res(1) = D(2) - D(1);
    end

    for i = 1:iter
        shg = I(:, end / 2 + 1);
        shg_ft = fftshift(ifft(shg, 'symmetric'));
        p = mypeak(FTconvert(F), abs(shg_ft).^2);
        Delay_time = (abs(p(1) - p(2)) + abs(p(1) - p(3))) / 2;

        p = mypeak(autocorrelation(I));
        Delay_N = (abs(p(1) - p(2)) + abs(p(1) - p(3))) / 2;
        dt_res(i + 1) = Delay_time / Delay_N;
        if abs(dt_res(i + 1) - dt_res(i)) < 1e-8
            break;
        end
    end
    dt = dt_res(i + 1);
end

function p = mypeak(x, y)
    if nargin == 1
        y = x;
        x = 1:length(y);
    end
    [~, pos_ind, ~, prominence] = findpeaks(y, "WidthReference", "halfheight", "MinPeakHeight", 0);
    [~, ind] = sort(prominence, 'descend');
    pos_ind = pos_ind(ind);
    p = zeros(1, 3);
    for ii = 1:3
        x1 = x(pos_ind(ii) - 1);
        x2 = x(pos_ind(ii));
        x3 = x(pos_ind(ii) + 1);
        y1 = y(pos_ind(ii) - 1);
        y2 = y(pos_ind(ii));
        y3 = y(pos_ind(ii) + 1);
        A = [ ...
                 x1^2, x1, 1; ...
                 x2^2, x2, 1; ...
                 x3^2, x3, 1; ...
             ];
        b = [y1; y2; y3];
        coff = A \ b;
        p(ii) = -coff(2) / 2 / coff(1);
    end
end
