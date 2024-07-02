function [P, P_sp, f0] = removeFirstOrderPhase(varargin)
% n: 输出长度，通过fft补零增加P, P_sp分辨率
% Usage:
%     [P] = removeFirstOrderPhase(P)
%     [P] = removeFirstOrderPhase(P,n)
%     [P] = removeFirstOrderPhase(T,P)
%     [P] = removeFirstOrderPhase(T,P,n)
%     [P, P_sp] = removeFirstOrderPhase(__)
%     [P, P_sp,f0] = removeFirstOrderPhase(T,P)
%     [P, P_sp,f0] = removeFirstOrderPhase(T,P,n)

switch nargin
    case 1
        P = varargin{1};
        N = numel(P);
        T = (1:N)';
        F = T;
    case 2
        if isscalar(varargin{2})
            P = varargin{1};
            n = varargin{2};
            N = max(n, numel(P));
            T = (1:N)';
            F = T;
        else
            T = varargin{1};
            P = varargin{2};
            N = numel(P);
            F = FTconvert(T);

        end
    case 3
        T = varargin{1};
        P = varargin{2};
        n = varargin{3};
        N = max(n, numel(P));
        F = FTconvert(T);
end
level = 0.5;
level_sp = 0.5;
P_sp = fftshift(fft(P, N));
if nargout == 3
    P = ifft(ifftshift(P_sp), N);
    [~, regionP] = PulseMainWidth(T, abs(P).^2, level);
    [P, f0] = removeLinearPhase(T, P, regionP);
    P_sp = fftshift(fft(P, N));

end
[~, regionP_sp] = PulseMainWidth(F, abs(P_sp).^2, level_sp);
P_sp = removeLinearPhase(F, P_sp,regionP_sp,[],N/2);
P = ifft(ifftshift(P_sp), N);
if nargout >= 2 
    [~,maxPInd]=max(abs(P).^2);
    if maxPInd<N/4||maxPInd>N*3/4
        P = fftshift(P);
    end
end
end

function [y, phi1] = removeLinearPhase(x, y, region, level, indCenter)
N = numel(y);
if isempty(region)
    ind = abs(y) > level * max(abs(y));
else
    ind = x > region(1) & x < region(2);
end
if nargin <5 ||isempty(indCenter)
    if isempty(region)
        [~, indCenter] =max(abs(y)); 
    else
        [~, indCenter]=min(abs(x-(region(1)+region(2))/2)); 
    end
end
ang = unwrap(angle(y));
fit_p = polyfit(x(ind), ang(ind), 1);
phi1 = 0;
while abs(fit_p(1)) >= 1 / N
    phi1 = phi1 + fit_p(1);
    ang = ang - fit_p(1) .* x - fit_p(2);
    ang = unwrap(ang);
    fit_p = polyfit(x(ind), ang(ind), 1);
end
% fit_p = polyfit(x(indCenter-1:indCenter+1),ang(indCenter-1:indCenter+1),1);
% ang = ang - fit_p(1) .* x - fit_p(2);
y = abs(y) .* exp(1i * ang);
phi1 = phi1 / 2 / pi;
end
