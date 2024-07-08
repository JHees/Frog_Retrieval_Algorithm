function [out, in] = FTconvert(in, N, k)
% N < length(in)时抽样, 否则向外以相同间隔延展
diffs = in(2) - in(1);
switch nargin
    case 1
        N = length(in);
    case 2
        if length(in) > N
            diffs = diffs * floor(length(in) / N);
        end
    case 3
        diffs = diffs * k;
end

if isrow(in)
    out = (- N / 2:N / 2 - 1) / diffs / N;
    in = (- N / 2:N / 2 - 1) * diffs;
else
    out = (- N / 2:N / 2 - 1)' / diffs / N;
    in = (- N / 2:N / 2 - 1)' * diffs;
end
end
