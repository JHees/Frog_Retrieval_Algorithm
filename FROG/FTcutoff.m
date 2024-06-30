function [in] = FTcutoff(in, max)
% 将输入的F/T在max处截断
diffs = in(2) - in(1);
in = in(in >- max & in < max - diffs);

end
