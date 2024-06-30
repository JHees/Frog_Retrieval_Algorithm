function width = PulseFullWidth(x, y, level)
w = PulseCrossWidth(x, max(y) - y, max(y) - level);
width = x(end) - x(1) - max(w);
end
