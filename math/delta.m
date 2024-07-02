function f = delta(x, x0)
    f = ismembertol(x, x0, 1e-8);
end
