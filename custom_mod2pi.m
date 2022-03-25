function [y] = custom_mod2pi(x) %helper function
    y = x-2*pi*floor(0.5+x/(2*pi));
end