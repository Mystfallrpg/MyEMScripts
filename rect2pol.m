function [r,a] = rect2pol(x)

r_temp = abs(x);
a_temp = rad2deg(angle(x));

r = simplify(r_temp);
a = simplify(a_temp);

end