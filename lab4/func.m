function [f] = func(i,j,k,hx,hy,hz)
x = (i-1)*hx;
y = (j-1)*hy;
z = (k-1)*hz;
f = x^2+y^2+z^2;
end

