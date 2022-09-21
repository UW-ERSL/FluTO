syms t a b real
n=1.5;
x=(((cos(t)))^(2/n))*a*(cos(t));
y=(((sin(t)))^(2/n))*b*(sin(t));
dx = diff(x,t);
dy = diff(y,t);
intg = sqrt( (dx)^2 + (dy)^2 );
L = 4*int(intg, t, 0, pi/2);
matlabFunction(L, 'file', 'super_ellipse_length');