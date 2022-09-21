syms t a real
n=4;
x=(((cos(t)))^(2/n))*a*(cos(t));
y=(((sin(t)))^(2/n))*a*(sin(t));
dx = diff(x,t);
dy = diff(y,t);
intg = sqrt( (dx)^2 + (dy)^2 );
L = 4*int(intg, t, 0, pi/2);
matlabFunction(L, 'file', 'squircle_length');