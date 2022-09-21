function [mstr, area, length] = squircle(nelx, nely, mstrsize,s)
% nelx, nely: num elems in x and y
% a, b: scalars in [0,1] that dictate size of microstr
% mstr: 0 where there is solid material and 1 where there is fluid
% area: region occupied by fluid
% length: curve length of the perimeter
a=1+(s-1)*(mstrsize*mstrsize);
aScaled = sqrt((1-a)/3.71); % to avoid parabolic scaling
center = [0.5, 0.5];
mstr = ones(nelx, nely);
for rw = 1:nelx
    x = (rw/nelx)-center(1);
    for col = 1:nely
        y = (col/nely) - center(2);
        if( abs(x/aScaled)^4 + abs(y/aScaled)^4 <= 1 )
            mstr(col, rw) = 0;
        end
    end
end

area = mean(mstr , 'all');
length = integral(@(t)sqrt(cos(t).*sin(t).^2.*(aScaled.^2).*(9.0./4.0)+cos(t).^2.*sin(t).*(aScaled.^2).*(9.0./4.0)),0.0,pi./2.0).*4.0;
end