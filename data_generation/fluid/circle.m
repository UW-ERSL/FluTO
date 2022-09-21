function [mstr, area, length] = circle(nelx, nely, mstrsize, s)
% nelx, nely: num elems in x and y
% a, b: scalars in [0,1] that dictate size of microstr
% mstr: 0 where there is solid material and 1 where there is fluid
% area: region occupied by fluid
% length: curve length of the perimeter
as_ratio = 1; %aspect ratio
a=1+(s-1)*(mstrsize*mstrsize);%volume
vfCoeff = (1-s)/(0.25*as_ratio);
aScaled = sqrt(((1-a)/(vfCoeff*as_ratio))); % to avoid parabolic scaling
center = [0.5, 0.5];
mstr = ones(nelx, nely);
for rw = 1:nelx
    x = (rw/nelx)-center(1);
    for col = 1:nely
        y = (col/nely) - center(2);
        if( (x/aScaled)^2 + (y/aScaled)^2 <= 1 )
            mstr(col, rw) = 0;
        end
    end
end

area = mean(mstr , 'all');
length = 2*pi*aScaled;
end