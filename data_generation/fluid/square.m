function [mstr, area, length] = square(nelx, nely, mstrsize, s)
% nelx, nely: num elems in x and y
% a, b: scalars in [0,1] that dictate size of microstr
% mstr: 0 where there is solid material and 1 where there is fluid
% area: region occupied by fluid
% length: curve length of the perimeter
    a=1+(s-1)*(mstrsize*mstrsize);
    aScaled = sqrt(((1-a)));
    bScaled=aScaled;
    mstr = ones(nelx, nely);
    center = [0.5, 0.5]; 
    for rw = 1:nelx
        x = (rw/nelx) - center(1);
        for col = 1:nely
            y = (col/nely) - center(2);
            if((abs(x) <= 0.5*aScaled) && (abs(y) <= 0.5*bScaled))
                mstr(col, rw) = 0;
            end
        end
    end
    area = mean(mstr, 'all');
    length = 2*(aScaled+bScaled);
end