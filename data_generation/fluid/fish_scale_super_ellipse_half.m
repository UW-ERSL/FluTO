function [mstr, area, length] = fish_scale_super_ellipse_half(nelx, nely, mstrsize,s)
% nelx, nely: num elems in x and y
% a, b: scalars in [0,1] that dictate size of microstr
% mstr: 0 where there is solid material and 1 where there is fluid
% area: region occupied by fluid
% length: curve length of the perimeter
a=1+2*(s-1)*(mstrsize*mstrsize*0.5);%volume
aScaled = sqrt(((1-a)/(2.737*0.5))); % to avoid parabolic scaling
bScaled = aScaled*0.5;
center = [0.5, 0.5];
mstr = ones(nelx, nely);
for rw = 1:nelx
    x = (rw/nelx)-center(1);
    for col = 1:nely
        y = (col/nely) - center(2);
        if( abs(x/aScaled)^1.5 + abs(y/bScaled)^1.5 <= 1 )
            mstr(col, rw) = 0;
        end
    end
end

area = mean(mstr , 'all');
length = 4*integral(@(t)sqrt(aScaled.^2.*cos(t).^(8.0./3.0).*sin(t).^2.*(4.9e+1./9.0)+bScaled.^2.*cos(t).^2.*sin(t).^(8.0./3.0).*(4.9e+1./9.0)),0.0,pi./2.0);


% length = pi*(  3*(aScaled+bScaled) - sqrt( (3*aScaled+bScaled)*(3*bScaled+aScaled) )  );
end