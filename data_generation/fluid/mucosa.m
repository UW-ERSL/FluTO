function [mstr, area, length] = mucosa(nelx, nely, mstrsize, s, nvilli)
% nelx, nely: num elems in x and y
% Ri : scalar in [0.1,0.95] %Inner Radius of mucosa
% Ro : scalar in [0.1,0.95] %Outer Radius of mucosa
% nvilli: %No. of microvilli projections
% eps:  non zero offset
% mstr: 0 where there is solid material and 1 where there is fluid
% area: region occupied by fluid
% length: curve length of the perimeter
as_ratio = 0.9; %aspect ratio
a=1+(s-1)*(mstrsize*mstrsize);%volume
vfCoeff = (1-s)/(0.25*as_ratio);
Ro = sqrt(((1-a)/(vfCoeff*as_ratio))); 
%Ro = 0.5; 
Ri= as_ratio*Ro;

mstr = ones(nelx, nely);
center = [0.5, 0.5]; 
    for rw = 1:nelx
        x = (rw/nelx) - center(1);
        for col = 1:nely
            y = (col/nely) - center(2);
            if (x==0 && y==0)
                theta=0;
            else
            theta=atan(y/x);
            end
            kx=(0.5*(Ro+Ri)+(Ro-Ri)*sin(nvilli*theta))*cos(theta);
            ky=(0.5*(Ro+Ri)+(Ro-Ri)*sin(nvilli*theta))*sin(theta);
            if((x^2+y^2)<(kx^2)+(ky^2))
                mstr(col, rw) = 0;
            end
        end
    end
area = mean(mstr , 'all');
    
length = integral(@(t)sqrt((sin(t).*(0.5*Ri+0.5*Ro-sin(nvilli.*t).*(Ri-Ro))+nvilli.*cos(t).*cos...
    (nvilli.*t).*(Ri-Ro)).^2+(cos(t).*(0.5*Ri+0.5*Ro-sin(nvilli.*t).*(Ri-Ro))-nvilli.*sin(t).*cos...
    (nvilli.*t).*(Ri-Ro)).^2),0.0,pi.*2.0);
