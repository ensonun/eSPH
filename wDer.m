function dwdr = wDer(r,h,type)
%WDER Computes dw/dr of the chosen kernel function
%   Detailed explanation goes here

if type == 1
    dwdr = -12*r.*(h-r).^2/(0.2*pi*h^6).*(r<=h);
    
elseif type == 2 %Wendland
    dwdr = -320*r.*(h-r).^3/h^5/(16/7*pi*h^2).*(r<=h);
    
elseif type == 3 %Cubic spline
    dwdr = ( (-8*r/h^2 + 12*r.^2/h^3).*(r<=h/2) ...
        + (-(2-2*r/h).^2/h).*(r>h/2).*(r<=h) )/(7/60*pi*h^2);
end

end

