function dwdr = wDer(r,kh,type)
%WDER Computes dw/dr of the chosen kernel function

if type == 5 %Wendland 5th order
    dwdr = -320*r.*(kh-r).^3/kh^5/(16/7*pi*kh^2).*(r<=kh);
    
elseif type == 3 %Cubic b-spline
    h = kh/2;
    dwdr = ( -0.75*(2-r/h).^2.*(r/h<=2) +3*(1-r/h).^2.*(r/h<=1) )*10/(7*pi*h^2);
    
elseif type == 4 %Quartic b-spline
    h = kh/2.5;
    dwdr = ( -4*(2.5-r/h).^3.*(r/h<=2.5) +20*(1.5-r/h).^3.*(r/h<=1.5) ...
        -40*(0.5-r/h).^3.*(r/h<=0.5) )*96/(1199*pi*h^2);
end

end

