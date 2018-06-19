function [IndOutOfEllipse] = OutOfEllipse(Comp1, Comp2, Comp1bis, Comp2bis)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
el = @(x,y,x0, y0, a, b, alpha) (((x-x0)*cos(-alpha) + (y-y0)*sin(-alpha)).^2/(a^2))...
    + (((x-x0)*sin(-alpha) - (y-y0)*cos(-alpha)).^2/(b^2)); 

axis_handle = gca;
ellipse_t = fit_ellipse( Comp1, Comp2, [] )

 % What are the points outside the ellipse for the seizure and pre-seizure
% data 

el = @(x,y,x0, y0, a, b, alpha) (((x-x0)*cos(-alpha) + (y-y0)*sin(-alpha)).^2/(a^2))...
    + (((x-x0)*sin(-alpha) - (y-y0)*cos(-alpha)).^2/(b^2)); 

coeff = 1; 
step = 0.1;
s = length(Comp1); 

while (s > 0.001*length(Comp1))
ind = (el(Comp1, Comp2, ellipse_t.X0, ...
    ellipse_t.Y0, coeff*ellipse_t.a, coeff*ellipse_t.b, ellipse_t.phi )> 1); 
coeff = coeff + step
s = sum(ind)
end

ellipse_t.a = coeff*ellipse_t.a; 
ellipse_t.b= coeff*ellipse_t.b; 

% rotation matrix to rotate the axes with respect to an angle phi

R = [ cos(ellipse_t.phi) sin(ellipse_t.phi); -sin(ellipse_t.phi) cos(ellipse_t.phi) ];

% the axes
ver_line        = [ [ellipse_t.X0 ellipse_t.X0]; ellipse_t.Y0+ellipse_t.b*[-1 1] ];
horz_line       = [ ellipse_t.X0+ellipse_t.a*[-1 1]; [ellipse_t.Y0 ellipse_t.Y0] ];
new_ver_line    = R*ver_line;
new_horz_line   = R*horz_line;

% the ellipse
theta_r         = linspace(0,2*pi);
ellipse_x_r     = ellipse_t.X0 + ellipse_t.a*cos( theta_r );
ellipse_y_r     = ellipse_t.Y0 + ellipse_t.b*sin( theta_r );
rotated_ellipse = R * [ellipse_x_r;ellipse_y_r];

IndOutOfEllipse = (el(Comp1bis, Comp2bis, ellipse_t.X0, ...
    ellipse_t.Y0, ellipse_t.a, ellipse_t.b, ellipse_t.phi )> 1); 
end

