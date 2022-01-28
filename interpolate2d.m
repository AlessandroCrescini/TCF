function result = interpolate2d(f,x,y,z,x_point,y_point,z_point)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the value of f(x_point, y_point, z_point)) by  %
% linear approximation on a grid of  x  and  y, weighted %
% according to 2-dimensional distance                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% x1 is the index such that: x(x1) < x_point < x(x1+1)
% (the same holds for y). All the possible cases are included.

if x_point < x(1)
    x1 = 1;
    omega_x = 1;
elseif x_point > x(end)
    x1 = length(x)-1;
    omega_x =  0;
else 
    [temp, x1] = max(x.*(x<x_point));
    omega_x = (x(x1+1)-x_point) / (x(x1+1)-x(x1));
end


if y_point < y(1)
    y1 = 1;
    omega_y = 1;
elseif y_point > y(end)
    y1 = length(y)-1;
    omega_y =  0;
else 
    [temp, y1] = max(y.*(y<y_point));
    omega_y = (y(y1+1)-y_point) / (y(y1+1)-y(y1));
end

z_pos = z==z_point;

result = omega_x * omega_y * f(x1,z_pos,y1) ...
    + omega_x * (1-omega_y) * f(x1,z_pos,y1+1) ...
    + (1-omega_x) * omega_y * f(x1+1,z_pos,y1) ...
    + (1-omega_x)*(1-omega_y)* f(x1+1,z_pos,y1+1) ;
end


