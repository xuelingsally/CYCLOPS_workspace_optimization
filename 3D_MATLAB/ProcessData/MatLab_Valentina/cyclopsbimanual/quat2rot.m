% Conversion from quaternion to rotation matrix

function [Rot]=quat2rot(Q)

s = Q(1);
x = Q(2);
y = Q(3);
z = Q(4);

Rot = [   1-2*(y^2+z^2)   2*(x*y-s*z) 2*(x*z+s*y)
        2*(x*y+s*z) 1-2*(x^2+z^2)   2*(y*z-s*x)
        2*(x*z-s*y) 2*(y*z+s*x) 1-2*(x^2+y^2)   ];