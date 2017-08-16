% Conversion from rotation matrix to quaternion

function [Q]=rot2quat(T)

rot=T(1:3,1:3);
qs = sqrt(trace(rot)+1)/2.0;
kx = rot(3,2) - rot(2,3);   % Oz - Ay
ky = rot(1,3) - rot(3,1);   % Ax - Nz
kz = rot(2,1) - rot(1,2);   % Ny - Ox

if (rot(1,1) >= rot(2,2)) && (rot(1,1) >= rot(3,3))
    kx1 = rot(1,1) - rot(2,2) - rot(3,3) + 1; % Nx - Oy - Az + 1
    ky1 = rot(2,1) + rot(1,2);          % Ny + Ox
    kz1 = rot(3,1) + rot(1,3);          % Nz + Ax
    add = (kx >= 0);
elseif (rot(2,2) >= rot(3,3))
    kx1 = rot(2,1) + rot(1,2);          % Ny + Ox
    ky1 = rot(2,2) - rot(1,1) - rot(3,3) + 1; % Oy - Nx - Az + 1
    kz1 = rot(3,2) + rot(2,3);          % Oz + Ay
    add = (ky >= 0);
else
    kx1 = rot(3,1) + rot(1,3);          % Nz + Ax
    ky1 = rot(3,2) + rot(2,3);          % Oz + Ay
    kz1 = rot(3,3) - rot(1,1) - rot(2,2) + 1; % Az - Nx - Oy + 1
    add = (kz >= 0);
end

if add
    kx = kx + kx1;
    ky = ky + ky1;
    kz = kz + kz1;
else
    kx = kx - kx1;
    ky = ky - ky1;
    kz = kz - kz1;
end

nm = norm([kx ky kz]);

if nm == 0
    Q = [1 0 0 0];
else
    s = sqrt(1 - qs^2) / nm;
    qv = s*[kx ky kz];
    
    Q = [qs qv];
end