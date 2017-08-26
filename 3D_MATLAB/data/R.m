function [ R ] = R( angle,around )
    % second argument is the coordinate around which the rotation goes,
    % x,y,z => 1,2,3 resp.
    % e.g. rotation around z coordinate: R(theta,3);
    
    switch(around)
        case 1
            R = [   1   0     0;
                    0   cos(angle) -sin(angle);
                    0   sin(angle)  cos(angle)];              
        case 2
            R = [  cos(angle)      0     sin(angle);
                    0               1    0  ;
                    -sin(angle)      0   cos(angle)];  
        case 3
            R = [  cos(angle)     -sin(angle)   0;
                    sin(angle)     cos(angle)   0;
                    0               0           1];  
    end
end

