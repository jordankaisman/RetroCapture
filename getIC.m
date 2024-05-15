function [x_0,v_0] = getIC(R, q, e, T)


    % Relevant Constants
    G = 6.6743e-11; 
    M = 2e30; % star mass
    mu = G * M;
    
    
    %Planet 2 IC
    x = [R, -q]; % position vector

    a = (e * q + sqrt(q^2 + R^2)) / (e^2 - 1);
    c = a * e;
    h = (c-q) * (e^2-1) / R;
    
    speed = sqrt(mu * ( 2 / norm(x) + 1 / a ) );

    v = speed / sqrt(1+h^2) * [h, 1]; % velocity vector
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Time
    tmax = T; 
    clockmax = 8000; 
    dt = tmax/clockmax; 
    
    
    % Main Loop
    for clock = 1:clockmax
    
        % Intermediate Values
        D2 = x / norm(x)^3;
    
        % Update velocities
        v = v + dt * G*(- M * D2);
    
        % Update Positions (using new velocities)
        x = x + dt * v;
    end
    
    % Return terminal position and velocity
    x_0 = x;
    v_0 = -v;


end

