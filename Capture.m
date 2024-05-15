clear all
close all

% Relevant Constants
G = 6.6743e-11; 
M = 2e30; % star mass
mu = G * M;
m1 = 3e28; % planet 1 mass
m2 = 5e20; % planet 2 mass
earth_year = 365.25*24*60*60; 
R = 0.5*(sqrt(2*G*M)*earth_year/pi)^(2/3); % earth radius


%Time
tmax = 5 * earth_year; 
dt = 4e-5 * earth_year; 
clockmax = floor(tmax / dt);

%Planet 1 IC
theta = pi; 
w_mag = sqrt(mu / R);
z = R * [cos(theta), sin(theta)]; % position vector
w = w_mag * [sin(theta), -cos(theta)]; % velocity vector

T = pi * sqrt(R^3 / mu); % half period 


%Planet 2 IC
q = .01 * R; % <----- FREE, GAP -------------------
ecc = 1.5; % <----- FREE, ECCENTRICITY (>1)-------------------
[x, v] = getIC(R, q, ecc, T);




% Arrays to store trajectory:
tsave = zeros(1,clockmax);
x1save = zeros(1,clockmax);
y1save = zeros(1,clockmax);
x2save = zeros(1,clockmax);
y2save = zeros(1,clockmax);


% Initialization for animation:
plot(0,0,'k*','linewidth',6) % Sun at center
hold on

% Planet Handles
hp1 = plot(z(1),z(2),'bo','linewidth',3); %blue 
hp2 = plot(x(1),x(2),'mo','linewidth',3); %magenta

% Trajectory Handles
ht1 = plot(z(1),z(2),'b','linewidth',1);
ht2 = plot(x(1),x(2),'m','linewidth',1); 


% Axes
ax = 2 * R;
axis equal 
axis([-ax,ax,-ax,ax])
axis manual 

% Video
title("Eccentricity: ", ecc)
writerObj1 = VideoWriter('E5.avi');
open(writerObj1);
currentFrame = getframe(gcf);
writeVideo(writerObj1, currentFrame);

% Main Loop
for clock = 1:clockmax
    t = clock * dt;

    % Intermediate Values
    D1 = z / norm(z)^3;
    D2 = x / norm(x)^3;
    D3 = (x - z) / norm(x - z)^3;
    
    % Update velocities
    w = w + dt * G*(m2 * D3 - M * D1);
    v = v + dt * G*(-m1 * D3 - M * D2);

    % Update Positions (using new velocities)
    z = z + dt * w;
    x = x + dt * v;
    
    % Save Data
    tsave(clock) = t;
    x1save(clock) = z(1);
    y1save(clock) = z(2);
    x2save(clock) = x(1);
    y2save(clock) = x(2);

  
    % Update Handles and Plot

    if (mod(clock,100) == 0)
        hp1.XData = z(1);
        hp1.YData = z(2);
        ht1.XData = x1save(1:clock);
        ht1.YData = y1save(1:clock);
    
        hp2.XData = x(1);
        hp2.YData = x(2);
        ht2.XData = x2save(1:clock);
        ht2.YData = y2save(1:clock);
        drawnow 
        currentFrame = getframe(gcf);
        writeVideo(writerObj1, currentFrame);
    end

end

close(writerObj1);



