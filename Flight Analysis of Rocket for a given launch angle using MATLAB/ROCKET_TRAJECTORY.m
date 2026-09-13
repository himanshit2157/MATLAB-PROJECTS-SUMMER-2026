clear;
clc;
close all;

%parameters all in SI units
g=9.81; %gravity
rho0=1.225; %air density at sea level
h=10000; %scale height 
m0=80; % initial mass
m1=30; %final mass
burntime=10; %burn time
cd=0.5; %drag coefficient assumed
farea=0.02; %frontal area
thrust=3000; %thrust
lang = input('Enter launch angle in degrees: ');
lang = deg2rad(lang); %converting deg to radians
simtimespan=[0,120]; % simulation time span of 120 seconds

%initial conditions of position and velocities x,y,vx,vy of the rocket 
ic=[0;0;0;0];

%prop mass flow rate constant
mp=(m0-m1)/burntime; 

%stop when ground hit setting, define touchdown fun
opts = odeset('Events', @touchdown);

%solving eq of motion using ode45 and anonymous fun because using parameters
[t,st]=ode45(@(t,st) rocketDY(t,st,m0,m1,thrust,burntime,cd,farea,g,rho0,h,lang,mp),simtimespan,ic,opts);

%extract positions and velocities
x=st(:,1); %range
y=st(:,2); %altitude
vx=st(:,3); %vel in x dir
vy=st(:,4); %vel in y dir

%total velocity, density, dynamic pressure
v=sqrt(vx.^2+vy.^2); 
rho=rho0.*exp(-y/h); 
q=0.5.*rho.*v.^2 ;

% max values
[maxAl, alindex] = max(y); %max altitude
[maxVel, veindex] = max(v); %max velocity
maxDP = max(q); %max dynamic pressure
range = x(end); % final range

%display
fprintf('SIMULATION RESULTS \n')
fprintf('launch angle: %.2f \n', rad2deg(lang));
fprintf('Maximum altitude : %.2f m \n', maxAl);
fprintf('Time at which Max altitude was reached : %.2f s \n', t(alindex));
fprintf('Maximum velocity : %.3f m/s \n', maxVel);
fprintf('Time at which Max velocity was reached : %.2f s \n', t(veindex));
fprintf('Final range: %.2f m \n', range);
fprintf('Maximum dynamic pressure: %.3f Pa \n', maxDP);

% Plotting the trajectory, dynamic pressure, altitude and velocity
figure;
subplot(2,1,1) %horizontal distance as ascending, will be minuscule id launch angle is 90deg
plot(x, y,'c',LineWidth=1.5); 
xlabel('Horizontal Distance m');
ylabel('Altitude m');
title('Rocket Trajectory [2D]');
grid on;

subplot(2,1,2);
plot(t,q,'g',LineWidth=1.5);
xlabel('Time s');
ylabel('Dynamic Pressure Pa');
title('Dynamic Pressure vs Time');
grid on;

figure;
subplot(2,1,1);
plot(t,y,'b',LineWidth=1.5);
xlabel('Time s');
ylabel('Altitude m');
title('Altitude vs Time');
grid on;
subplot(2,1,2);
plot(t,v,'w',LineWidth=1.5);
xlabel('Time s');
ylabel('Velocity m/s');
title('Velocity vs Time');
grid on;

%this part is independent of the given launch angle
%with different launch angles from 30 to 90deg with intervals of 5
langs=30:5:90;
%making placeholder arrays
maxAls=zeros(size(langs));
ranges=zeros(size(langs));
for i=1:length(langs)
    rang=deg2rad(langs(i));
    [ti,sti] = ode45(@(t,st) rocketDY(t,st,m0,m1,thrust,burntime,cd,farea,g,rho0,h,rang,mp), simtimespan,ic,opts);
    %position
    xi = sti(:,1); %ranges
    yi = sti(:,2); %altitudes
    %maximum altitude
    maxAls(i) = max(yi);
    %range
    ranges(i) = xi(end); 
end

% Display results for all the different launch angles
figure;
subplot(2,1,1);
plot(langs, maxAls,'y',LineWidth=1.5);
xlabel('Launch Angle deg');
ylabel('Max Altitude m');
title('Maximum Altitude vs Launch Angle');
grid on;

subplot(2,1,2);
plot(langs, ranges,'m',LineWidth=1.5);
xlabel('Launch Angle deg');
ylabel('Range m');
title('Range vs Launch Angle');
grid on;


%function definitions for rocketDY and touchdown

function dst = rocketDY(t, st, m0, m1, thrust, burntime, cd, farea, g, rho0, h, lang, mp)
    %extract states and put them in st matrix
    x=st(1);
    y=st(2);
    vx=st(3);
    vy=st(4);
    
    %current mass
    if t < burntime %if time is less than burning time i.e. motor is burning
        m = m0 - mp * t;
        currentthrust=thrust;
    else
        m = m1;% after burn is over, mass is the final mass
        currentthrust=0;
    end
    
    %drag force
    v = sqrt(vx^2 + vy^2);
    rho = rho0 * exp(-y / h);
    drag = 0.5 * rho * v^2 * cd * farea;

    %Drag in components
    if v > 0
        Dragx = drag * (vx / v);
        Dragy = drag * (vy / v);
    else
        Dragx = 0;
        Dragy = 0;
    end
    
    %thrust into components
    if t < burntime
        Thrustx = currentthrust * cos(lang);
        Thrusty = currentthrust * sin(lang);
    else
        Thrustx = 0;
        Thrusty = 0;
    end
    
    %eq of motion
    ax = (Thrustx - Dragx) / m;
    ay = (Thrusty - Dragy) / m - g;
    
    %derivatives
    dst = [vx; vy; ax; ay];
end

function [value, terminal, direction] = touchdown(t,st)
    value = st(2);  %watch the altitude and stop when altitude y is 0m
    terminal = 1; %true (yes stop)    
    direction = -1; %only trigger when descending down past 0 and -ve
end

