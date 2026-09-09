clear;
clc;
close all;

% parameters all in SI units
ro=0.1;
rp=0.03;
gl=0.5;
prho=1900;
a=5e-5; % burn rate coeff 
n=0.35; % pressure coeff
gamma= 1.20; % ratio of specific heats
R=300; % specific gas constant
tc= 3000; % combustion temp
at= 1e-4;  % throat area 
ae= 3e-3; % exit area 
pa= 101325; % ambient pressure
ipvol= pi*(ro^2-rp^2)*gl; % initial prop vol
ipmass= prho*ipvol; % initial prop mass

fprintf('initial prop mass= %.3f kg\n', ipmass);
fprintf('initial port radius= %.4f m\n', rp);
fprintf('grain length= %.3f m\n', gl);

istate= rp;
tspan=[0, 30];

options = odeset('Events', @(t, r) burn_events(t, r, ro));
[t, rport] = ode45(@(t,r) motorDY(t,r,prho,a,n,gamma,R,tc,at,gl), tspan, istate, options);

% placeholder arrays
burnrate= zeros(size(t));
burnarea= zeros(size(t));
massgen= zeros(size(t));
chamberpa= zeros(size(t));
massflownozzle= zeros(size(t));
exitvel= zeros(size(t));
thrust = zeros(size(t));
specificimp= zeros(size(t));

% ode45 results processing loop
for i=1:length(t)
    r=rport(i); 
    ab=2*pi*r*gl; % burn area
    burnarea(i)=ab;

    pc=chambpa(r,prho,a,n,gamma,R,tc,at,gl); % chamber pressure
    chamberpa(i)=pc;

    rb=a*pc^n; % burnrate
    burnrate(i)=rb;

    mgen=prho*ab*rb; % prop mass generation
    massgen(i)=mgen;

    mchoked=pc*at*sqrt(gamma/R/tc)*(2/(gamma+1))^((gamma+1)/(2*(gamma-1))); % Fixed choked mass flow 
    massflownozzle(i)=mchoked;

    expratio=ae/at; 
    me=emach(expratio,gamma); % exit mach no.
    
    te=tc/(1+(gamma-1)/2*me^2); % exit temp

    ve=me*sqrt(gamma*R*te);% exit velocity
    exitvel(i)=ve;

    pe=pc*(1+(gamma-1)/2*me^2)^(-gamma/(gamma-1)); % exit pressure

    ft=mchoked*ve+(pe-pa)*ae; % thrust
    thrust(i)=ft;

    g0=9.806;
    isp=ft/(mchoked*g0); % specific impulse
    specificimp(i) = isp; 
end

burnend=find(rport>=ro,1);
if isempty(burnend)
    burnend=length(t);
end
burntime=t(burnend);

% imp results
[maxPa,maxPaind]=max(chamberpa); 
[maxthrust,maxtind]=max(thrust);
[maxburnate,maxburnrateind]=max(burnrate);
[maxisp,maxispind]=max(specificimp);

% Display results
fprintf('Burn time: %.3f s\n', burntime);
fprintf('Max chamber pressure: %.3f Pa at index %d\n', maxPa, maxPaind);
fprintf('Max thrust: %.3f N at index %d\n', maxthrust, maxtind);
fprintf('Max burn rate: %.4f m/s at index %d\n', maxburnate, maxburnrateind); 
fprintf('Max specific impulse: %.3f s at index %d\n', maxisp, maxispind);

% Plot results
figure;
subplot(2,3,1); plot(t, chamberpa,'b'); xlabel('Time (s)'); ylabel('Chamber Pressure (Pa)'); title('Chamber Pressure vs Time');
subplot(2,3,2); plot(t, thrust,'y'); xlabel('Time (s)'); ylabel('Thrust (N)'); title('Thrust vs Time');
subplot(2,3,3); plot(t, burnrate,'g'); xlabel('Time (s)'); ylabel('Burn rate (m/s)'); title('Burn Rate vs Time');
subplot(2,3,4); plot(t, rport,'r'); xlabel('Time (s)'); ylabel('Port radius (m)'); title('Port Radius vs Time');
subplot(2,3,5); plot(t, burnarea,'c'); xlabel('Time (s)'); ylabel('Burn Area (m^2)'); title('Burn Area vs Time');
subplot(2,3,6); plot(t, specificimp,'w'); xlabel('Time (s)'); ylabel('Specific Impulse (s)'); title('Specific Impulse vs Time');

function drdt = motorDY(~, r, prho, a, n, gamma, R, tc, at, gl)
    pc = chambpa(r, prho, a, n, gamma, R, tc, at, gl);
    drdt = a * (pc)^n; 
end

function pc = chambpa(r, prho, a, n, gamma, R, tc, at, gl)
    ab = 2 * pi * r * gl; 
    c_star = sqrt(R * tc) / (sqrt(gamma) * (2 / (gamma + 1))^((gamma + 1) / (2 * (gamma - 1))));
    pc = ((a * prho * ab) / (at / c_star))^(1 / (1 - n));
end

function me = emach(expratio, gamma)  % using Newton-Raphson method
    me = 2.5; % Initial guess supersonic flow
    for iter = 1:50
        f = (1/me) * ((2/(gamma+1)) * (1 + (gamma-1)/2 * me^2))^((gamma+1)/(2*(gamma-1))) - expratio;
        df = -1/(me^2) * ((2/(gamma+1)) * (1 + (gamma-1)/2 * me^2))^((gamma+1)/(2*(gamma-1))) + ...
             (1/me) * ((gamma+1)/(2*(gamma-1))) * ((2/(gamma+1)) * (1 + (gamma-1)/2 * me^2))^(((gamma+1)/(2*(gamma-1))) - 1) * (2/(gamma+1)) * (gamma-1)*me;
        me = me - f/df;
        if abs(f) < 1e-6
            break;
        end
    end
end

function [value, isterminal, direction] = burn_events(~, r, ro)
    value = r - ro;      
    isterminal = 1;      
    direction = 1;       
end


