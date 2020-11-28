%%
%forward stroke
%parameters
clear all; close all; % clean workspace before beginning

syms t y

U = 1.57 * 10^-2*3*10^6;   %rad/s
w = 10*pi; %rad/s
mu = 1*10^(-6);  %Pa-s in g/micro*s 10^-2 to 10^3
p = 10^-12; %grams/micro^3
a = 1/sqrt(2)*sqrt(w*p/mu);
H = 50;    %micro
V_0 = 2.5*10^6; %microm/s

%constants
C2 = @(t) (U*cos(a*H+w*t-pi/2)*exp(a*H)-V_0)/(cos(a*H+w*t-pi/2)*exp(a*H)-cos(a*H-w*t+pi/2)*exp(-a*H));

C1 = @(t) (V_0-C2(t)*exp(-a*H)*cos(a*H-w*t+pi/2))/(cos(a*H+w*t-pi/2)*exp(a*H));

%equation
fun = @(y,t) C1(t)*exp(a*y)*cos(a*y+w*t-pi/2)+C2(t)*exp(-a*y)*cos(a*y-w*t+pi/2);

%%
%plot
time = [0:0.001:0.1];
space = [0:0.5:H];
u = zeros(101);

    for t = 1:101
        for y = 1:101
            v(t,y) = fun((y-1)/2,(t-1)/1000);
        end
    end
%%
[Y,T] = meshgrid(space,time);
figure(1)
surf(Y,T, real(v))
title('Velocity vs Y vs Time');
xlabel('Y')
ylabel('Time')
zlabel('Velocity')
colorbar
hold on;

%%
%recovery stroke
syms t y

U = 7.9 * 10^-3*3*10^6;   %rad/s changed
w = 10*pi; %rad/s
mu = 1*10^(-6);  %Pa-s in kg/micro*s 10^-2 to 10^3
p = 10^-12; %grams/micro^3
a = 1/sqrt(2)*sqrt(w*p/mu);
H = 50;    %micro
V_0 = 2.5*10^6; %microm/s

%constants
C2 = @(t) (U*cos(a*H+w*t-pi/2)*exp(a*H)-V_0)/(cos(a*H+w*t-pi/2)*exp(a*H)-cos(a*H-w*t+pi/2)*exp(-a*H));

C1 = @(t) (V_0-C2(t)*exp(-a*H)*cos(a*H-w*t+pi/2))/(cos(a*H+w*t-pi/2)*exp(a*H));

%equation
fun = @(y,t) C1(t)*exp(a*y)*cos(a*y+w*t-pi/2)+C2(t)*exp(-a*y)*cos(a*y-w*t+pi/2);

%%
%plot
time = [0.1:0.001:0.2]; %changed
space = [0:0.5:H];
u = zeros(101);

    for t = 1:101
        for y = 1:101
            u(t,y) = fun((y-1)/2,(t-1)/1000+0.1);
        end
    end
%%
ru = real(u);
[Y,T] = meshgrid(space,time);
surf(Y,T, real(u))
title('Velocity vs Y vs Time in Healthy Bronchioles');
xlabel('Y (micrometers)')
ylabel('Time s')
zlabel('Velocity (micrometers/s)')
