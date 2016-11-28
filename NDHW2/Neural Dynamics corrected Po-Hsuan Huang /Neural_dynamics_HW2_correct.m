% Neural Dynamics Exercise 2
% Problem 3(1)
% Po-Hsuan Huang 2014.11.11
% Simulation of single compartment model.

%% specs of modeled neuro compartment
close all;  % close all figures
clc;       % clear command area
clear;     % clear variables
% a good hobby I learned today.
% control +C kills the process

Em = 0;            % in volt
l = 100E-6 ;     % in meter
diameter = 2E-6  ;    % in meter
sRm = 1  ;         % specific memebrane resistance in ohm*sqrmeter
sRa = 1  ;        % specific axial resistivity in ohm*mete
sCm = 0.01 ;       % specific memberane capacitance in Farad/sqrmeter
f =10;              % ingective current frequency
%% calculating the parameter in analogic circuit.
Rm = sRm/(pi*diameter*l);     % membrane resistance in ohm.
Cm =  sCm*pi*diameter*l;                       % capacitance in Farad.
Ie = -50e-12;                                % injective constant current in pA.

%% backward Euler Method
dt= 0.0001;   % zie of timestep in second.
time=0:dt:1;
v= zeros(1,length(time));
v(1,1) = Em; % construct a matrix to store datapoints.

for i =1:(length(time)-1)
Ie = (100e-12)*sin(2*pi*f*(i+1)*dt);     % In cases injective current is time dependent, take Ie in n+1 timestep. 
v(1,i+1) = (v(1,i)+dt*(Em/Rm+Ie)/Cm)/(1+dt/(Rm*Cm));
end

%% Plot the V-curve
figure(1)
%clf;
plot(time,v);
legend('voltage');
title('intracellular potential of theda over t(s) with sinisoidal injective current.  ');
xlabel('time (s)');
ylabel('potential (V)');











