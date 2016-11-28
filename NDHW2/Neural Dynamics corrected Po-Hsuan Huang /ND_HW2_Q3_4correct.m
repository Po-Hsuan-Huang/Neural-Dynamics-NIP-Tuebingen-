% Neural Dynamics Exercise 2
% Problem 3(4)
% Po-Hsuan Huang 2014.11.11
% Simulation of single compartment model.

%% specs of modeled neuro compartment
% close all;  % close all figures
% clc;       % clear command area
% clear;     % clear variables
% a good hobby I learned today.
% control +C kills the process

Em = 0;            % in volt
length = 100E-6 ;     % in meter
diameter = 2E-6  ;    % in meter
sRm = 1  ;         % specific memebrane resistance in ohm*sqrmeter
sRa = 1  ;        % specific axial resistivity in ohm*mete
sCm = 0.01 ;       % specific memberane capacitance in Farad/sqrmeter
f =1000;              % ingective current frequency
%% calculating the parameter in analogic circuit.
Rm =  sRm/(pi*diameter*length);     % membrane resistance in ohm.
Cm =  sCm*pi*diameter*length;                       % capacitance in Farad.
Ie = -50e-12;                                % injective constant current in A.


%% backward Euler Method
n = 100000;   % total timesteps.
dt= 0.0001;   % timestep in second. must small enough to plot correctly.
Vamp= zeros(1,1000);
v= zeros(1,n);

for f = 1:1000-1
    v(1,1) = Em; % construct a matrix to store datapoints.
for i =1:n
Ie = (100e-12)*sin(2*pi*f*(i+1)*dt);     % In cases injective current is time dependent, take Ie in n+1 timestep. 
v(1,i+1) = (v(1,i)+dt*(Em/Rm+Ie)/Cm)/(1+dt/(Rm*Cm));
end

%% Plot the Bode diagram for Q3.4.

%retract the sub matrix of converged voltage values.
conV = v(1, 80000:100000);
Vmax = (max(conV));
Vmin = (min(conV));
Vamp(f)= (Vmax-Vmin)/2;

end
figure(1)
fig = loglog(Vamp);
title(  'loglog voltage-frequency');
xlabel('frequency');
ylabel('amplitude');