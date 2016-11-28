% Neural dynamics Homework3 
% Question 3-3
% Po-Hsuan Huang 2014,11,17
% The program use campartment model solving a patch of neuron with two
% compartment. One can designate different ingective current Ie between three
% time interval 0:Te Te:Ts Ts:T 
% The program print out the plot of the numerical result of the simulation.
% In this quietion, you just play with the coefficient of the frequency of
% the injective sinusoidal current.


%% clear mess 
close all;  % close all figures
%clc;       % clear command area
clear;     % clear variables
% a good hobby I learned today.
% control +C kills the process

%% 

% Initail condition given by the qeustion. All units in MKS.
Em= 0;
Rm= 265E6;
Ra=300E6;
Cm= 75E-12;

% Euler's Method
dt= 0.1E-3;

% designating time intervals
Te = 0.4;
Ts = 10.4;
T = 10.4;
time= 0:dt:T;   

V1 = zeros(1,length(time));
V2 = zeros(1,length(time));

% initial conditions of the dynamic system
V1(1)= Em;
V2(1)= Em;
f = [1 2 5 10 20 50 100 200 500 1e3 2e3 5e3];  % frequency you can play with ^^
fq = f(1);

for j = 1:length(f)
%% for time from 0 sec to Te sec

fq = f(j);

Ie = 0;

%Euler's method
for i= 1:Te/dt-1

V1(i+1) = V1(i)+ dt*((Em-V1(i))/Rm + (V2(i) -V1(i))/Ra  - Ie)/Cm;
V2(i+1) = V2(i)+ dt*((Em-V2(i))/Rm + (V1(i) -V2(i))/Ra  )/Cm;
end

%% for time from Te s to Ts s

for i= Te/dt:Ts/dt
Ie = -(100E-12)*sin(2*pi*fq*i*dt); % injective current = -100 pA* sin(2pi*f*t)

V1(i+1) = V1(i)+ dt*((Em-V1(i))/Rm + (V2(i) -V1(i))/Ra  - Ie)/Cm;
V2(i+1) = V2(i)+ dt*((Em-V2(i))/Rm + (V1(i) -V2(i))/Ra  )/Cm;
end

%% for ttme form Ts to Ts

Ie = 0; % injective current = -100 pA

for i= Ts/dt:T/dt-1

V1(i+1) = V1(i)+ dt*((Em-V1(i))/Rm + (V2(i) -V1(i))/Ra  - Ie)/Cm;
V2(i+1) = V2(i)+ dt*((Em-V2(i))/Rm + (V1(i) -V2(i))/Ra  )/Cm;
end


%% plotting the result
V1max = max(V1(0.8*Ts/dt:Ts/dt));
V1min = min(V1(0.8*Ts/dt:Ts/dt));
Amp(1,j) = (V1max-V1min)/2;
V2max = max(V2(0.8*Ts/dt:Ts/dt));
V2min = min(V2(0.8*Ts/dt:Ts/dt));
Amp(2,j) = (V2max-V2min)/2;

end
figure(1)
FqAmp(1,:)= f;
FqAmp(2,:)= Amp(1,:);
FqAmp(3,:)= Amp(2,:);

loglog(FqAmp(1,:),FqAmp(2,:),'*-');
hold on
loglog(FqAmp(1,:),FqAmp(3,:),'o-');
legend('AmpV1', 'AmpV2');
title({'Amplitudes change over freqency of Vm of two compartment model' ,'when injecting sinusoidal current on V1 '});
xlabel('frequency (Hz)');
ylabel('Amplitude of Vm (V)');
hold off





