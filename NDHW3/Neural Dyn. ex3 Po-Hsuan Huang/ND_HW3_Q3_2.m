% Neural dynamics Homework3 
% Question 3-2
% Po-Hsuan Huang 2014,11,17
% The program use campartment model solving a patch of neuron with two
% compartment. One can designate different ingective current Ie between three
% time interval 0:Te Te:Ts Ts:T 
% The program print out the plot of the numerical result of the simulation.

close all;  % close all figures
%clc;       % clear command area
clear;     % clear variables
% a good hobby I learned today.
% control +C kills the process

%% Try different coefficient Ra1 Ra2 Ra3 = [7E6 265E6 30E9]

% Initail condition given by the qeustion. All units in MKS.
Em= 0;
Rm= 265E6;
Ra=30E9;
Cm= 75E-12;

% Euler's Method
dt= 0.1E-3;

% designating time intervals
Te = 0.4;
Ts = 0.44;
T = 1;
time= 0:dt:T;   

V1 = zeros(1,length(time));
V2 = zeros(1,length(time));

% initial conditions of the dynamic system
V1(1)= Em;
V2(1)= Em;

%% for time from 0 sec to Te sec


Ie = 0;

%Euler's method
for i= 1:Te/dt-1

V1(i+1) = V1(i)+ dt*((Em-V1(i))/Rm + (V2(i) -V1(i))/Ra  - Ie)/Cm;
V2(i+1) = V2(i)+ dt*((Em-V2(i))/Rm + (V1(i) -V2(i))/Ra  )/Cm;
end

%% for time from Te s to Ts s
Ie = -100E-12; % injective current = -100 pA

for i= Te/dt:Ts/dt

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

figure()
plot(time,V1,time,V2);
legend('V1','V2');
title('membrance potential of two compartment model, with injecting currnet on V1  ');
xlabel('time (s)');
ylabel('potential (V)');





