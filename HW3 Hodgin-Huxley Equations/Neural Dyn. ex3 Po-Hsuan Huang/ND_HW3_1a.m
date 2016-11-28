% Neural dynamics Homework3 
% Question 1
% Po-Hsuan Huang 2014,11,17



close all;  % close all figures
%clc;       % clear command area
clear;     % clear variables
% a good hobby I learned today.
% control +C kills the process


% initial conditions
sRm=1;              
sRa = 1;
sCm= 0.01;       

d =   [12E-6, 5E-6,4E-6];      % diameters of three neuron compartments
l = [180E-6,100E-6,100E-6];   % lengths of three neuron compartemnts.


Rm= sRm./(pi*d);
Ra= 4*sRa./(pi*d.^2);

lamda= sqrt(Rm./Ra);     % length constant
Rinf = sqrt(Rm.*Ra);     % semi-infinite resistance.

% Input resistance
 Rin = Rinf.*coth(l./lamda);


