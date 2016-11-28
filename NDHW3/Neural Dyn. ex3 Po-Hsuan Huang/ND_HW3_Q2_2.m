% Neural dynamics Homework3 
% Question 2-3
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

d =   [16E-6, 5.37E-6,2.92E-6,3.84E-6];      % diameters of three neuron compartments
l = [210E-6,150E-6,100E-6,120E-6];   % lengths of three neuron compartemnts.


Rm= sRm./(pi*d);
Ra= 4*sRa./(pi*d.^2);

lamda= sqrt(Rm./Ra);     % length constant
Rinf = sqrt(Rm.*Ra);     % semi-infinite resistance.

% Input resistance
 Rin = Rinf.*coth(1./lamda);
 
d2= ( 16^(3/2)/((1+l(3)/l(2))^3+(l(4)/l(2))^3))^(2/3);
d3= d2*(l(3)/l(2))^(3/2)
d4= d2*(l(4)/l(2))^(3/2)

equil2= l(2)*lamda(1)/lamda(2)   % equivalent length of compartment2

