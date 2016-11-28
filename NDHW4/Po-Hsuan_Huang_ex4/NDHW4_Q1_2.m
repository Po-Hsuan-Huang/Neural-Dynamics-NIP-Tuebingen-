% Neural dynamics Homework4
% Question 1-2
% Po-Hsuan Huang 2014,11,23
% The program use campartment model solving a neuron with 50 
% compartment.
% the Ingecting current is a step function of time and a delta function of
% space.
% The program print out the plot of the numerical result of the simulation.
% all units are in MKS system.

%% clear mess 
close all;  % close all figures
%clc;       % clear command area
clear;     % clear variables
%% specs of the neuron

Cm = 62.8e-12;  % capacity of each neuron segment
Rm = 1.59e9  ;  % membrane Resistance of each neuron segment
Ra =  0.0318e9; % axial Resistance of each neuron segment
Em = 0         ;% rest potential/ reversal potential of the membrane
dt = 0.02e-3 ;% time step must be smaller than 0.02e-3
Num  = 50;        % number of compartment
T_t = 1000e-3;     % total modeling time 
V=  Em*ones(1,50); % initial voltage of each compartment eaquals to Em 
%% specs of the injecting current

I0 = 10e-12;     % amplitude of the injecting current
je = 20;         % the index of compartment being stimulated
te = 20e-3;      % the time when stimuation is applied

Delta= @(j) logical(j==je);
Ie= @(t,j) I0*Delta(j)*heaviside(t-te);


%% Use farward Euler method to simulate each compartent.
 Num_t= 0: dt :T_t;   % number of time steps 
for step = 1: length(Num_t)-1   % time step iteration
for j = 1:Num                 % neuron iteration from origin j=1 to terminal j=50
    
    t= dt*step;
    
  switch j
   %if j==1   
      case 1  % the first neuron is sealed end
        V(step+1,j)=  V(step,j)+ dt* (-V(step,j)/Rm + (2*V(step,j+1)-2*V(step,j))/Ra+ Ie(t,j))/Cm;
      
%       case 49   % encounter the killed end
%          V(step+1,j)=  V(step,j)+ dt* (-V(step,j)/Rm + (V(step,j+1)+V(step,j-1)-2*V(step,j))/Ra+ Ie(t,j))/Cm;
     
 %  elseif j==50
        case 50   % the last neuron is killed end, always 0 volt
          V(step,j) = 0;
          V(step+1,j)= V(step,j);
   %else  
        otherwise % compartment in between
        V(step+1,j)=  V(step,j)+ dt* (-V(step,j)/Rm + (V(step,j+1)+V(step,j-1)-2*V(step,j))/Ra+ Ie(t,j))/Cm; 
    
  end
end  
end

% By now we should get a  matrix V of  Step x Num

figure
Num_x = 1:Num;
surf(Num_x,Num_t,V);   %(X(j), Y(i), Z(i,j))
colormap summer;
colorbar
axis tight;
title({'Multicompartment model' ,'when injecting step current on j=20 '});
xlabel('compartment');
ylabel('time (s)');
shading interp


figure(2)
plot(V(step+1,:),'o')
title('steady state')
xlabel('compartment');
ylabel('voltage');



