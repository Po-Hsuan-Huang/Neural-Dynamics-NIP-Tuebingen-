% Neural dynamics Homework4
% Question 2-2
% Po-Hsuan Huang 2014,11,23
% The program model single compartment Hudgin-Huxley model, 
% Plot the relation of gating particle opening rate and the clamping
% voltage.
% all units are in MKS system.
% clamping voltage increase incrementally.

function NDHW4_Q2_3_2
%% clear mess 
close all;  % close all figures
%clc;       % clear command area
clear;     % clear variables
%% specs of the neuron
tic;

% check the dynamic differential eqaution of HH model, you will find out 
% all units cencel out if we use ms,mV,microF, micro amp, and micro
% conductance. Same for the gating particle and its rate factors alpha
% beta.

Cm  = 1;   % e-6
Ena= 115;  % e-3
Ek = -12;   % e-3
El = 10.6;  %e-3
gNa = 120;  %e-3
gK = 36;    %e-3
gL = 0.3;   %e-3

Em =0;       %e-3 initial voltage
dt = 0.025 ;   %e-3 time step size
Num  = 1;    %    number of compartment
V= Em*ones(1,1);   % initial voltage 
% setGlobalv(V)      % initial voltage of each compartment eaquals to Em
v=V(1,1);
%% specs of the injecting current
% oberving the relationship between firing rate and injecting current.
% Inject current from 0~40mA with steps of 4 mA.
% we spend 20 ms for each mesurement, and only use the post 40ms for 
% firing rate calculation.

I_start = 0; % in mA
I_end = 20;%in mA
I_step = 0.5; % mA
I_int= (I_end-I_start)/I_step; % how many time intervals we need. 
T_t = 300;   %e-3 total modeling time 
te= 50;
ts =300; 
je = 1;                              % the index of compartment being stimulated
I_record(1,1)=0;

thres = 80; %  threshold for detecting fire.
FR(1)=0;    % firing rate for no injecting current is zero.
T_int= ts-te;  %length of each interval ms   
I_space=I_start:I_step:I_end;
toc
for p =1:length(I_space)
%% specs of channel dynamics
% create alpha for Na channal subunit m,h , K channel n.  
% alpha= {alpha_m,alpha_h, alpha_n}

% because you shifted the resting potential by +65mV, the initial condition
% for m,h,n are no longer zero.
m(1)=  Alpha(1,v)/(Alpha(1,v)+Beta(1,v));
h(1)=  Alpha(2,v)/(Alpha(2,v)+Beta(2,v));
n(1)=  Alpha(3,v)/(Alpha(3,v)+Beta(3,v));
% alpha= { @(v) aDelta(v-25)*0.1*(v-25)/(1-exp(-((v-25)/10)))+Delta(v-25)*1 ...
%           ,@(v) 0.07*exp(-v/20)...
%           ,@(v) aDelta(v-10)*0.01*(v-10)/(1-exp(-((v-10)/10)))+Delta(v-10) };
% beta ={@(v) 4*exp(-v/18),@(v) 1/(1-exp(-(v-10)/10)),@(v) 0.125*exp(-v/80)};


Ie= @(t,j) p*0.5*Delta(j-je)*heaviside(t-te)*heaviside(ts-t);  % stimulating current

%%Use farward Euler method to simulate each compartent.
 Num_t= 0: dt :T_t;   % number of time steps 
for step = 1: length(Num_t)-1   % time step iteration
    for j=1:Num                    % compartment iteration
    t= dt*step;
       v=V(step,j);
      % take the integer part of t/20 to know which of the 40 interval is
      % it, and 4*fix(t/20) is the injecting current amp in that interval.
    m(step+1)= m(step)+ dt*(Alpha(1,v)*(1-m(step))-Beta(1,v)*m(step));
    h(step+1)= h(step)+ dt*(Alpha(2,v)*(1-h(step))-Beta(2,v)*h(step));
    n(step+1)= n(step)+ dt*(Alpha(3,v)*(1-n(step))-Beta(3,v)*n(step));
    V(step+1,j)= V(step,j)+(dt/Cm)*(Ie(t,j)-gL*(V(step,j)-El)-gNa*m(step)^3*h(step)*(V(step,j)-Ena)-gK*n(step)^4*(V(step,j)-Ek));
    end
end  
    
    %% Calculate firing rate in each time interval.



% Boarder detecting techniques
Window1 = V( (10+te)/dt  : ts/dt,j) >= thres;
Window2 = V(1+(10+te)/dt  : 1+ts/dt,j)>= thres;
FR(p)=  0.5*sum(abs(Window1-Window2))/(T_int);

% or you can use funciton findpeak, and find the average period, more
% accurate.

toc
end
toc
%% plot
figure
subplot(3,1,1)
plot(Num_t,V);   

title({'single compartment Hodgin-Huxley model' ,'with step injection current '});
xlabel('time (ms)');
ylabel('Voltage (mV)');
ylim([-20 120]);

subplot(3,1,2)
plot(Num_t,m,Num_t,h,Num_t,n);
title({'channel subunit vs time'});
xlabel('time (ms)');
ylabel('open posibility');
legend('m','h','n');

subplot(3,1,3)
plot(Num_t,I_record)
xlabel('time (ms)');
ylabel('Injecting current(mA)');


% plot the firing rate, current relation.
figure(2)
plot(I_space,FR,'-o');
title('firing_rate')
xlabel('injecting current intensity (uA)')
ylabel('firing rate  (1/ms)')

toc
end



%% funtions

% setGlobalv
% function setGlobalv(value)
%   global v;
%   v = value;
% end  
% 
% function value = getGlobalv
%   global v;
% value = v;
% end
% 

% delta
function  value=Delta(j)
    value =logical(j==0);
end
% adelta
function value= aDelta(j)
    value= logical(Delta(j)==0);
end
% alpha factor for channel dynamics
function value=Alpha(channeltype,v)
  
   
   
  switch channeltype
      case 1   % channel m
         value= aDelta(v-25)*0.1*(v-25)/(1-exp(-(v-25)/10))+Delta(v-25)*1;
      case 2   % channel h
          value= 0.07*exp(-v/20);
      case 3   % channel n
          value= aDelta(v-10)*0.01*(v-10)/(1-exp(-(v-10)/10))+Delta(v-10)*0.1 ;
  end        
end
% beta factor for channel dynamics
function value=Beta(channeltype,v)


  switch channeltype
      case 1 % channel m
           value=4*exp(-v/18);
      case 2 % channel h
          value= 1/(1+exp(-(v-30)/10));
      case 3  % channel n
           value= 0.125*exp(-v/80);
  end
end