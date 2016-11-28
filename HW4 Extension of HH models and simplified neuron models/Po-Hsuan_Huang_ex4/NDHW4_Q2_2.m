% Neural dynamics Homework4
% Question 2-2
% Po-Hsuan Huang 2014,11,23
% The program model single compartment Hudgin-Huxley model, 
% with square injecting current
% all units are in MKS system.
function NDHW4_Q2_2
%% clear mess 
close all;  % close all figures
%clc;       % clear command area
clear;     % clear variables
%% specs of the neuron

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
T_t = 200;   %e-3 total modeling time 
V= Em*ones(1,1);   % initial voltage 
setGlobalv(V)      % initial voltage of each compartment eaquals to Em

%% specs of the injecting current

I0 = [0,3,6,8];           % amplitude of the injecting current in e-6 A;
je = 1;                              % the index of compartment being stimulated
te = 50;                           % the time when the stimulation starts
ts = 300;                         % the time when theh stimulation ends

Ie= @(t,j) I0(3)*Delta(j-je)*heaviside(t-te)*heaviside(ts-t);  % stimulating current


%% specs of channel dynamics
% create alpha for Na channal subunit m,h , K channel n.  
% alpha= {alpha_m,alpha_h, alpha_n}

% because you shifted the resting potential by +65mV, the initial condition
% for m,h,n are no longer zero.
m(1)=  Alpha(1)/(Alpha(1)+Beta(1));
h(1)=  Alpha(2)/(Alpha(2)+Beta(2));
n(1)=  Alpha(3)/(Alpha(3)+Beta(3));
% alpha= { @(v) aDelta(v-25)*0.1*(v-25)/(1-exp(-((v-25)/10)))+Delta(v-25)*1 ...
%           ,@(v) 0.07*exp(-v/20)...
%           ,@(v) aDelta(v-10)*0.01*(v-10)/(1-exp(-((v-10)/10)))+Delta(v-10) };
% beta ={@(v) 4*exp(-v/18),@(v) 1/(1-exp(-(v-10)/10)),@(v) 0.125*exp(-v/80)};



%%Use farward Euler method to simulate each compartent.
 Num_t= 0: dt :T_t;   % number of time steps 
for step = 1: length(Num_t)-1   % time step iteration
    for j=1:Num                    % compartment iteration
    t= dt*step;
    m(step+1)= m(step)+ dt*(Alpha(1)*(1-m(step))-Beta(1)*m(step));
    h(step+1)= h(step)+ dt*(Alpha(2)*(1-h(step))-Beta(2)*h(step));
    n(step+1)= n(step)+ dt*(Alpha(3)*(1-n(step))-Beta(3)*n(step));
    V(step+1,j)= V(step,j)+(dt/Cm)*(Ie(t,j)-gL*(V(step,j)-El)-gNa*m(step)^3*h(step)*(V(step,j)-Ena)-gK*n(step)^4*(V(step,j)-Ek));
    setGlobalv(V(step+1,j));
    end
end

%% plot
figure
subplot(2,1,1)
plot(Num_t,V);   

title({'single compartment Hodgin-Huxley model' ,'with step injection current '});
xlabel('time (ms)');
ylabel('Voltage (mV)');
% ylim([-20 120]);

subplot(2,1,2)
plot(Num_t,m,Num_t,h,Num_t,n);
title({'channel subunit vs time'});
xlabel('time (s)');
ylabel('open posibility');
legend('m','h','n');


end



%% funtions

% setGlobalv
function setGlobalv(value)
  global v;
  v = value;
end  

function value = getGlobalv
  global v;
value = v;
end


% delta
function  value=Delta(j)
    value =logical(j==0);
end
% adelta
function value= aDelta(j)
    value= logical(Delta(j)==0);
end
% alpha factor for channel dynamics
function value=Alpha(channeltype)
  
   v = getGlobalv;
   
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
function value=Beta(channeltype)

   v = getGlobalv;

  switch channeltype
      case 1 % channel m
           value=4*exp(-v/18);
      case 2 % channel h
          value= 1/(1+exp(-(v-30)/10));
      case 3  % channel n
           value= 0.125*exp(-v/80);
  end
end