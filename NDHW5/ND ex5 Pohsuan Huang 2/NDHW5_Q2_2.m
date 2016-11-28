% Neural dynamics Homework4
% Question 2-2
% Po-Hsuan Huang 2014,11,23
% The program model A type current, single compartment Hudgin-Huxley model , 
% A-type current incorporate potassium channel. 
% all units are in MKS system.
function NDHW5_Q2_2
%% clear mess 
close all;  % close all figures
%clc;       % clear command area
clear;     % clear variables
%% specs of the neuron

Cm  = 1;   % e-6
Ena= 50;  % e-3
Ek = -77;   % e-3
Ea = -80;
El = -22;  %e-3
gNa = 120;  %e-3
gK = 20;    %e-3
gL = 0.3;   %e-3
gA = 47.7;
Em =-73;       %e-3 initial voltage
dt = 0.025 ;   %e-3 time step size
Num  = 1;    %    number of compartment
T_t = 460;   %e-3 total modeling time 
V= Em*ones(1,1);   % initial voltage 

%% specs of the injecting current

I0 = [1,20];           % amplitude of the injecting current in e-6 A;
je = 1;                              % the index of compartment being stimulated
te = 60;                           % the time when the stimulation starts
ts = 460;                         % the time when theh stimulation ends


Inject_x = I0(1):0.5:I0(2);

for trail = 1: length(Inject_x)
Ie= @(t,j) (I0(1)*trail*0.5)*Delta(j-je)*heaviside(t-te)*heaviside(ts-t);  % stimulating current

v = V(1);
%% specs of channel dynamics
% create alpha for Na channal subunit m,h , K channel n.  
% alpha= {alpha_m,alpha_h, alpha_n}

% because you shifted the resting potential by +65mV, the initial condition
% for m,h,n are no longer zero.
m(1)=  Alpha(1,v)/(Alpha(1,v)+Beta(1,v));
h(1)=  Alpha(2,v)/(Alpha(2,v)+Beta(2,v));
n(1)=  Alpha(3,v)/(Alpha(3,v)+Beta(3,v));
a(1) = Alpha(4,v)/(Alpha(4,v)+Beta(4,v));
b(1) = Alpha(5,v)/(Alpha(5,v)+Beta(5,v));
% alpha= { @(v) aDelta(v-25)*0.1*(v-25)/(1-exp(-((v-25)/10)))+Delta(v-25)*1 ...
%           ,@(v) 0.07*exp(-v/20)...
%           ,@(v) aDelta(v-10)*0.01*(v-10)/(1-exp(-((v-10)/10)))+Delta(v-10) };
% beta ={@(v) 4*exp(-v/18),@(v) 1/(1-exp(-(v-10)/10)),@(v) 0.125*exp(-v/80)};



%%Use farward Euler method to simulate each compartent.
 Num_t= 0: dt :T_t;   % number of time steps 
for step = 1: length(Num_t)-1   % time step iteration
    for j=1:Num                    % compartment iteration
    t= dt*step;
    v= V(step,j);
    
    m(step+1)= m(step)+ dt*(Alpha(1,v)*(1-m(step))-Beta(1,v)*m(step));
    h(step+1)= h(step)+ dt*(Alpha(2,v)*(1-h(step))-Beta(2,v)*h(step));
    n(step+1)= n(step)+ dt*(Alpha(3,v)*(1-n(step))-Beta(3,v)*n(step));
    a(step+1)= a(step)+ dt*(Alpha(3,v)*(1-a(step))-Beta(3,v)*a(step));
    b(step+1)= b(step)+ dt*(Alpha(3,v)*(1-b(step))-Beta(3,v)*b(step));
    
    V(step+1,j)= V(step,j)+(dt/Cm)*(Ie(t,j)-gL*(V(step,j)-El)-gNa*...
        m(step)^3*h(step)*(V(step,j)-Ena)-gK*n(step)^4*(V(step,j)-Ek)-gA*...
        a(step)^3*b(step)*(V(step,j)-Ea));
    
    %V(step+1,j)= V(step,j)+(dt/Cm)*(Ie(t,j)-gL*(V(step,j)-El)-gNa*...
      %m(step)^3*h(step)*(V(step,j)-Ena)-gK*n(step)^4*(V(step,j)-Ek));%-gA*a(step)^3*b(step)*(V(step,j)-Ea));

    end
end
%% find peaks

 [pks,locs]=findpeaks(V(60/dt:460/dt,1));
 if length(locs)~=0;
 freq(trail) = 1000*length(locs)/((locs(end)-locs(1))*dt);  % firing rate in Hz

 elseif length(locs)==0
 
 freq(trail)=0;
 end
end
%% plot
figure(1)
plot(Inject_x,freq);
title({'single compartment Hodgin-Huxley model with A type channel' ,'firing rate- injecting current '});
xlabel('Injecting current (microA)');
ylabel('Firing rate (Hz)');


figure(2)
subplot(2,1,1)
plot(Num_t,V);   

title({'single compartment Hodgin-Huxley model with A type channel' ,'with step injection current '});
xlabel('time (ms)');
ylabel('Voltage (mV)');
% ylim([-20 120]);

subplot(2,1,2)
plot(Num_t,m,Num_t,h,Num_t,n,Num_t,a,Num_t,b);
title({'channel subunit vs time'});
xlabel('time (s)');
ylabel('open posibility');
legend('m','h','n','a','b');


end



%% funtions



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
         value= aDelta(v+34.7)*3.8*0.1*(v+34.7)/(1-exp(-(v+34.7)/10))+Delta(v+34.7)*3.8;
      case 2   % channel h
          value= 3.8*0.07*exp(-(v+53)/20);
      case 3   % channel n
          value= (3.8/2)*aDelta(v+50.7)*0.01*(v+50.7)/(1-exp(-(v+50.7)/10))+Delta(v+50.7)*0.1*(3.8/2) ;
      case 4   % channel a    
           inf4=( 0.0761*exp((v+99.22)/31.84)/(1+exp((v+6.17)/28.93)))^(1/3);
           tau4= 0.3632+1.158/(1+exp((v+60.96)/20.12));
           value = inf4/tau4;
      case 5    % channel b
          inf5= (1+exp((v+58.3)/14.54))^-4;
          tau5 = 1.24+ 2.678/(1+exp((v-55)/16.027));
          value = inf5/tau5;
  end        
end
% beta factor for channel dynamics
function value=Beta(channeltype,v)


  switch channeltype
      case 1 % channel m
           value=3.8*4*exp(-(v+59.7)/18);
      case 2 % channel h
          value= 3.8/(1+exp(-(v+23)/10));
      case 3  % channel n
           value= (3.8/2)*0.125*exp(-(v+60.7)/80);
      case 4   % channel a   
            inf4=( 0.0761*exp((v+99.22)/31.84)/(1+exp((v+6.17)/28.93)))^(1/3);
            tau4= 0.3632+1.158/(1+exp((v+60.96)/20.12));
            value = (1-inf4)/tau4;
      case 5    % channel b
            inf5= (1+exp((v+58.3)/14.54))^-4;
            tau5 = 1.24+ 2.678/(1+exp((v-55)/16.027));
            value = (1-inf5)/tau5;
  end
end