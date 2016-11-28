% Neural dynamics Homework4
% Question 2-2
% Po-Hsuan Huang 2014,11,23
% The program model single compartment Hudgin-Huxley model, 
% Plot the relation of gating particle opening rate and the clamping
% voltage.
% all units are in MKS system.
% clamping voltage increase incrementally. 
function NDHW4_Q2_3
%% clear mess 
close all;  % close all figures
%clc;       % clear command area
clear;     % clear variables
%% specs of the neuron

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
setGlobalv(V)      % initial voltage of each compartment eaquals to Em

%% specs of the injecting current
% oberving the relationship between firing rate and injecting current.
% Inject current from 0~40mA with steps of 4 mA.
% we spend 20 ms for each mesurement, and only use the post 40ms for 
% firing rate calculation.

I_start = 0; % in mA
I_end = 20;%in mA
I_step = 0.5; % mA
I_int= (I_end-I_start)/I_step; % how many time intervals we need. 
T_t = 40*I_int;   %e-3 total modeling time 

je = 1;                              % the index of compartment being stimulated
I_record(1,1)=0;



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
   
      % take the integer part of t/20 to know which of the 40 interval is
      % it, and 4*fix(t/20) is the injecting current amp in that interval.
     Ie= @(t,j)0.5*(fix(t/I_int))*Delta(j-je)*heaviside(t-fix(t/I_int)*I_int);  % stimulating current
    I_record(1,step+1)=Ie(t,j);
    m(step+1)= m(step)+ dt*(Alpha(1)*(1-m(step))-Beta(1)*m(step));
    h(step+1)= h(step)+ dt*(Alpha(2)*(1-h(step))-Beta(2)*h(step));
    n(step+1)= n(step)+ dt*(Alpha(3)*(1-n(step))-Beta(3)*n(step));
    V(step+1,j)= V(step,j)+(dt/Cm)*(Ie(t,j)-gL*(V(step,j)-El)-gNa*m(step)^3*h(step)*(V(step,j)-Ena)-gK*n(step)^4*(V(step,j)-Ek));
    setGlobalv(V(step+1,j));
    end
end


%% Calculate firing rate in each time interval.
thres = 80; %  threshold for detecting fire.
FR(1)=0;    % firing rate for no injecting current is zero.
for p = 1: I_int
T_int= T_t/I_int;  %length of each interval ms   

% bad algorithm here------
% %firing rate  (0.5* points between [threshold threshold+0.01])
% % for element value >= thres, return 1, otherwise 0
% Upper= V((1+(2*p-1)*  0.5*T_int/dt  ):1+(2*p*  0.5*T_int )/dt,j) >= thres;
% % for element value<= thres+0.01, return 1, otherwise 0
% Lower= V((1+(2*p-1)*  0.5*T_int/dt  ):1+(2*p*  0.5*T_int )/dt,j) <= thres+2; 
% 
% % if element between Upper and Lower, return 1, otherwise 0.
% % divide by 2 since two points consist a spike.
% FR(p+1) = 0.5*sum(Upper+Lower==2)/(T_int);    



% Boarder detecting techniques
Window1 = V((2*p-1)*  0.5*T_int/dt  :(2*p*  0.5*T_int )/dt,j) >= thres;
Window2 = V(1+(2*p-1)* 0.5*T_int/dt  :1+(2*p* 0.5*T_int )/dt,j)>= thres;
FR(p+1)=  0.5*sum(abs(Window1-Window2))/(T_int);

end
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
I_space= I_start:I_step:I_end;
plot(I_space,FR,'-o');
title('firing_rate')
xlabel('injecting current intensity (mA)')
ylabel('firing rate  (1/ms)')
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