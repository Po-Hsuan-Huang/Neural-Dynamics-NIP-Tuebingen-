% Neural dynamics Homework5
% Question 1-2
% Po-Hsuan Huang 2014,11,29
% The program use campartment model solving a neuron with 100 
% compartment.
% the Ingecting current is a step function of time and a delta function of
% space.
% The program print out the plot of the numerical result of the simulation.
% all units are in MKS system.
function NDHW5_Q1_2
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
gax = 0.5;   %e-3 
Em =0;       %e-3 initial voltage
dt = 0.025 ;   %e-3 time step size
Num  = 100;    %    number of compartment
T_t = 300;   %e-3 total modeling time 
V= Em*ones(1,Num);   % initial voltage 
    
%% specs of the injecting current

I0 = [6,8,15,20];           % amplitude of the injecting current in e-6 A;
je = 14;                              % the index of compartment being stimulated
te = 60;                           % the time when the stimulation starts
ts = 260;                         % the time when theh stimulation ends

Ie= @(t,j) I0(4)*Delta(j-je)*heaviside(t-te)*heaviside(ts-t);  % stimulating current


%% specs of channel dynamics
% create alpha for Na channal subunit m,h , K channel n.  
% alpha= {alpha_m,alpha_h, alpha_n}

% because you shifted the resting potential by +65mV, the initial condition
% for m,h,n are no longer zero.
Num_t= 0: dt :T_t;   % number of time steps 
m=zeros(length(Num_t),Num);
h=zeros(length(Num_t),Num);
n=zeros(length(Num_t),Num);
for j=1:Num
 % initial voltage of each compartment eaquals to Em
 v = V(1,j);
m(1,:)=  Alpha(1,v)/(Alpha(1,v)+Beta(1,v));
h(1,:)=  Alpha(2,v)/(Alpha(2,v)+Beta(2,v));
n(1,:)=  Alpha(3,v)/(Alpha(3,v)+Beta(3,v));
end
% alpha= { @(v) aDelta(v-25)*0.1*(v-25)/(1-exp(-((v-25)/10)))+Delta(v-25)*1 ...
%           ,@(v) 0.07*exp(-v/20)...
%           ,@(v) aDelta(v-10)*0.01*(v-10)/(1-exp(-((v-10)/10)))+Delta(v-10) };
% beta ={@(v) 4*exp(-v/18),@(v) 1/(1-exp(-(v-10)/10)),@(v) 0.125*exp(-v/80)};



%%Use farward Euler method to simulate each compartent.
 Num_t= 0: dt :T_t;   % number of time steps 
for step = 1: length(Num_t)-1   % time step iteration
    for j=1:Num                    % compartment iteration
        
    t= dt*step;
    v=V(step,j);
    m(step+1,j)= m(step,j)+ dt*(Alpha(1,v)*(1-m(step,j))-Beta(1,v)*m(step,j));
    h(step+1,j)= h(step,j)+ dt*(Alpha(2,v)*(1-h(step,j))-Beta(2,v)*h(step,j));
    n(step+1,j)= n(step,j)+ dt*(Alpha(3,v)*(1-n(step,j))-Beta(3,v)*n(step,j));
    
   if j ==1 
    %switch j
        %case 1
         V(step+1,j)= v-(dt/Cm)*(gL*(v-El)+gNa*m(step,j)^3*h(step,j)...
         *(v-Ena)+gK*n(step,j)^4*(v-Ek)+gax*(2*v-V(step,j+1)-V(step,j+1))-Ie(t,j));
    

%          V(step+1,j)=v+(dt/Cm)*(Ie(t,j)-gL*(v-El)-gNa*m(step)^3*h(step)...
%          *(v-Ena)-gK*n(step)^4*(v-Ek));
    
   elseif j==Num     
    %case Num
           v=0;
            V(step+1,j)=v;
       else     
        %otherwise
            V(step+1,j)= v-(dt/Cm)*(-Ie(t,j)+gL*(v-El)+gNa*m(step,j)^3*h(step,j)...
                *(v-Ena)+gK*n(step,j)^4*(v-Ek)+gax*(2*v-V(step,j+1)-V(step,j-1)));
    end    
    end
end

%% plot
figure(1)
Num_x = 1:Num;
 surf(Num_x,Num_t,V);   %(X(j), Y(i), Z(i,j))
%plot(Num_t, V)
colorbar
axis tight;
%zlim([-20 120]);
title({'Multicompartment model' ,'when injecting step current on j=14 '});
xlabel('compartment');
ylabel('time (ms)');
zlabel('Voltage(mV)');
shading interp;

% ylim([-20 120]);

figure(2)


subplot(1,3,1)
 surf(Num_x,Num_t,m);
 title({'channel subunit vs time'});
xlabel('compartment');
ylabel('time (ms)');
legend('m');
shading interp;
subplot(1,3,2) 
 surf(Num_x,Num_t,h);
 title({'channel subunit vs time'});
xlabel('compartment');
ylabel('time (ms)');
legend('h');
shading interp;
subplot(1,3,3)
 surf(Num_x,Num_t,n);
title({'channel subunit vs time'});
xlabel('compartment');
ylabel('time (ms)');
legend('n');
shading interp;

end



%% funtions

% % setGlobalv
% function setGlobalv(value)
%   global v;
%   v = value;
% end  
% 
% function value = getGlobalv
%   global v;
% value = v;
% end


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