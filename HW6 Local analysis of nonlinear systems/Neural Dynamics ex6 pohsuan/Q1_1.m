% Neural Dynamics Exercise 6 
% Po-Hsuan Huang
% 2014/12/6
% Q1_1 
% modelling homogenious dynamic systems


%%

%% 
clear
close all

%% constants
A = [-0.5 -0.5 0; -0.5 -0.5 0; 0 0 2];  %  dynamics metrix for dx/dt = Ax


Initial = [1 1 0 ; 1 0 0 ; 0 1 0; 0 0 1e-6]; % four rows of initial conditions



T_t = 10 ; %   total modeling time.
dt = 1e-1 ; % time stepsize

figure(1)

for  trial = 1:4

x(:,1)= Initial(trial,:)';

[V D] = eig(A); % calculate the right eigenvectors and Diagonalized matrix.

%% Dynamical system solution

time = 1:dt:T_t;


for iter =1 : length(time)-1

    x(:,iter+1) = x(:,iter)+A*x(:,iter)*dt;
    
    
end

plot3(x(1,:),x(2,:),x(3,:),'o');
hold on

end

%% ploting solution

grid on
xlabel('x1');
ylabel('x2');
zlabel('x3');
legend('inital1','initial2','initial3','initial4')
title({'homogenious 3-D dynamic systems'})
   
hold off


