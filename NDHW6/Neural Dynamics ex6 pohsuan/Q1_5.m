% Neural Dynamics Exercise 6 
% Po-Hsuan Huang
% 2014/12/17
% Q1_5
% plot vector field of the matrix

%% constants
close
clear

A = [-0.5 -0.5 0; -0.5 -0.5 0; 0 0 2];  %  dynamics metrix for dx/dt = Ax

x0 = [0; 0; 0 ]; % four rows of initial conditions
s0 = [1; 2 ;0];

T_t = 10 ; %   total modeling time.
dt = 1e-2 ; % time stepsize
%%  Trasform everything into Eigenspace

[V D] = eig(A); % calculate the right eigenvectors and Diagonalized matrix.

y0 = V'*x0;   %transformation of x
r0 = V'*s0;   %transfomration of s
y(:,1) = y0;
w(:,1)= x0;
display(y0,'y0');
display(s0,'s0');
%% Dynamical system solution

time = 1:dt:T_t;


for iter =1 : length(time)-1
    w(:,iter+1) = w(:,iter)+(A*w(:,iter)+s0)*dt;  
    % check if x= inv(V)*y transform can produce correct trajactory in realspace 
    y(:,iter+1) = y(:,iter)+(D*y(:,iter)+r0)*dt;
end
i = 1:1:length(time);
  x(:,i)= V*y(:,i);
% transform y back to real space.

figure(1)
scatter3(y(1,:),y(2,:),y(3,:),'filled'); hold on
%scatter3(x(1,:),x(2,:),x(3,:),'filled');  % trajactory in real space
scatter3(w(1,:),w(2,:),w(3,:),'filled');  % trajactory in real space

title('tragectory in eigen space with x = [0 0 0] ')
legend('eigenspace', 'real sapce')
xlabel('v1 [0.7 0.7 0]');
ylabel('v2 [0.7 -0.7 0]');
zlabel('v3 [0 0 1]');

hold off

figure(2)
scatter(y(2, :,:),y(3,:,:), 'filled');
title('tragectory in eigen space with x = [0 0 0] ')
xlabel('v2 [0.7 -0.7 0]');
ylabel('v3 [0 0 1]');

