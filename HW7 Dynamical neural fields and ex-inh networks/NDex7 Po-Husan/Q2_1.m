% Neural Dynamics Exercise 7
% Po-Hsuan Huang
% 2015/1/13
% Q2_1
% Find the fixed points of the matrix in auto-assoicative memory Neuron Network

%%

%% 
clear
close all

%% constants

M = [ 1 -0.1 -0.1; -0.1 1 -0.1; -0.1 -0.1 -0.1];  %  dynamics metrix for du/dt = -u + [Mu]
K = [ 0 -0.1 -0.1; -0.1 0 -0.1; -0.1 -0.1 -1.1];  %  dynamics metrix for du/dt = Au

[V ,D, W] = eig(K); % calculate the right eigenvectors and Diagonalized matrix.

display (W,'right eigenvectors')

display (D,'diagonalized metrix')


tau = 1;



% Q2_1

% 
% % 3d quiver
% 
% 
% % faster way to plot, but didn't figure out how to do it correctly
% % fix on axis, and plot the rest two.
% point= 0:2.5:10;
% 
% clear z
% 
% for j = 1: length(point)
% for i = 1: length(point)
% for k = 1: length(point)
%     
%     
% thres =   M*[point(i);point(j);point(k)];
% 
% 
% zen(i,j,k,:) =(1/tau)*(-[point(i);point(j);point(k)]+ (thres > 0).*thres);
% 
% 
% end
% end
% end
% [x,y,w]=meshgrid(0:2.5:10);
% 
% figure(10)
% quiver3(x,y,w,zen(:,:,:,1),zen(:,:,:,2),zen(:,:,:,3)); % projection on x3 =0
% xlabel('x1');
% ylabel('x2');
% zlabel('x3');
% title({'velocity field'})
% axis tight
% 
% Plane x-y
[x,y]=meshgrid(-1:1:10);

point= -1:1:10;

for i = 1: length(point)
for j = 1: length(point)
    
thres =     M*[point(i);point(j);0];
z(j,i,:) =(1/tau)*(-[point(i);point(j);0]+ (thres > 0).*thres);

Energy (j,i) = sqrt(z(j,i,1)^2+ z(j,i,2)^2 + z(j,i,3)^2 );

end
end

figure(2)
quiver(x,y,z(:,:,1),z(:,:,2));% proection on x3 =0
hold on
surf(x,y,Energy);
hold off
xlabel('x1');
ylabel('x2');
title({'velocity field x3 = 0'})
axis tight

% clear z;
% 
% % Plane z-x
% 
% for i = 1: length(point)
% for j = 1: length(point)
%     
% thres =     M*[point(i); 0 ; point(j)];
% z(j,i,:) =(1/tau)*(-[point(i);0;point(j)]+ (thres > 0).*thres);
% end
% end
% 
% figure(3)
% quiver(x,y,z(:,:,3),z(:,:,1));% proection on x2 =0
% xlabel('x1');
% ylabel('x3');
% title({'velocity field x2 = 0'})
% axis tight
% 
% clear z
% 
% % Plane y-z
% 
% for i = 1: length(point)
% for j = 1: length(point)
%     
% thres =     M*[0;point(i);point(j)];
% z(j,i,:) =(1/tau)*(-[0;point(i);point(j)]+ (thres > 0).*thres);
% end
% end
% 
% figure(4)
% quiver(x,y,z(:,:,2),z(:,:,3));% proection on x1 =0
% xlabel('x2');
% ylabel('x3');
% title({'velocity field x1 = 0'})
% axis tight

%% Study of trajectory



T_t = 4 ; %   total modeling time.

dt = 1e-2 ; % time stepsize

tau = 1;



time = 0:dt:T_t;

Initial = [ 1 1 0 ; -1 -1 0 ; -1 1 0; 1 -1 0];


Initial(5:8,:)= [-2 -3 0 ; -3 -2 0 ; 2 3 0 ; 3 2 0];   

%Initial (1:4,:)= zeros(4,3);


M = [ 1 -0.1 -0.1;
    -0.1 1 -0.1; -0.1 -0.1 -0.1];  %  dynamics metrix for dx/dt = Ax



figure(1)

for trial =1 :8

u(:,1)= Initial(trial,:)';

for iter =1 : length(time)-1

    u(:,iter+1) = u(:,iter)+ dt*( -u(:,iter) + (M*u(:,iter)> 0)'*M*u(:,iter))/tau;
    
    
end

plot3(u(1,:),u(2,:),u(3,:),'o');


hold on

end


grid on
xlabel('u1');
ylabel('u2');
zlabel('u3');
legend('inital1','initial2','initial3','initial4','initial5','initial6','initial7','initial8' )
title({'homogenious 3-D dynamic systems'})
   
hold off


