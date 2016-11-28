% Neural Dynamics Exercise 6 
% Po-Hsuan Huang
% 2014/12/6
% Q1_3
% plot vector field of the matrix


%%

%% 
clear
close all

%% constants

A = [-0.5 -0.5 0; -0.5 -0.5 0; 0 0 2];  %  dynamics metrix for dx/dt = Ax
T_t = 10 ; %   total modeling time.
dt = 1e-1 ; % time stepsize
[V ,D, W] = eig(A); % calculate the right eigenvectors and Diagonalized matrix.

display (W,'right eigenvectors')

display (D,'diagonalized metrix')
%% Q1_3
% easy way.


% [x,y]=meshgrid(-10:2:10);
% 
% point= -10:2:10;
% 
% for j = 1: length(point)
% for i = 1: length(point)
% z(i,j,:) = A*[point(i);point(j);0];
% end
% end
% 
% figure(2)
% quiver(x,y,z(:,:,1),z(:,:,2));% proection on x3 =0
% xlabel('x1');
% ylabel('x2');
% title({'velocity field x3 = 0'})
% axis tight
% 
% [x,y]=meshgrid(-10:2:10);
% 
% clear z
% for i = 1: length(point)
% for k = 1: length(point)
% z(i,k,:) = A*[point(i);0 ;point(k)];
% end
% end
% figure(3)
% quiver(x,y,z(:,:,3),z(:,:,1)); % proection on x3 =0
% xlabel('x3');
% ylabel('x1');
% title({'velocity field'})
% axis tight
% 
% 
% clear z
% for j = 1: length(point)
% for k = 1: length(point)
% z(j,k,:) = A*[0;point(j);point(k)];
% end
% end
% figure(4)
% quiver(x,y,z(:,:,3),z(:,:,2)); % proection on x3 =0
% xlabel('x3');
% ylabel('x2');
% title({'velocity field'})
% axis tight



%% faster way to plot, but didn't figure out how to do it correctly
% fix on axis, and plot the rest two.
point= -10:2:10;
clear z
for j = 1: length(point)
for i = 1: length(point)
for k = 1: length(point)
z(i,j,k,:) = A*[point(i);point(j);point(k)];

u(i,j,k,:)= D*[point(i) ;point(j); point(k)]; % y is the vector field is eigen space.
end
end
end
[x,y,w]=meshgrid(-10:2:10);

figure(1)
quiver3(x,y,w,z(:,:,:,1),z(:,:,:,2),z(:,:,:,3)); % proection on x3 =0
xlabel('x1');
ylabel('x2');
zlabel('x3');
title({'velocity field'})
axis tight

figure(2)
quiver(x(:,:,6),y(:,:,6),z(:,:,6,1),z(:,:,6,2)); % proection on x2 =0
xlabel('x1');
ylabel('x2');
title({'velocity field for x3=0'})
axis tight

figure(3)
z2 = permute(z,[1 3 2 4]);
quiver(x(:,:,6),y(:,:,6),z2(:,:,6,3),z2(:,:,6,1)); % proection on x2 =0
xlabel('x3');
ylabel('x1');
title({'velocity field for x2=0'})
axis tight

figure(4)
z3= permute(z,[2,3,1,4]);
quiver(x(:,:,6),y(:,:,6),z3(:,:,6,3),z3(:,:,6,2)); % proection on x1 =0
xlabel('x3');
ylabel('x2');
title({'velocity field for x1=0'})
axis tight

%% Q1_4 projection on eigen space.

figure(5)
quiver3(x,y,w,u(:,:,:,1),u(:,:,:,2),u(:,:,:,3)); % proection on x3 =0
xlabel('v1');
ylabel('v2');
zlabel('v3');
title({'velocity field'})
axis tight

figure(6)
quiver(x(:,:,6),y(:,:,6),u(:,:,6,2),u(:,:,6,1));% proection on x3 =0
xlabel('v2');
ylabel('v1');
title({'velocity field e1=[0.7 0.7 0], e2=[0.7 -0.7 0]'})
axis tight

figure(7)
u2 = permute(u,[1 3 2 4]);
quiver(x(:,:,6),y(:,:,6),u2(:,:,6,3),u2(:,:,6,1)); % proection on x2 =0
xlabel('v3');
ylabel('v1');
title({'velocity field for e1=[0.7 0.7 0], e3=[0 0 1]'})
axis tight

figure(8)
u3= permute(u,[2,3,1,4]);
quiver(x(:,:,6),y(:,:,6),u3(:,:,6,3),u3(:,:,6,2)); % proection on x1 =0
xlabel('v3');
ylabel('v2');
title({'velocity field for e2=[0.7 -0.7 0], e3=[0 0 1]'})
axis tight










