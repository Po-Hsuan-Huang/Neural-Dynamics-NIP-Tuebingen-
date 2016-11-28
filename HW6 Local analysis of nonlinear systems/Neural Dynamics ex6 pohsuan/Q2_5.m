% Neural Dynamics Exercise 6 
% Po-Hsuan Huang
% 2014/12/17
% Q2_5
% plot vector field of the matrix

%% 
close
clear

c = 2;

s = [1.2;1];

A = [0 1 ; 1 0];

point = -10:1:10;

thres =@(g) (g<=0)*0+(g>0)*g;

for i =1:length(point)
for j = 1:length (point)
    
    u(i,j,:)=-[point(i);point(j)]-c*A*[thres(point(i));thres(point(j))]+s;
    
end
end


S3= u(:,:,1)+u(:,:,2)==0;
stable = find(S3); % index of non-zero element, which is the stable point
% r =mod(stable,length(point));
% c =(stable-r)/length(point); 

[x ,y ]= meshgrid(-10:1:10);

row = S3.*x;
col = S3.*y;

figure(1)
quiver(x,y,u(:,:,1),u(:,:,2)); hold on
%plot(row,col,'ro');
xlabel('x1');
ylabel('x2');
title({'velocity field s = [1.2 1] '})
axis tight
hold off

%% Tragectory

%% Dynamical system solution

dt = 1e-2;
T_t = 2;
x0 = [ 1, -1 ;-1, 1; 0, 0 ];
time = 0:dt:T_t;


figure(3)
for  trail = 1:3
    
    
w(:,1)= x0(trail,:)';

for iter =1 : length(time)-1
    w(:,iter+1) = w(:,iter)+(-[w(1,iter);w(2,iter)]-c*A*w(:,iter)+s)*dt;  
    % check if x= inv(V)*y transform can produce correct trajactory in realspace 
    
end
i = 1:1:length(time);
% transform y back to real space.

scatter(w(1,:),w(2,:));  % trajactory in real space

hold on
end
title('tragectory with initial points x0 = [0, 0 ] [1, -1] [-1, 1 ], s= [1.2 1]')
xlabel('x1');
ylabel('x2');
grid on
hold off
%% initial condition 2


clear u

c = 2;

s = [1;1.2];

A = [0 1 ; 1 0];

point = -10:1:10;

thres =@(g) (g<=0)*0+(g>0)*g;

for i =1:length(point)
for j = 1:length (point)
    
    u(i,j,:)=-[point(i);point(j)]-c*A*[thres(point(i));thres(point(j))]+s;
    
end
end


S3= u(:,:,1)+u(:,:,2)==0;
stable = find(S3); % index of non-zero element, which is the stable point
% r =mod(stable,length(point));
% c =(stable-r)/length(point); 

[x ,y ]= meshgrid(-10:1:10);

row = S3.*x;
col = S3.*y;

figure(2)
quiver(x,y,u(:,:,1),u(:,:,2));
%plot(row,col,'ro');
xlabel('x1');
ylabel('x2');
title({'velocity field s = [1 1.2] '})
axis tight
hold off

figure(4)
for  trail = 1:3
    
    
w(:,1)= x0(trail,:)';

for iter =1 : length(time)-1
    w(:,iter+1) = w(:,iter)+(-[w(1,iter);w(2,iter)]-c*A*w(:,iter)+s)*dt;  
    % check if x= inv(V)*y transform can produce correct trajactory in realspace 
    
end
i = 1:1:length(time);
% transform y back to real space.

scatter(w(1,:),w(2,:));  % trajactory in real space

hold on
end
title('tragectory with initial points x0 = [0, 0 ] [1, -1] [-1, 1 ], s= [1 1.2]')
xlabel('x1');
ylabel('x2');
grid on
hold off