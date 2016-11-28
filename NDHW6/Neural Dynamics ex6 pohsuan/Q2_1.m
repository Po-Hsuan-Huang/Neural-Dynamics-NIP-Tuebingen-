% Neural Dynamics Exercise 6 
% Po-Hsuan Huang
% 2014/12/17
% Q2_1
% plot vector field of the matrix

%% 
close
clear

c = 2;

s = [1;1];

A = [0 1 ; 1 0];

point = -10:1:10;

thres =@(g) (g<=0)*0+(g>0)*g;

for i =1:length(point)
for j = 1:length (point)
    
    u(i,j,:)=-[point(i);point(j)]-c*A*[thres(-10+point(i));thres(-10+point(j))]+s;
    
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
title({'velocity field '})
axis tight
hold off




%% for second quadrat


point = 0:1:10;
clear u
for i = 1:length(point)
for j = 1:length (point)
    
    u(i,j,:)=-[-10+point(i);point(j)]-c*A*[thres(-10+point(i));thres(point(j))]+s;
    
end
end
[p ,q]= meshgrid(-10:1:0,0:1:10);


figure(2)
quiver(p,q,u(:,:,1),u(:,:,2)); hold on
xlabel('x1');
ylabel('x2');
title({'2nd quadrant velocity field '})
axis tight
hold off


%% for 3rd quadrant
point = 0:1:10;
clear u
for i = 1:length(point)
for j = 1:length (point)
    
    u(i,j,:)=-[-10+point(i);-10+point(j)]-c*A*[thres(-10+point(i));thres(-10+point(j))]+s;
    
end
end
[p ,q]= meshgrid(-10:1:0);


figure(3)
quiver(p,q,u(:,:,1),u(:,:,2)); 
xlabel('x1');
ylabel('x2');
title({'3rd quadrant velocity field '})
axis tight

%% for 4th quadrant
point = 0:1:10;
clear u
for i = 1:length(point)
for j = 1:length (point)
    
    u(i,j,:)=-[point(i);-10+point(j)]-c*A*[thres(point(i));thres(-10+point(j))]+s;
    
end
end
[p ,q]= meshgrid(0:1:10,-10:1:0);


figure(4)
quiver(p,q,u(:,:,1),u(:,:,2)); 
xlabel('x1');
ylabel('x2');
title({'4th quadrant velocity field '})
axis tight
%% Q2_4

% for 1st quadrant, threshold(u) = u;

% system metrix

S = -c*A-eye(2);
[V ,D]=eig(S);
display(S,'system metrix')
display(D,'dioganalized metrix')
display(V,'right eigenvector')

% for 2nd quadrant, threashold(u1)= 0

S2 = [-1 -c; 0 -1];
[V2 ,D2]=eig(S2);
display(S2,'system metrix')
display(D2,'dioganalized metrix')
display(V2,'right eigenvector')
% for 4th quadrant, threshold(u2)= 0
S3 = [-1 0; -c -1];
[V3 ,D3]=eig(S3);
display(S3,'system metrix')
display(D3,'dioganalized metrix')
display(V3,'right eigenvector')
% for 3rd quadrant, threashold(u1)= 0,threshold(u2)= 0
S4 = [-1 0; 0 -1];
[V4 ,D4]=eig(S4);
display(S4,'system metrix')
display(D4,'dioganalized metrix')
display(V4,'right eigenvector')
