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


clear;

signal = [0, 3/4, 3/4];

fix = [ 0 0 ; 1/2 1/2 ;-3/2 -3/2] ; % fixed points for two cases.

for s = 1: length(signal)
    
range= {  -0.5:0.1:0.5  ; 0:0.1:1 ;  -1.8: 0.05 :-1.3  }; 

point = range{s};

[x,y]=meshgrid(point);
 
    
for i = 1: length(point)
for j = 1: length(point)
 

  
M = [ -1  -signal(s)*(1+point(j))^-2;  -signal(s)*(1+point(i))^-2 -1 ];  %  dynamics metrix for du/dt = -u + [Mu]




zero = [fix(s,1)-signal(s)/(1+fix(s,2)); fix(s,2)-signal(s)/(1+fix(s,1))];

% zero = f0 = 0 at fixed points.

% linearization by Taylor expansion to 1st order (only accurate near fixed point)

% du/dt = f0 + df/du * ( x- x0);

z(j,i,:) = zero +M*([point(i)-fix(s,1);point(j)-fix(s,2)]);

if point(i)==-1 || point(j)==-1    
   
    z(j,i,:) = [0;0];
    
end  
    


end
end


figure(s)
quiver(x,y,z(:,:,1),z(:,:,2));% proection on x3 =0
ax = gca;
ax.XTick = -2:0.1:1;
ax.YTick = -2:0.1:1;
grid(ax)

xlabel('x1');
ylabel('x2');

str = sprintf(' linearized system input field s = %d,  ', signal(s));
str2 = sprintf('fixed point = %d, %d  ', fix(s,1),fix(s,2));
title  ({str, str2 });

axis tight

end

