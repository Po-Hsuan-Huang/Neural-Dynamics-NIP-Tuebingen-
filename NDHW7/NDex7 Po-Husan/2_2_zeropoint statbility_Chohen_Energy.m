clear;
s =[1 2];

% points of plane1
pA{1} = [0 0 0 ];
pB{1} = [-10 -90 -10];
pC{1} = [10 90 10];
N{1} = [ 10 -1 -1]; % normal of plane1
% points of plane2

pA{2} = [0 0 0 ];
pB{2} = [10 2 10];
pC{2} = [-10 -2 -10];
N{2} = [-1 10 -1]; % normal of plane2


figure(1);


for i= 1:length(s)
pointA = pA{i};
pointB = pB{i};
pointC = pC{i}; 
    
    
normal = N{s}; %# Calculate plane normal
%# Transform points to x,y,z
x = [pointA(1) pointB(1) pointC(1)];  
y = [pointA(2) pointB(2) pointC(2)];
z = [pointA(3) pointB(3) pointC(3)];

%Find all coefficients of plane equation    
A = normal(1); B = normal(2); C = normal(3);
D = -dot(normal,pointA);
%Decide on a suitable showing range
xLim = [min(x) max(x)];
zLim = [min(z) max(z)];
[X,Z] = meshgrid(xLim,zLim);
Y = (A * X + C * Z + D)/ (-B);
reOrder = [1 2  4 3];



% plot
subplot(2,1,i)
patch(X(reOrder),Y(reOrder),Z(reOrder),'b');
grid on;
alpha(0.3);

hold on


end

hold off

U = [ 12 12 100];
K= -(1/50)*(U(1)^2+U(2)^2)+(49/50)*U(3)^2 +(2/5)*U(1)*U(2)-(83/100)*(U(1)+U(2))*U(3)


