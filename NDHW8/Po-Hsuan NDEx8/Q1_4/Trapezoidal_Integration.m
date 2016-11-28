%% Name:
% Numerical Method Ex7 2014/12/16
% Trapeziodal Rule Integration

% First-order Newton-Cotes. 
% with two point form f1+ f2 in each panel.



%% input requirement:

% Integrating functon, argument func should be a function handle, inter is the
% inteval of integration in form of [lowerlimit  upperlimit], and n is the number fo 
% points in the interval.
% Minimum 3 data points are required .

%% limitation:
% This  method has limitation that inside integration interval the function
% should not oscillate, and its 2th derivatie had better be constant.
 

function [Integral, Error, h]=Trapezoidal_Integration(func, n, inter)

a = inter(1);
b=  inter(2);

h =( b-a)/(n-1);  % stepsize for equally spaced n points


%Use trapazoidal rule to calculate the integral.

Integral = 0;

Error= 0;

if   mod(n,1)== 0 && n>= 3 % n odd

x1 = 1:1:n-2; % construct an array of inteval after point a

sum1 = sum(2*func(a + h*x1));

x2 = [0 n-1];

sum2 = sum(func(a + h*x2));

Integral =  (h/2)*(sum1+sum2);

else 
    
    display(' n is not a integer >= 3' );
end
       
%%   calculate Error

x0 =0:1:n-1;

array = func(a+h*x0);  

dev = diff(array,2) ;
% 2th differentaion in Newton-Cotes method. Use arbituary element of dev,
% since we assume they should be nearly the same if there is at most only one 
% extremum in the inteval of integration .
middle = dev(round(length(dev)/2))/ h^2;

Error = -(b-a)*h^3*middle/12;

