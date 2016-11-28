%% Neural Dynamics Ex8
% Q 1_4. Simulate the stationary solution of Neural Field.
% in version 2, instead of integrate from -k0 to k0, we integrate from -100
% to 100
%% clear mess 
clear all;  % close all figures


a = 1; 
b=0.6;
d=2;
k0=4;
tau = 10;
sd = 0.01; % standard deviation of Gaussian noise.
inter = [-10 10];
num = 101;           % number of intervals between neurons

n = num-1;
h = (inter(2)-inter(1))/n;  % interval between neurons

funW = @(x) (a/(sqrt(pi)*b))*exp(-x^2 /(4*b^2))*cos(k0*x);
funS = @(x)  1/(2*sqrt(pi)*d)*exp(-x^2/(4*d^2));

dt = 0.1e-3 ; % time step size
T_t = 80e-3;     % total modeling time 
Num_t= 0: dt :T_t;   % number of time steps 

%% assigning initial values at stationary point.
% U= zeros(num,1);   % initially, all neurons has 0 potential in the first timestep.
for j= 1 : num
   
    a_po = inter(1)+ (j-1)*h;       %position of the neuron 
  
    k = -100:1:100;
    
    
    % wrong way to integrate over field input.
    func3 = (1-a*(exp(-0.6^2*(k-k0).^2)+exp(-0.6^2*(k0+k).^2))).^-1;
    func4 =  exp(-d^2 * k.^2).*func3;
    summ =  sum (func4.*exp(1i*k*a_po));
    V(j,1)= real(summ);

     summ2 =0;
for i = 1:length(k)     
    alpha = exp(-b^2*(k(i)+k0)^2);
    beta = exp(-b^2*(k(i)-k0)^2);
    denom =(1-a*(beta+alpha))^-1;
    summ2 =  summ2 + exp(-d^2*k(i)^2)*exp(1i*k(i)*a_po)*denom;
    count(j,i)=summ2;
end
    U(j,1)= real(summ2);

end

%%

for step = 2: length(Num_t)  % time step iteration
    
    for j = 1: n+1         % number of neuron in a layer.
     a_po = inter(1)+ (j-1)*h;       %position of the neuron 
    
     
     
    Int = 0;    
    for k = 1:n    
    b_po = inter(1)+ (k-1)*h   ;   % position of presynapic neuron 
    Int = Int + funW(b_po-a_po)*U(k,step-1)*h; 
    end
        
    S = funS(a_po) ;% sd*randn(1,1);
    
    U(j,step) = U(step-1) + dt*(-U(j,step-1) +Int + S)/tau;  % Backward-Euler.
    
  
    end
end

%%




figure(1)

Num_x = -10:h:10;
sf = surf(Num_t,Num_x,U);   %(X(j), Y(i), Z(i,j))s
set(sf,'LineStyle','none');
colormap summer;
colorbar
axis tight;
%view(90,0);camlight;
title({'Multicompartment model' ,' neural field with stationary signal '});
ylabel('X');
xlabel('time (s)');
%%
Num_x = -10:h:10;


figure(2)
subplot(2,1,1)
plot(Num_x,U(:,1));
subplot(2,1,2)
plot(Num_x,V(:,1));


title({'initial value of  stationary field '});
ylabel('field value (U)');
xlabel('neurons');