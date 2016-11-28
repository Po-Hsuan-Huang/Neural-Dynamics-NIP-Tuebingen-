% p_HW8_Neural Dynamics
%
% Santiago Cadena 
% sa.cadnea721@gmail.com
%
clear; close all;clc
%% 

 % Parameter:
deltaT = 0.5;
T = 200;
A = -10;
B = 10;
tau= 10;
a = 1;
b = 0.6;
d = 2;
k0 = 4;
varS = 0.01^2;
eta = deltaT/tau;
M = 200;
f_s = @(x)(1/(2*sqrt(pi)*d))*exp(-(x.^2)/(4*d^2))';
f_w = @(x) ((a/(sqrt(pi)*b))*(exp(-(x.^2)/(4*b^2)).*cos(k0*x)))';
deltax = (B-A)/M;
x = A:deltax:B;
%%
% Spatial Fourier transform
k = -2:0.001:2;
v_a = 0.5:0.05:1;
u_k = zeros(length(a),length(k));
for i = 1:length(v_a) 
    u_k(i,:) = (exp(-(d^2)*k.^2))./(1-v_a(i)*(exp(-(b^2)*(k-k0))+exp(-(b^2)*(k+k0))));
end

figure,subplot(1,2,1);imagesc(k,v_a,u_k);
xlabel('$k$','interpreter','latex','Fontsize',18);
ylabel('$a$','interpreter','latex','Fontsize',18);
title('$\tilde{u}(k)$','interpreter','latex','Fontsize',18);
subplot(1,2,2);plot(k,u_k(end,:))
xlabel('$k$','interpreter','latex','Fontsize',18);
ylabel('$\tilde{u}(k)$','interpreter','latex','Fontsize',18);
title('$a=1$','interpreter','latex','Fontsize',18);
%%

% Euler approximation
time = 0:deltaT:T;
u = zeros(length(time),M+1);
z = zeros(length(time),M+1);
for tind = 1:length(time)
    for xind = 1:length(x)  
        z(tind,xind) = sum(f_w(x(xind)-x).*u(tind,:)'*deltax);
    end
    u(tind+1,:) = u(tind,:) + eta*(-u(tind,:)+z(tind,:)+(f_s(x)'));
end

figure,subplot(1,2,1);imagesc(time,x,u')
xlabel('time','FontSize',16)
ylabel('$x$','interpreter','latex','FontSize',20)
title('No noise','FontSize',14)

% Plott with noise:
un = zeros(length(time),M+1);
z = zeros(length(time),M+1);
for tind = 1:length(time)
    for xind = 1:length(x)  
        z(tind,xind) = sum(f_w(x(xind)-x).*un(tind,:)'*deltax);
    end
    un(tind+1,:) = un(tind,:) + eta*(-un(tind,:)+z(tind,:)+(f_s(x)'+sqrt(varS)*randn(size(f_s(x)))'));
end

subplot(1,2,2);imagesc(time,x,un')
xlabel('time','FontSize',16)
ylabel('$x$','interpreter','latex','FontSize',20)
title('With Gaussian noise','FontSize',14)

% Plott with a = 0.7 and 1.5
figure,
a = 0.7;
f_w = @(x) ((a/(sqrt(pi)*b))*(exp(-(x.^2)/(4*b^2)).*cos(k0*x)))';
u = zeros(length(time),M+1);
z = zeros(length(time),M+1);
for tind = 1:length(time)
    for xind = 1:length(x)
        z(tind,xind) = sum(f_w(x(xind)-x).*u(tind,:)'*deltax);
    end
    u(tind+1,:) = u(tind,:) + eta*(-u(tind,:)+z(tind,:)+(f_s(x)'));
end
subplot(1,2,1)
imagesc(time,x,u')
xlabel('time','FontSize',16)
ylabel('$x$','interpreter','latex','FontSize',20)
title('$a = 0.7$','interpreter','latex','FontSize',16)
%%
a = 1.5;
time = 0:deltaT:T;

f_w = @(x) ((a/(sqrt(pi)*b))*exp(-(x.^2)/(4*b^2)).*cos(k0*x))';
u = zeros(length(time),M+1);
z = zeros(length(time),M+1);
for tind = 1:length(time)
    for xind = 1:length(x)
        z(tind,xind) = sum(f_w(x(xind)-x).*u(tind,:)'*deltax);
    end
    u(tind+1,:) = u(tind,:) + eta*(-u(tind,:)+z(tind,:)+(f_s(x)'));
end
subplot(1,2,2)
imagesc(time,x,u')
xlabel('time','FontSize',16)
title('$a = 1.5$','interpreter','latex','FontSize',16)
%%
% When k0 = 8
a =1.5;
k0 = 8;
f_w = @(x) ((a/(sqrt(pi)*b))*(exp(-(x.^2)/(4*b^2)).*cos(k0*x)))';
u = zeros(length(time),M+1);
z = zeros(length(time),M+1);
for tind = 1:length(time)
    for xind = 1:length(x)  
        z(tind,xind) = sum(f_w(x(xind)-x).*u(tind,:)'*deltax);
    end
    u(tind+1,:) = u(tind,:) + eta*(-u(tind,:)+z(tind,:)+(f_s(x)'));
end
figure,
imagesc(time,x,u')
xlabel('time','FontSize',16)
title('$k_0 = 8$','interpreter','latex','FontSize',16)

%% part 1.6

% Parameter:
deltaT = 0.5;
T = 200;
A = -10;
B = 10;
tau= 10;
c =1;
d1 = 0.5;
v = 0.1;
eta = deltaT/tau;
M = 200;
f_s = @(x,t) (c/(2*sqrt(pi)*d1))*exp(-((x-v*t).^2)/(4*d1^2))';
f_w = @(x) sign(x).*exp(-c*abs(x));
deltax = (B-A)/M;
x = A:deltax:B;

% Euler approximation
time = 0:deltaT:T;
u = zeros(length(time),M+1);
z = zeros(length(time),M+1);
for tind = 1:length(time)
    for xind = 1:length(x)  
        z(tind,xind) = sum(f_w(x(xind)-x).*u(tind,:)*deltax);
    end
    u(tind+1,:) = u(tind,:) + eta*(-u(tind,:)+z(tind,:)+(f_s(x,time(tind))'));
end
figure,subplot(1,2,1)
imagesc(time,x,u')
xlabel('time','FontSize',16)
ylabel('$x$','interpreter','latex','FontSize',20)
title('$v = 0.1$','interpreter','latex','FontSize',16)

% When v = -0.1
v= -0.1;
f_s = @(x,t) (c/(2*sqrt(pi)*d1))*exp(-((x-v*t).^2)/(4*d1^2))';
u = zeros(length(time),M+1);
z = zeros(length(time),M+1);
for tind = 1:length(time)
    for xind = 1:length(x)  
        z(tind,xind) = sum(f_w(x(xind)-x).*u(tind,:)*deltax);
    end
    u(tind+1,:) = u(tind,:) + eta*(-u(tind,:)+z(tind,:)+(f_s(x,time(tind))'));
end
subplot(1,2,2)
imagesc(time,x,u')
xlabel('time','FontSize',16)
ylabel('$x$','interpreter','latex','FontSize',20)
title('$v = -0.1$','interpreter','latex','FontSize',16)

%% Part 2.1
A = 3;
B = 2;
C = 0.6;
a = 1;
b = 3;
d = 4;
h =1;
tau = 10;
deltaT = 0.5;
T = 200;
eta = deltaT/tau;
M = 200;
Xlow= -10;
Xhigh = 10;
deltax = (Xhigh-Xlow)/M;
x = Xlow:deltax:Xhigh;
f_s = @(x)(C*(1-abs(x)/d)).*(abs(x)<d);
f_w = @(x) (A*(abs(x)<=a))+(-B*((abs(x)>a)&(abs(x)<=b)));

% Euler approximation
time = 0:deltaT:T;
u = zeros(length(time),M+1);
z = zeros(length(time),M+1);
u(1,:) = -h*ones(1,M+1);
for tind = 1:length(time)
    for xind = 1:length(x)  
        z(tind,xind) = sum(f_w(x(xind)-x).*(u(tind,:)>0)*deltax);
    end
    u(tind+1,:) = u(tind,:) + eta*(-u(tind,:)+z(tind,:)+(f_s(x))-h);
end
figure,subplot(1,2,1)
imagesc(time,x,u')
xlabel('time','FontSize',16)
ylabel('$x$','interpreter','latex','FontSize',20)
title('$u(x,0) = -h$','interpreter','latex','FontSize',16)

% Initial condition: u(x,0) = 3.3s(x)-h

% Euler approximation
time = 0:deltaT:T;
u2 = zeros(length(time),M+1);
z = zeros(length(time),M+1);
u2(1,:) = 3.3*f_s(x)-h;
for tind = 1:length(time)
    for xind = 1:length(x)  
        z(tind,xind) = sum(f_w(x(xind)-x).*(u2(tind,:)>0)*deltax);
    end
    u2(tind+1,:) = u2(tind,:) + eta*(-u2(tind,:)+z(tind,:)+(f_s(x))-h);
end
subplot(1,2,2)
imagesc(time,x,u2')
xlabel('time','FontSize',16)
ylabel('$x$','interpreter','latex','FontSize',20)
title('$u(x,0) = 3.3s(x)-h$','interpreter','latex','FontSize',16)


% Stable solutions
figure,subplot(1,2,1)
plot(x,u(end,:),'LineWidth',2)
ylabel('$u^{*}(x)$','interpreter','latex','FontSize',20)
xlabel('$x$','interpreter','latex','FontSize',20)
title('$u(x,0) = -h$','interpreter','latex','FontSize',16)
subplot(1,2,2)
plot(x,u2(end,:),'LineWidth',2)
ylabel('$u^{*}(x)$','interpreter','latex','FontSize',20)
xlabel('$x$','interpreter','latex','FontSize',20)
title('$u(x,0) = 3.3s(x)-h$','interpreter','latex','FontSize',16)

% For w = 0.1w

f_w = @(x) 0.1*(A*(abs(x)<=a))+(-B*((abs(x)>a)&(abs(x)<=b)));

% Euler approximation
time = 0:deltaT:T;
u = zeros(length(time),M+1);
z = zeros(length(time),M+1);
u(1,:) = -h*ones(1,M+1);
for tind = 1:length(time)
    for xind = 1:length(x)  
        z(tind,xind) = sum(f_w(x(xind)-x).*(u(tind,:)>0)*deltax);
    end
    u(tind+1,:) = u(tind,:) + eta*(-u(tind,:)+z(tind,:)+(f_s(x))-h);
end
figure,subplot(1,2,1)
imagesc(time,x,u')
ylim([-1,-0,4])
xlabel('time','FontSize',16)
ylabel('$x$','interpreter','latex','FontSize',20)
title('$u(x,0) = -h$','interpreter','latex','FontSize',16)

% Initial condition: u(x,0) = 3.3s(x)-h

% Euler approximation
time = 0:deltaT:T;
u2 = zeros(length(time),M+1);
z = zeros(length(time),M+1);
u2(1,:) = 3.3*f_s(x)-h;
for tind = 1:length(time)
    for xind = 1:length(x)  
        z(tind,xind) = sum(f_w(x(xind)-x).*(u2(tind,:)>0)*deltax);
    end
    u2(tind+1,:) = u2(tind,:) + eta*(-u2(tind,:)+z(tind,:)+(f_s(x))-h);
end
subplot(1,2,2)
imagesc(time,x,u2')
xlabel('time','FontSize',16)
ylabel('$x$','interpreter','latex','FontSize',20)
title('$u(x,0) = 3.3s(x)-h$','interpreter','latex','FontSize',16)


% Stable solutions
figure,subplot(1,2,1)
plot(x,u(end,:),'LineWidth',2)
ylabel('$u^{*}(x)$','interpreter','latex','FontSize',20)
xlabel('$x$','interpreter','latex','FontSize',20)
title('$u(x,0) = -h$','interpreter','latex','FontSize',16)
subplot(1,2,2)
plot(x,u2(end,:),'LineWidth',2)
ylabel('$u^{*}(x)$','interpreter','latex','FontSize',20)
xlabel('$x$','interpreter','latex','FontSize',20)
title('$u(x,0) = 3.3s(x)-h$','interpreter','latex','FontSize',16)

% Necessary value of C > 1 for active regions
C =1.001;
f_s = @(x)(C*(1-abs(x)/d)).*(abs(x)<d);
f_w = @(x) (A*(abs(x)<=a))+(-B*((abs(x)>a)&(abs(x)<=b)));
% Euler approximation
time = 0:deltaT:T;
u = zeros(length(time),M+1);
z = zeros(length(time),M+1);
u(1,:) = -h*ones(1,M+1);
for tind = 1:length(time)
    for xind = 1:length(x)  
        z(tind,xind) = sum(f_w(x(xind)-x).*(u(tind,:)>0)*deltax);
    end
    u(tind+1,:) = u(tind,:) + eta*(-u(tind,:)+z(tind,:)+(f_s(x))-h);
end
figure,
plot(x,u(end,:),'LineWidth',2)
ylabel('$u^{*}(x)$','interpreter','latex','FontSize',20)
xlabel('$x$','interpreter','latex','FontSize',20)
title('$u(x,0) = -h$','interpreter','latex','FontSize',16)

%% part 2.4
A = 3;
B = 2;
C = 0.6;
a = 1;
b = 3;
d = 4;
h =1;
tau = 10;
deltaT = 0.5;
T = 200;
eta = deltaT/tau;
M = 200;
Xlow= -10;
Xhigh = 10;
deltax = (Xhigh-Xlow)/M;
x = Xlow:deltax:Xhigh;
f_s = @(x)(C*(1-abs(x)/d)).*(abs(x)<d);
f_w = @(x) (A*(abs(x)<=a))+(-B*((abs(x)>a)&(abs(x)<=b)));

time = 0:deltaT:T;
u2 = zeros(length(time),M+1);
z = zeros(length(time),M+1);
u2(1,:) = 3.3*f_s(x)-h;
for tind = 1:length(time)
    for xind = 1:length(x)  
        z(tind,xind) = sum(f_w(x(xind)-x).*(u2(tind,:)>0)*deltax);
    end
    u2(tind+1,:) = u2(tind,:) + eta*(-u2(tind,:)+z(tind,:)+(0)-h); %s =0
end
figure,subplot(1,2,1)
imagesc(time,x,u2')
xlabel('time','FontSize',16)
ylabel('$x$','interpreter','latex','FontSize',20)
title('$u(x,0) = 3.3s(x)-h$ but s = 0','interpreter','latex','FontSize',16)
subplot(1,2,2)
plot(x,u2(end,:),'LineWidth',2)
ylabel('$u^{*}(x)$','interpreter','latex','FontSize',20)
xlabel('$x$','interpreter','latex','FontSize',20)
title('$u(x,0) = 3.3s(x)-h$ but s = 0','interpreter','latex','FontSize',16)

