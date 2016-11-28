



%% study of stationary state k-mode weight distribution with different a.
clear
figure()

a=1;

 k0=-10:1:10;
 k = -10:0.1:10;
d =1;
for j = 1:length(k0)
for i = 1:length(k)
func2(j,i) = exp(-d^2 * k(i)^2) / (1-a*(exp(-0.6^2*(k(i)-k0(j))^2)+exp(-0.6^2*(k0(j)+k(i))^2)));

end
end
sf =surf(k,k0,func2);
set(sf,'LineStyle','none');
colormap summer;
camlight('headlight')
title('k0 versus k-th mode peak value');
ylabel('k th mode');
xlabel('ko')
figure(4)
plot(k,func2);
title('a versus k-th mode peak value');
xlabel('k th mode');
ylabel('value')


%%
figure(2)
subplot(2,1,1)
a = 1;
k=-20:0.1:20;
k0 =8;
func3 = (1-a*(exp(-0.6^2*(k-k0).^2)+exp(-0.6^2*(k0+k).^2))).^-1;  % denominator
func4 =  exp(-d^2 * k.^2).*cos(k*1) .* func3;  % whole factor
%plot(k, func3,k,func4);
plot(k,func4);

title('a versus k-th mode peak value');
xlabel('k th mode');
ylabel('a')
legend('k0 =8')



subplot(2,1,2)

k=-20:0.1:20;
k0 =4;
func3 = (1-a*(exp(-0.6^2*(k-k0).^2)+exp(-0.6^2*(k0+k).^2))).^-1;
func4 =  exp(-d^2 * k.^2).*cos(k*1).*func3;
func5 = func4 .* real(exp(1i*k*9))
%plot(k, func3,k,func5)
plot(k,func5);

title('a versus k-th mode peak value');
xlabel('k th mode');
ylabel('a')
legend('k0 =4')
