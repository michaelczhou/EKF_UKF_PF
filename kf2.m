clc
clear
close all
 
% 初始化参数
delta_t=0.1;
t=0:delta_t:5;
g=10;%加速度值
n_iter = length(t); % 序列的长度
sz = [n_iter, 1]; % 信号需开辟的内存空间大小
x=1/2*g*t.^2;
x=x';
z = x + sqrt(10).*randn(sz); % 测量时加入测量白噪声
Q = 0.9; % 过程激励噪声方差   
         %注意Q值得改变  待会增大到2，看看效果。对比看效果时，修改代码不要改变z的值
R = 10; % 测量方差估计，可以改变它来看不同效果
 
% 分配空间
xhat=zeros(sz);      % x的后验估计
P=zeros(sz);         % 后验方差估计
xhatminus=zeros(sz); % x的先验估计
Pminus=zeros(sz);    % 先验方差估计
K=zeros(sz);         % Kalman增益
 
% 估计的初始值
xhat(1) = 0.0;
P = 1.0;
for k = 2:n_iter   %
    % 时间更新过程
    xhatminus(k) = xhat(k-1);
    Pminus(k) = P(k-1)+Q;
    
    % 测量更新过程
    K(k) = Pminus(k)/( Pminus(k)+R );
    xhat(k) = xhatminus(k)+K(k)*(z(k)-xhatminus(k));
    P(k) = (1-K(k))*Pminus(k);
end
 
figure
plot(t,z);
hold on
plot(t,xhat,'r-')
plot(t,x,'g-');
legend('含有噪声的测量', '后验估计', '真值');
xlabel('Iteration');
