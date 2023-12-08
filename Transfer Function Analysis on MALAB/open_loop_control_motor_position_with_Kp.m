% This script is for motor position control using PID controller

% 1. Given physical parameters

J = 4.84E-6;   % Moment of inertia (kg.m^2/s^2)
b = 5.57E-6;   % Damping ratio of the mechanical system (Nms)
Kt = 0.0374;   % Electromotive force constant (Nm/Amp)
R = 3;         % Electric resistance (Ohms)
L = 6.75E-6;   % Electric inductance (H)
% L = 0;         % system order reduction
K = Kt;        % Electromotive force constant (assumed equal to Kt)

% Define the numerator and denominator polynomials
% numerator = Kt;
% denominator = [J*L, J*R+b*L, b*R+Kt*K, 0];



% 2. Create the transfer function model 'plant'

% plant = tf(numerator, denominator);
s = tf('s');
P = K / (s*((J*s+b)*(L*s+R)+K^2));



% 3. Create system model

% Define the parameters for the proportional C
Kp = 0.05;  % Proportional gain

% Create different PID controllers

% C = pid(Kp); % only P
% C = Kp/s; % only I
% C = Kp * (s + 50) * (s + 20) / s / (s + 500);
% C = Kp * (s+60) * (s+70) / s; % PID control
C = Kp * (s^2 + 200*s + 13600) / s;

% Open loop transfer function
G = series(C, P);
% Feedback
H = 1;

% Connect the plant and C in a negative feedback configuration
sysClosedLoop = feedback(G, H);

% Take disturbance as input
dist_sys = feedback(P, C);



% 4. Plot 

% Plot the step response due to a unit step input
% figure;
subplot(2,1,1);
step_response = step(sysClosedLoop);
plot(step_response);
title('Step Response due to Unit Step Input');
xlabel('Time');
ylabel('Position');
grid on;
hold on;

% Plot the step response due to a unit step disturbance
subplot(2,1,2);
step_response = step(dist_sys);
plot(step_response);
title('Step Response due to Unit Step Disturbance');
xlabel('Time');
ylabel('Position');
grid on;
hold on;

% Plot Bode diagram
figure;
margin(sysClosedLoop); % closed loop TF margin



% 5. Disp system info

stepInfoData = stepinfo(sysClosedLoop);   % 获取阶跃响应的信息，包括过冲

overshootValue = stepInfoData.Overshoot;   % 从信息中提取过冲值
disp(['系统的过冲值: ', num2str(overshootValue), ' %']); % 显示过冲值
SettlingTime = stepInfoData.SettlingTime; % 从信息中提取settling time
disp(['系统的稳定时间: ', num2str(SettlingTime), ' s or ', num2str(1000 * SettlingTime), ' ms']); % 显示settling time

[wn,zeta] = damp(sysClosedLoop); % find closed loop damping ratio
disp(['系统的一阶阻尼比: ', num2str(zeta(1))]); % 显示阻尼比
disp(['系统的一阶自然频率: ', num2str(wn(1)), ' rad/s']); % 显示对应的自然频率

bw = bandwidth(sysClosedLoop); % open loop system bandwidth
disp(['系统的带宽(-3dB): ', num2str(bw), ' rad/s or ', num2str(bw/2/pi), 'Hz']); % 显示当前系统对应的带宽

[Gm,Pm,Wcg,Wcp] = margin(sysClosedLoop);
disp(['系统的幅值裕度: ', num2str(Gm), 'dB']); % 显示GM
disp(['系统的相位裕度: ', num2str(Pm), 'deg']); % 显示PM



% 6. Calculate Steady State Error

% 获取系统对单位阶跃输入的稳态误差
t = 0:0.01:500;
[y, ~] = step(sysClosedLoop, t);
steadyStateError = y(end) - 1;  % 期望的稳态值是1（单位阶跃输入）
disp(['系统的稳态误差: ', num2str(steadyStateError)]); % 显示结果

% 获取系统对单位阶跃和脉冲扰动的稳态误差
% [y, ~] = step(dist_sys, t);
[y_d_s, ~] = step(dist_sys, t);
steadyStateError_dist = y_d_s(end);  % 期望的稳态值是0（单位阶跃扰动）
disp(['系统对单位阶跃扰动的稳态误差: ', num2str(steadyStateError_dist)]); % 显示结果
[y_d_i, ~] = impulse(dist_sys, t);
steadyStateError_dist_i = y_d_i(end);  % 期望的稳态值是0（单位脉冲扰动）
disp(['系统对单位脉冲扰动的稳态误差: ', num2str(steadyStateError_dist_i)]); % 显示结果



% 7. 绘制根轨迹图
figure;
rlocus(sysClosedLoop);
title('Root Locus Diagram');
grid on;
