clc;clear;close all


% Gather References to Compare Against
[vertical_circle, horizontal_circle, fig8_horz, fig8_vert_short, fig8_vert_tall] = get_reference_trajectories();

% Gather Logged Flight Data from .csv file
M = readmatrix('log_files/nr_vert_circ.csv');
[xs, ys, zs, yaws, thrusts, roll_rates, pitch_rates, yaw_rates, des_xs, des_ys, des_zs, des_yaws]=fix_data(M);


%% Vertical Circle Data Comparison:
% to prove they're the same and the path is therefore independent
% of the reference being different in either program

% Vertical Circle Matlab Reference
r = vertical_circle;

% Set Up Figure
n=2;m=1;
figure(1)
sgtitle('Newton-Raphson Tracker Vertical Circle')


subplot(n,m,1);
hold on;
plot(r(1,:),r(3,:), '--')
M = readmatrix('log_files/nr_vert_circ.csv');
[xs, ys, zs, yaws, thrusts, roll_rates, pitch_rates, yaw_rates, des_xs, des_ys, des_zs, des_yaws] = fix_data(M);
plot(xs,zs);
title('Trajectory vs Matlab Reference')
xlabel('$x$','Interpreter','latex');
ylabel('$z$','Interpreter','latex');
xlim([-2,2]);
ylim([0,3]);

subplot(n,m,2);
hold on;
plot(des_xs, des_zs)
plot(xs,zs)
title('Trajectory vs Gazebo Reference')
xlim([-2,2]);
ylim([0,3]);


% Thrusts and Rates over Time
n=4;m=1;
% n=2;m=2;

figure(2)
sgtitle('Thrusts and Rates over Time - Vertical Circle')
subplot(n,m,1);
plot(thrusts)
title('Thrust Over Time')
xlabel('$t$','Interpreter','latex');
ylabel('$thrust$','Interpreter','latex');
ylim([0.0,1.0]);

subplot(n,m,2);
plot(roll_rates)
title('Roll Rate Over Time')
xlabel('$t$','Interpreter','latex');
ylabel('$roll rates$','Interpreter','latex');
ylim([-0.9,0.9]);

subplot(n,m,3);
plot(pitch_rates)
title('Pitch Rate Over Time')
xlabel('$t$','Interpreter','latex');
ylabel('$pitch rates$','Interpreter','latex');
ylim([-0.9,0.9]);

subplot(n,m,4);
plot(yaw_rates)
title('Yaw Rate Over Time')
xlabel('$t$','Interpreter','latex');
ylabel('$roll rates$','Interpreter','latex');
ylim([-0.9,0.9]);


function [xs, ys, zs, yaws, thrusts, roll_rates, pitch_rates, yaw_rates, des_xs, des_ys, des_zs, des_yaws]=fix_data(M)
    MM = M(:,2);
    ind = ~isnan(MM);
    MMM=MM(ind);
    MMM(end+1) = 0;
    num_rows = 13;
    numCols = length(MMM) / num_rows;
    tempMatrix = reshape(MMM, num_rows, numCols);
    tempMatrix(num_rows, :) = [];
    resultMatrix = tempMatrix;
    
    
    xs=resultMatrix(1,:);
    ys=resultMatrix(2,:);
    zs=resultMatrix(3,:);
    yaws=resultMatrix(4,:);
    thrusts=resultMatrix(5,:);
    roll_rates=resultMatrix(6,:);
    pitch_rates=resultMatrix(7,:);
    yaw_rates=resultMatrix(8,:);
    des_xs=resultMatrix(9,:);
    des_ys=resultMatrix(10,:);
    des_zs=resultMatrix(11,:);
    des_yaws=resultMatrix(12,:);

    
end



function [vertical_circle, horizontal_circle, fig8_horz, fig8_vert_short, fig8_vert_tall] = get_reference_trajectories()

    t = 0:.02:20;
    w=1;
    vertical_circle = [cos(w.*t); zeros(1,length(t)); sin(w.*t)+2; (pi/2)*ones(1,length(t))];
    horizontal_circle = [cos(w.*t); sin(w*t); 3*ones(1,length(t)); (pi/2)*ones(1,length(t))];
    fig8_horz = [sin(t./2); sin(2*t/2); 3*ones(1,length(t)); (pi/2)*ones(1,length(t))];
    fig8_vert_short = [sin(t./2); zeros(1,length(t)); sin(2*t/2)+3*ones(1,length(t)); (pi/2)*ones(1,length(t))];
    fig8_vert_tall = [sin(2*t/2); zeros(1,length(t)); sin(t./2)+2*ones(1,length(t)); (pi/2)*ones(1,length(t))];

end
