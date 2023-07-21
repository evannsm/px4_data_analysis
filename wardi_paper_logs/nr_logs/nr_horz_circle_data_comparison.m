clc;clear;close all


% Gather References to Compare Against
[vertical_circle, horizontal_circle, fig8_horz, fig8_vert_short, fig8_vert_tall] = get_reference_trajectories();

% Gather Logged Flight Data from .csv file
M = readmatrix('log_files/nr_horz_circ.csv');
[xs, ys, zs, yaws, thrusts, roll_rates, pitch_rates, yaw_rates, des_xs, des_ys, des_zs, des_yaws]=fix_data(M);


%% Horizontal Circle Data Comparison:
% to prove they're the same and the path is therefore independent
% of the reference being different in either program

% Vertical Circle Matlab Reference
r = horizontal_circle;

% Set Up Figure for Reference vs Path Comparisons
n=2;m=1;
figure(1)
sgtitle('Newton-Raphson Tracker - Horizontal Circle')

% Subplot 1: matlab reference vs true path
subplot(n,m,1);
hold on;
plot(r(1,:),r(2,:), '--')
plot(xs,ys);
title('Trajectory vs Matlab Reference')
xlabel('$x$','Interpreter','latex');
ylabel('$y$','Interpreter','latex');
xlim([-2,2]);
ylim([-1.5,1.5]);

% Subplot 2: gazebo reference vs true path
subplot(n,m,2);
hold on;
plot(des_xs, des_ys)
plot(xs,ys)
title('Trajectory vs Gazebo Reference')
xlabel('$x$','Interpreter','latex');
ylabel('$y$','Interpreter','latex');
xlim([-2,2]);
ylim([-1.5,1.5]);

%% X,Y,Z OVer Time
%Set Up Figure for x,y,z over time
n=4;m=1;
figure(3)
sgtitle('X,Y,Z Over time - Horizontal Circle')

% Subplot 1: xs over time
subplot(n,m,1);
plot(xs)
title('X(m) Over Time')
xlabel('$t$','Interpreter','latex');
ylabel('$X(m)$','Interpreter','latex');
ylim([-3.5,3.5])


% Subplot 2: ys over time
subplot(n,m,2);
plot(ys)
title('Y(m) Over Time')
xlabel('$t$','Interpreter','latex');
ylabel('$Y(m)$','Interpreter','latex');
ylim([-3.5,3.5])


% Subplot 3: zs over time
subplot(n,m,3);
plot(zs)
title('Z(m) Over Time')
xlabel('$t$','Interpreter','latex');
ylabel('$Z(m)$','Interpreter','latex');
ylim([0.0,3.5])
% Subplot 4: yaw over time
subplot(n,m,4);
plot(yaws)
title('Yaw(deg) Over Time')
xlabel('$t$','Interpreter','latex');
ylabel('$yaw(deg)$','Interpreter','latex');


%% Thrusts and Rates over Time

% Set Up Figure for Thrusts and Rates over Time
n=4;m=1;
figure(2)
sgtitle('Thrusts and Rates over Time - Horizontal Circle')

% Subplot 1: thrust over time
subplot(n,m,1);
plot(thrusts)
title('Thrust Over Time')
xlabel('$t$','Interpreter','latex');
ylabel('$thrust$','Interpreter','latex');
ylim([0.0,1.0]);

% Subplot 2: roll rate over time
subplot(n,m,2);
plot(roll_rates)
title('Roll Rate Over Time')
xlabel('$t$','Interpreter','latex');
ylabel('$roll rates$','Interpreter','latex');
ylim([-0.9,0.9]);


% Subplot 3: pitch rate over time
subplot(n,m,3);
plot(pitch_rates)
title('Pitch Rate Over Time')
xlabel('$t$','Interpreter','latex');
ylabel('$pitch rates$','Interpreter','latex');
ylim([-0.9,0.9]);

% Subplot 4: yaw rate over time
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