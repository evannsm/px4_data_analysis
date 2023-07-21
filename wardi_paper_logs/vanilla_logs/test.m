clc;clear;close all


% Gather References to Compare Against
[vertical_circle, horizontal_circle, fig8_horz, fig8_vert_short, fig8_vert_tall] = get_reference_trajectories();

% Gather Logged Flight Data from .csv file
M = readmatrix('log_files/vanilla_vert_tall_fig8.csv');
[xs, ys, zs, yaws, des_xs, des_ys, des_zs, des_yaws]=fix_data(M);


%% Vertical Tall Figure8 Data Comparison:
% to prove they're the same and the path is therefore independent
% of the reference being different in either program

% Horizontal Figure 8 Matlab Reference
r = fig8_vert_tall;

% Set Up Figure
n=2;m=1;
figure()
sgtitle('PX4 PID Reference Tracking - Vertical Tall Fig8')


% Subplot 1: matlab reference vs true path
subplot(n,m,1);
hold on;
plot(r(1,:),r(3,:), '--') % matlab reference
plot(xs,zs); %true path
title('Trajectory vs Matlab Reference')
xlabel('$x$','Interpreter','latex');
ylabel('$z$','Interpreter','latex');
xlim([-2,2]);
ylim([0,4.5]);

% Subplot 2: gazebo reference vs true path
subplot(n,m,2);
hold on;
plot(des_xs, des_zs) % gazebo reference
plot(xs,zs) % true path
title('Trajectory vs Gazebo Reference')
xlim([-2,2]);
ylim([0,4.5]);


function [xs, ys, zs, yaws, des_xs, des_ys, des_zs, des_yaws]=fix_data(M)
    MM = M(:,2);
    ind = ~isnan(MM);
    MMM=MM(ind);
    MMM(end+1) = 0;
    num_rows = 9;
    numCols = length(MMM) / num_rows;
    tempMatrix = reshape(MMM, num_rows, numCols);
    tempMatrix(num_rows, :) = [];
    resultMatrix = tempMatrix;
    
    
    xs=resultMatrix(1,:);
    ys=resultMatrix(2,:);
    zs=resultMatrix(3,:);
    yaws=resultMatrix(4,:);
    des_xs=resultMatrix(5,:);
    des_ys=resultMatrix(6,:);
    des_zs=resultMatrix(7,:);
    des_yaws=resultMatrix(8,:);

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