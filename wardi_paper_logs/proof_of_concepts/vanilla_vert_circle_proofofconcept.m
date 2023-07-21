clc;clear;close all
% PLOT THE HOVER SEE IF IT EVENTUALLY GETS TO THE TRUE VALUE



% Gather References to Compare Against
[vertical_circle, horizontal_circle, fig8_horz, fig8_vert_short, fig8_vert_tall] = get_reference_trajectories();

% Vertical Circle
n=2;m=2;
figure()
sgtitle('PID Control Vertical Circle')

subplot(n,m,1);
hold on;
r = vertical_circle;
plot(r(1,:),r(3,:), '--')
M = readmatrix('new.csv');
[xs, ys, zs, yaws, des_xs, des_ys, des_zs, des_yaws]=fix_data(M);
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
