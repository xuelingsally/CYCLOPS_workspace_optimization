close all;

addpath data;
load('./data/L.mat');
load('./data/R.mat');
relative_tp;

radius_tool = 1.75;
radius_scaffold = 26.5491;
%radius_scaffold = 28.64;
%radius_scaffold = 30.51;

%load('data_r_tp.mat');
%taskspace = data_r_tp;
taskspace = dataL_r';
taskspace(1,:) = taskspace(1,:) + abs(min(taskspace(1,:))) + 5;
%taskspace(2,:) = taskspace(2,:) + radius_scaffold/2 - 7.4;
taskspace(2,:) = taskspace(2,:) - radius_scaffold/2;
taskspace(3,:) = taskspace(3,:) + abs(min(taskspace(3,:))) - radius_scaffold - 3;
%taskspace(3,:) = taskspace(3,:) - 13;
data1 = taskspace;


for iter=34:34
    
    eaB = eaBv(:,iter);

    ea(1:6) = eaB(1:6);

    % Attachment Points
    a(:,1) = [eaB(7), radius_tool * cos(ea(1)), radius_tool * sin(ea(1))];
    a(:,2) = [eaB(8), radius_tool * cos(ea(2)), radius_tool * sin(ea(2))];
    a(:,3) = [eaB(9), radius_tool * cos(ea(3)), radius_tool * sin(ea(3))];
    a(:,4) = [eaB(10), radius_tool * cos(ea(4)), radius_tool * sin(ea(4))];
    a(:,5) = [eaB(11), radius_tool * cos(ea(5)), radius_tool * sin(ea(5))];
    a(:,6) = [eaB(12), radius_tool * cos(ea(6)), radius_tool * sin(ea(6))];

    % Base Frame 
    B(:,1) = [eaB(13), radius_scaffold * cos(ea(1)), radius_scaffold * sin(ea(1))];
    B(:,2) = [eaB(13), radius_scaffold * cos(ea(2)), radius_scaffold * sin(ea(2))];
    B(:,3) = [eaB(13), radius_scaffold * cos(ea(3)), radius_scaffold * sin(ea(3))];
    B(:,4) = [eaB(14), radius_scaffold * cos(ea(4)), radius_scaffold * sin(ea(4))];
    B(:,5) = [eaB(14), radius_scaffold * cos(ea(5)), radius_scaffold * sin(ea(5))];
    B(:,6) = [eaB(14), radius_scaffold * cos(ea(6)), radius_scaffold * sin(ea(6))];


    add_cables = (size(eaB,1) - 18) / 3;

    for i=1:add_cables
        ea(end+1) = eaB(19+(i-1)*3);
        a(:,end+1) = [eaB(20+(i-1)*3), radius_tool * cos(ea(end)), radius_tool * sin(ea(end))];
        B(:,end+1) = [eaB(21+(i-1)*3), radius_scaffold * cos(ea(end)), radius_scaffold * sin(ea(end))];
    end

    gamma_y = eaB(16);
    gamma_z = eaB(17);

    R_y = [cos(gamma_y), 0, sin(gamma_y); 0, 1, 0; -sin(gamma_y), 0, cos(gamma_y)];
    R_z = [cos(gamma_z), -sin(gamma_z), 0; sin(gamma_z), cos(gamma_z), 0; 0, 0, 1];
    R = R_z * R_y;

    curve_length_x = eaB(15) - eaB(18);

    r_curve = [eaB(18); 0; 0];
    
    
    r_ee_dir = R * [1;0;0];
    r_ee_dir_x = r_ee_dir(1);
    r_ee_dir = r_ee_dir / r_ee_dir_x * curve_length_x;

    r_ee = r_ee_dir + [eaB(18); 0; 0];

    curve_length = norm(r_ee_dir);

    curve_x = eaB(18);

    t_min = ones(6,1) * 1;
    t_max = ones(6,1) * 60;

    W = zeros(1,6);

    %f_ee = [0,0,-0.1];
    feas = [];
    unfeas = [];
    
    for i=1:size(data1,2)
        alpha_y = data1(5,i);
        alpha_z = data1(6,i);

        R_y = [cos(alpha_y), 0, sin(alpha_y); 0, 1, 0; -sin(alpha_y), 0, cos(alpha_y)];
        R_z = [cos(alpha_z), -sin(alpha_z), 0; sin(alpha_z), cos(alpha_z), 0; 0, 0, 1];
        R = R_z * R_y;

        temp = R * [curve_length;0;0];

        data_t1 = -temp + data1(1:3,i);

        beta_y = alpha_y - gamma_y;
        beta_z = alpha_z - gamma_z;

        R_y_2 = [cos(beta_y), 0, sin(beta_y); 0, 1, 0; -sin(beta_y), 0, cos(beta_y)];
        R_z_2 = [cos(beta_z), -sin(beta_z), 0; sin(beta_z), cos(beta_z), 0; 0, 0, 1];
        R_2 = R_z_2 * R_y_2;

        temp2 = R_2 * r_curve;
        data_t2(1:3) = -temp2 + data_t1;
        data_t2(4:6) = [0;beta_y;beta_z];
        
        P = data_t2(:);
        P(1:3) = P(1:3)/1000;
        P(4:5) = P(5:6);
        P(6) = [];
        
        f_ee = taskspace(7:9,i)';
        
        P = P';
        
        [~, feasible] = feasible_pose(P, a/1000, B/1000, W, f_ee, r_ee'/1000, t_min, t_max);

        if feasible == 1
            feas(:,end+1) = data1(:,i);
        else
            unfeas(:,end+1) = data1(:,i);
        end
    end
    draw_cyclops_curved(eaB, feas , radius_tool, radius_scaffold);
    if (size(unfeas) >0)
        plot3(unfeas(1,:), unfeas(2,:), unfeas(3,:), 'r.');
    end;
    view(135,20);
    grid on;
    
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 10 10]);
    
    temp_string = strcat('./images/',num2str(iter),'.png');
    saveas(gcf, temp_string);
    close all;
    
end;

