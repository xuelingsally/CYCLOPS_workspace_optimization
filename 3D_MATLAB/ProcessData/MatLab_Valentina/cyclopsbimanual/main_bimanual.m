% Cyclops forward and inverse kinematics

% Clear all the variables
clear all
clc
% Close all the figures
close all

%%
% THESE ARE THE PARAMETERS YOU CAN SET
% Length of the tool [mm]
L_tool=86;

% Diameter of tool [mm]
d_tool=4;

% Distance between front end of tool and first set of tendon attachment
% points [mm]
x_f=1;

% Distance between first and second set of tendon attachment points on tool
% [mm]
L=40;

% Diameter of the baloon at the tendon feeding point [mm]
d_bal=55;


% Distance between the two sets of tendon feeding points on the baloon [mm]
d_feed=78;

%%
% THESE ARE JUST FOR ME TO DO THE PLOTTING
% Number of slices along the tool's length
no_slice=10;

% Number of verteces on the tool's base
no_arc=50;

%%
% THIS IS THE FUNCTION CALCULATING THE PARAMETERS NEEDED FOR THE KINEMATICS
% Geometrical configuration of the tool at the homing position
[r_tool,s,h1,h2,p_in,T_in,r_ee,p_ee,r1,r2,r3,r4,r5,r5_left,r6,r6_left,B,B_left,P,l,p_x_b]=geometry(L_tool,d_tool,x_f,L,d_bal,d_feed);
% Homing position for left tool
[p_in_left,T_in_left,p_ee_left,P_left,l_left]=homing_left(r_tool,r_ee,r1,r2,r3,r4,r5_left,r6_left,B_left);
% Homing position for right tool
[p_in_right,T_in_right,p_ee_right,P_right,l_right]=homing_right(r_tool,r_ee,r1,r2,r3,r4,r5,r6,B);

%%
% INVERSE KINEMATICS
% These are the values you can set, e.g. the desired motion of the tool tip

% Desired rotation of the tool around Z (rad)
% Need to check which rotation angle of the stylus this corresponds to
alpha=0.2;
alpha_left=-0.2;

% Desired rotation of the tool around Y (rad)
% Need to check which rotation angle of the stylus this corresponds to
beta=0.2;
beta_left=-0.2;

% NB Right now the tool cannot rotate along its own axis but let's keep this
% in for future development
% Desired rotation of the tool around X (rad)
gamma=0;
gamma_left=0;

% Desired new position of the tool end-effector [mm]
% This is the (x,y,z) position of the stylus
p_ee_new=[20,-5,10];
p_ee_new_left=[20,5,10];

% Interpolation step
dr=0.01;

% This function gives you the tendon lengths generating the desired tool tip motion
% NB l_diff is a 6XN matrix with 6 tendon length values at N steps of
% motion; the number of steps depends on the interpolation step dr you set
% above, so if dr=0.01 then N=101 (including step 0 when you start the
% motion)
[l_diff]=cyclops_IK(alpha,beta,gamma,p_ee,p_ee_new,dr,p_in,T_in,r_ee,r1,r2,r3,r4,r5,r6,B,l);
[l_diff_left]=cyclops_IK(alpha_left,beta_left,gamma_left,p_ee_left,p_ee_new_left,dr,p_in_left,T_in_left,r_ee,r1,r2,r3,r4,r5_left,r6_left,B_left,l_left);
[l_diff_right]=cyclops_IK(alpha,beta,gamma,p_ee_right,p_ee_new,dr,p_in_right,T_in_right,r_ee,r1,r2,r3,r4,r5,r6,B,l_right);

%%
% FORWARD KINEMATICS
% In this case you can send to the function the current tendon lengths and
% it will give you the corresponding tool tip pose (position and
% orientation matrix) T_ee
% NB T_interp and P_new are just for me to do the plotting!
[T_ee,T_interp,P_new]=cyclops_FK(l_diff,l,P,B,d_tool,L,s,h1,h2,r_tool,p_x_b,r_ee);
[T_ee_left,T_interp_left,P_new_left]=cyclops_FK(l_diff_left,l_left,P_left,B_left,d_tool,L,s,h1,h2,r_tool,p_x_b,r_ee);
[T_ee_right,T_interp_right,P_new_right]=cyclops_FK(l_diff_right,l_right,P_right,B,d_tool,L,s,h1,h2,r_tool,p_x_b,r_ee);

%%
% FORCE CONTROL
% Read data from file
[p_test,R_test,det_N_test,tau_test]=read_data('test_2.txt',1089);
% Calculation of internal force vector 
for i=1:size(p_test,2)
    [tau(:,i),W(:,:,i),det_N(i)]=internal_force(B,R_test(:,:,i),p_test(:,i),r1,r2,r3,r4,r5,r6);
    % Calculation of tendon tension vector
    [t(:,i),A(:,:,i),rank_A(i)]=tension_vector(B,R_test(:,:,i),p_test(:,i)',r1,r2,r3,r4,r5,r6,r_ee);
end

%%
% WORKSPACE
% Set the minimum allowable tendon tension [N]
t_min=0.01*ones(6,1);

% Set the maximum allowable tendon tension [N]
t_max=120*ones(6,1);

% Set the force applied at the tool end-effector [N]
f_ee=[0,0,0];

% Set the mass of the tool [kg]
m=0.000;

% Set the acceleration of gravity value [m/(s^2)]
g=9.81;

% Set the gravity force vector [N]
G=[0,0,-m*g,0,0,0];

% Set of points to check
n=11;
% Right tool
[p_range_x,p_range_y,p_range_z]=ndgrid(linspace(0,40,n),linspace(-25,-r_tool,n),linspace(-25,25,n));
p_x=reshape(p_range_x,n^3,1);
p_y=reshape(p_range_y,n^3,1);
p_z=reshape(p_range_z,n^3,1);
p=[p_x,p_y,p_z]';

% Left tool
[p_range_x_left,p_range_y_left,p_range_z_left]=ndgrid(linspace(0,40,n),linspace(r_tool,25,n),linspace(-25,25,n));
p_x_left=reshape(p_range_x_left,n^3,1);
p_y_left=reshape(p_range_y_left,n^3,1);
p_z_left=reshape(p_range_z_left,n^3,1);
p_left=[p_x_left,p_y_left,p_z_left]';

% Orientation of the tool at each point
% alpha_range=[-pi/4,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,pi/4];
% beta_range=[-pi/4,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,pi/4];
alpha_range=linspace(-30*pi/180,30*pi/180,13);
beta_range=linspace(-30*pi/180,30*pi/180,13);
count=1;
for i=1:length(alpha_range)
    for j=1:length(beta_range)
        R_p(:,:,count)=rpy2rot(alpha_range(i),beta_range(j),gamma);
        alpha_beta_comb(:,count)=[alpha_range(i);beta_range(j)];
        count=count+1;
    end
end

% Calculate the controllable workspace of the system
% % Find highest and lowest tension solutions
% for j=1:count-1
%     for i=1:size(p,2)
%         [A(:,:,i),t_high(:,i),t_low(:,i),exitflag_high(:,i),exitflag_low(:,i)]=cyclops_WS(t_min,t_max,f_ee,G,B,R_p(:,:,j),p(:,i),r1,r2,r3,r4,r5,r6,r_ee);
%     end
%     % Points with feasible tension solution
%     I_high=find(exitflag_high>0);
%     I_low=find(exitflag_low>0);
%     % Controllable workspace of the centre of mass of the tool
%     WS=p(:,I_high);
%     % Controllable workspace of the tip of the tool
%     tip_WS=WS+repmat(R_p(:,:,j)*r_ee,1,length(I_high));
%     total_tip_WS(j).WS=tip_WS;
% end
% save 'total_tip_WS.mat' total_tip_WS
load 'total_tip_WS.mat'

% % Right tool
% % Find highest and lowest tension solutions
% for j=1:count-1
%     for i=1:size(p,2)
%         [A(:,:,i),t_high(:,i),t_low(:,i),exitflag_high(:,i),exitflag_low(:,i)]=cyclops_WS(t_min,t_max,f_ee,G,B,R_p(:,:,j),p(:,i),r1,r2,r3,r4,r5,r6,r_ee);
%     end
%     % Points with feasible tension solution
%     I_high=find(exitflag_high>0);
%     I_low=find(exitflag_low>0);
%     % Controllable workspace of the centre of mass of the tool
%     WS=p(:,I_high);
%     % Controllable workspace of the tip of the tool
%     tip_WS=WS+repmat(R_p(:,:,j)*r_ee,1,length(I_high));
%     total_tip_WS_right(j).WS=tip_WS;
% end
% 
% save 'overall_tip_ws_no_force_v1_right.mat' total_tip_WS_right;
load 'overall_tip_ws_no_force_v1_right.mat'

% % Left tool
% % Find highest and lowest tension solutions
% for j=1:count-1
%     for i=1:size(p_left,2)
%         [A_left(:,:,i),t_high_left(:,i),t_low_left(:,i),exitflag_high_left(:,i),exitflag_low_left(:,i)]=cyclops_WS(t_min,t_max,f_ee,G,B_left,R_p(:,:,j),p_left(:,i),r1,r2,r3,r4,r5_left,r6_left,r_ee);
%     end
%     % Points with feasible tension solution
%     I_high_left=find(exitflag_high_left>0);
%     I_low_left=find(exitflag_low_left>0);
%     % Controllable workspace of the centre of mass of the tool
%     WS_left=p_left(:,I_high_left);
%     % Controllable workspace of the tip of the tool
%     tip_WS_left=WS_left+repmat(R_p(:,:,j)*r_ee,1,length(I_high_left));
%     total_tip_WS_left(j).WS=tip_WS_left;
% end
% 
% save 'overall_tip_ws_no_force_v1_left.mat' total_tip_WS_left;
load 'overall_tip_ws_no_force_v1_left.mat'



