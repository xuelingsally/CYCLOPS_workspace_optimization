
clear all
close all
clc

addpath MatLab_Valentina/cyclopsbimanual
%%

M =  csvread('ESD_singleInstrument.csv',7,2);
M = M(:,1:7);      % Cut of last column, is empty

Translation_taskspace = 1000*M(:,5:7)';
Quaternion_taskspace = M(:,1:4)';

Trans_frame = [-1 0 0; 0 0 1; 0 1 0];
Translation_taskspace = Trans_frame*Translation_taskspace;

figure
plot3(Translation_taskspace(1,:),Translation_taskspace(2,:),Translation_taskspace(3,:),'.b')
axis equal

for i = 1:length(Translation_taskspace(1,:))
    if Translation_taskspace(1,i) < -91
        Translation_taskspace(:,i) = zeros(3,1);
        Quaternion_taskspace(:,i) = zeros(4,1);
    end
    if Translation_taskspace(3,i) > 22.5
        Translation_taskspace(:,i) = zeros(3,1);
        Quaternion_taskspace(:,i) = zeros(4,1);
    end
end

Translation_taskspace(:,any(Translation_taskspace,1)==0)=[];
Quaternion_taskspace(:,any(Quaternion_taskspace,1)==0)=[];

hold on    
plot3(Translation_taskspace(1,:),Translation_taskspace(2,:),Translation_taskspace(3,:),'.r')
axis equal    



for i = 1:length(Quaternion_taskspace(1,:))
%     norm quaternion
    Rotation_taskspace(:,:,i) = quat2rot([Quaternion_taskspace(end,i); Quaternion_taskspace(1:3,i)]);
    temprot2rpy = rot2rpy(Rotation_taskspace(:,:,i));
    %[Rotation_angle(3,i),Rotation_angle(2,i),Rotation_angle(1,i)] = temprot2rpy;
    Rotation_angle(3,i) = temprot2rpy(1,1);
    Rotation_angle(2,i) = temprot2rpy(1,2);
    Rotation_angle(1,i) = temprot2rpy(1,3);
    
end
    
Rotation_angle = Trans_frame*Rotation_angle;
Rotation_angle(1,:) = zeros(size(Rotation_angle(1,:)));
p_end_effector = [45; 0; -6.5];


for i =1:length(Quaternion_taskspace(1,:))
    R(:,:,i) = rpy2rot(Rotation_angle(3,i),Rotation_angle(2,i),Rotation_angle(1,i));
    p(:,i) = R(:,:,i)*p_end_effector + Translation_taskspace(:,i);
    
end

plot3(p(1,:),p(2,:),p(3,:),'.r')
axis equal 

% Relative Taskspace
mean_p = mean(p,2);
p_relative = [];
for i=1:size(p,2)
    p_relative(1:3,i) = p(:,i) - mean_p;
    p_relative(4:6,i) = Rotation_angle(:,i);
end

data_r_tp = p_relative;
save('data_r_tp.mat', 'data_r_tp');

% Relative trans and rot
% max_x = max(Translation_taskspace(1,:));
% min_x = min(Translation_taskspace(1,:));
% max_y = max(Translation_taskspace(2,:));
% min_y = min(Translation_taskspace(2,:));
% max_z = max(Translation_taskspace(3,:));
% min_z = min(Translation_taskspace(3,:));
% 
% cx = min_x + (max_x - min_x)/2;
% cy = min_y + (max_y - min_y)/2;
% cz = min_z + (max_z - min_z)/2;

% Relative_Translation = Translation_taskspace - repmat([cx;cy;cz],1,length(Translation_taskspace(1,:)));
% figure 
% plot3(Relative_Translation(1,:),Relative_Translation(2,:),Relative_Translation(3,:),'.b')
% axis equal

% minry = min(Rotation_angle(2,:));
% minrz = min(Rotation_angle(3,:));
% 
% maxry = max(Rotation_angle(2,:));
% maxrz = max(Rotation_angle(3,:));
% 
% Relative_Rotation = Rotation_angle - repmat([0;minry;minrz],1,length(Rotation_angle(1,:)));
% 
% 
% Matrix_trans_rot = [Translation_taskspace; Rotation_angle];
% Matrix_relative_trans_rot = [Relative_Translation; Relative_Rotation];
% 
% Matrix_tracking_data(1:5,:) = Matrix_relative_trans_rot(1:5,:);
% Matrix_tracking_data(6,:) = Matrix_trans_rot(6,:);
% 

