load('L.mat');
load('R.mat');

c_x = ( max(dataR(:,1)) - min(dataR(:,1)) )/2;
c_y = ( max(dataR(:,2)) - min(dataR(:,2)) )/2;
c_z = ( max(dataR(:,3)) - min(dataR(:,3)) )/2;

dataR_r(:,1) = dataR(:,1) + c_x - max(dataR(:,1));
dataR_r(:,2) = dataR(:,2) + c_x - max(dataR(:,2));
dataR_r(:,3) = dataR(:,3) + c_x - max(dataR(:,3));
dataR_r(:,4) = 0;
dataR_r(:,5) = dataR(:,5);
dataR_r(:,6) = dataR(:,4);

c_x = ( max(dataL(:,1)) - min(dataL(:,1)) )/2;
c_y = ( max(dataL(:,2)) - min(dataL(:,2)) )/2;
c_z = ( max(dataL(:,3)) - min(dataL(:,3)) )/2;

dataL_r(:,1) = dataL(:,1) + c_x - max(dataL(:,1));
dataL_r(:,2) = dataL(:,2) + c_x - max(dataL(:,2));
dataL_r(:,3) = dataL(:,3) + c_x - max(dataL(:,3));
dataL_r(:,4) = 0;
dataL_r(:,5) = dataL(:,5);
dataL_r(:,6) = dataL(:,4);

dataL_r(:,7) = dataL(:,9); % x is z
dataL_r(:,8) = dataL(:,7); % y is x
dataL_r(:,9) = dataL(:,8); % z is y
