
c_x = ( max(dataR(:,1)) - min(dataR(:,1)) )/2;
c_y = ( max(dataR(:,2)) - min(dataR(:,2)) )/2;
c_z = ( max(dataR(:,3)) - min(dataR(:,3)) )/2;

dataR_r(:,1) = dataR(:,1) + c_x - max(dataR(:,1));
dataR_r(:,2) = dataR(:,2) + c_x - max(dataR(:,2));
dataR_r(:,3) = dataR(:,3) + c_x - max(dataR(:,3));
dataR_r(:,4) = 0;
dataR_r(:,5) = dataR(:,5);
dataR_r(:,6) = dataR(:,4);