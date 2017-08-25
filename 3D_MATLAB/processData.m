function output = processData(DataMat)

% DataMat_x = DataMat(1,:);
% 
% min_x = min(DataMat_x);
% 
% translate_x = abs(min_x) + T_Dist(1);
% translate_z = abs(min_z)

fileID = fopen('taskspace.txt','w');
fprintf(fileID, 'TP_TYPE 2\n');

DataMat(4:5,:) = DataMat(5:6,:);
DataMat(6,:) = [];

output = [];

counter2 = 1;
counter = 5;
for i=1:size(DataMat,2)
    if counter == 5
        %DataMat(1,i) = DataMat(1,i) + translate_x;
        output(:,counter2) = DataMat(:,i);
        fprintf(fileID, '%4.3f %4.3f %4.3f %4.3f %4.3f \n', DataMat(:,i));
        counter = 0;
        counter2 = counter2 + 1;
    end
    counter = counter + 1;
end

fclose(fileID);

end