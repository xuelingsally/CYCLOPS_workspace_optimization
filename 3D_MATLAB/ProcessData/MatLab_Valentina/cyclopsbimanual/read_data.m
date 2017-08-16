% Read data from text file

function [p_test,R_test,det_N,tau_test]=read_data(filename,N)

% Open text file
fid=fopen(filename);

% Read text file
n=1;
for i=1:N
    C=fgetl(fid);
    C_num=str2num(C);
    if ~isempty(C_num)          
        p_test(1,n)=C_num(1);
        p_test(2,n)=C_num(2);
        p_test(3,n)=C_num(3);
        R_test(1,1,n)=C_num(4);
        R_test(1,2,n)=C_num(5);
        R_test(1,3,n)=C_num(6);
        R_test(2,1,n)=C_num(7);
        R_test(2,2,n)=C_num(8);
        R_test(2,3,n)=C_num(9);
        R_test(3,1,n)=C_num(10);
        R_test(3,2,n)=C_num(11);
        R_test(3,3,n)=C_num(12);
        det_N(n)=C_num(13);
        tau_test(1,n)=C_num(14);
        tau_test(2,n)=C_num(15);
        tau_test(3,n)=C_num(16);
        tau_test(4,n)=C_num(17);
        tau_test(5,n)=C_num(18);
        tau_test(6,n)=C_num(19);
        n=n+1;
    end
end

% Close text file
fclose(fid);
