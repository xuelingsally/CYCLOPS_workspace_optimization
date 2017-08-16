% Forward kinematics of the cyclops

function [T_ee,T_interp,P_new]=cyclops_FK(l_diff,l,P,B,d_tool,L,s,h1,h2,r_tool,p_x_b,r_ee)

r_pt=size(l_diff,2);

l_new=l_diff+repmat(l',1,r_pt);

x0=P;
options.MaxFunEvals=30000;
options.MaxIter=2000;
options.ScaleProblem='Jacobian';
for i=1:r_pt
    x=fsolve(@(x) myfun(x,B(:,1),B(:,2),B(:,3),B(:,4),B(:,5),B(:,6),l_new(1,i),l_new(2,i),l_new(3,i),l_new(4,i),...
        l_new(5,i),l_new(6,i),d_tool,L,s,h1,h2),x0,options);
    P_new(:,:,i)=x;
    x0=x;
    p_mid_f=x(:,2)+r_tool*(x(:,1)-x(:,2))/norm(x(:,1)-x(:,2));
    p_mid_b=x(:,4)+r_tool*(x(:,3)-x(:,4))/norm(x(:,3)-x(:,4));
    p_new(:,i)=p_mid_b+p_x_b*(p_mid_f-p_mid_b)/norm(p_mid_f-p_mid_b);
    q=(p_new(:,i)-p_mid_f);
    q_new=q/norm(p_new(:,i)-p_mid_f);
    theta_y=-asin(q_new(3));
    theta_z=asin(q_new(2)/cos(theta_y));
    R_new(:,:,i)=rpy2rot(theta_z,theta_y,0);
    T_interp(:,:,i) = [R_new(:,:,i) p_new(:,i); 0 0 0 1];
end

for i=1:r_pt
    % Position of the end-effector with respect to the inertial base frame[mm]
    p_ee(:,i)=p_new(:,i)+R_new(:,:,i)*r_ee;
    % Pose of the end-effector with respect to the inertial base frame
    T_ee(:,:,i)=[R_new(:,:,i) p_ee(:,i); 0 0 0 1];
end