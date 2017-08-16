function draw_cyclops_eaB(eaB, taskspace, radius_tool, radius_scaffold)

% Euler Angles
ea(1:6) = eaB(1:6);

% Attachment Points
a(:,1) = [eaB(7), radius_tool * sin(ea(1)), radius_tool * cos(ea(1))];
a(:,2) = [eaB(8), radius_tool * sin(ea(2)), radius_tool * cos(ea(2))];
a(:,3) = [eaB(9), radius_tool * sin(ea(3)), radius_tool * cos(ea(3))];
a(:,4) = [eaB(10), radius_tool * sin(ea(4)), radius_tool * cos(ea(4))];
a(:,5) = [eaB(11), radius_tool * sin(ea(5)), radius_tool * cos(ea(5))];
a(:,6) = [eaB(12), radius_tool * sin(ea(6)), radius_tool * cos(ea(6))];

% Base Frame 
B(:,1) = [eaB(13), radius_scaffold * sin(ea(1)), radius_scaffold * cos(ea(1))];
B(:,2) = [eaB(13), radius_scaffold * sin(ea(2)), radius_scaffold * cos(ea(2))];
B(:,3) = [eaB(13), radius_scaffold * sin(ea(3)), radius_scaffold * cos(ea(3))];
B(:,4) = [eaB(14), radius_scaffold * sin(ea(4)), radius_scaffold * cos(ea(4))];
B(:,5) = [eaB(14), radius_scaffold * sin(ea(5)), radius_scaffold * cos(ea(5))];
B(:,6) = [eaB(14), radius_scaffold * sin(ea(6)), radius_scaffold * cos(ea(6))];

dist_tooltip = eaB(15);

x_middle = (min(B(1,:)) + max(B(1,:)))/2;
p = [x_middle; 0; 0];

length_scaffold = abs(eaB(13) - eaB(14));
length_overtube = max(a(1,:)) - min(a(1,:));

dist_tool_b_cg = min(a(1,:));

draw_cyclops_full(a,B, taskspace, radius_tool, radius_scaffold, length_scaffold, length_overtube, dist_tool_b_cg, dist_tooltip);

end
