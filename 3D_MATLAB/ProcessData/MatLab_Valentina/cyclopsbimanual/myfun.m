% Forward kinematics of the cyclops

function F=myfun(x,b1,b2,b3,b4,b5,b6,l1,l2,l3,l4,l5,l6,d,L,s,h1,h2)

p1=x(1:3,1);
p2=x(1:3,2);
p3=x(1:3,3);
p4=x(1:3,4);
p5=x(1:3,5);
p6=x(1:3,6);

F=[norm(p1-b1)-l1;
    norm(p2-b2)-l2;
    norm(p3-b3)-l3;
    norm(p4-b4)-l4;
    norm(p5-b5)-l5;
    norm(p6-b6)-l6;
    norm(p1-p2)-d;
    norm(p3-p4)-d;
    norm(p1-p3)-L;
    norm(p2-p4)-L;
    norm(p5-p6)-L;
    norm(p1-p5)-s;
    norm(p2-p5)-s;
    norm(p3-p6)-s;
    norm(p4-p6)-s;
    norm(p1-p4)-h1;
    norm(p2-p3)-h1;
    norm(p1-p6)-h2;
    norm(p4-p5)-h2;
    norm(p2-p6)-h2;
    norm(p3-p5)-h2];

