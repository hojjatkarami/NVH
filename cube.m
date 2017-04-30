function P = cube(T,eul)

A = [-.5 -.5 -.5];
B = [.5 -.5 -.5];
C = [-.5 .5 -.5];
D = [-.5 -.5 .5];
E = [-.5 .5 .5];
F = [.5 -.5 .5];
G = [.5 .5 -.5];
H = [.5 .5 .5];
P = [A;B;F;H;G;C;A;D;E;H;F;D;E;C;G;B];
R = [cos(eul(3))*cos(eul(2)) -sin(eul(3))*cos(eul(1))+cos(eul(3))*sin(eul(2))*sin(eul(1)) sin(eul(3))*sin(eul(1))+cos(eul(3))*sin(eul(2))*cos(eul(1));
       sin(eul(3))*cos(eul(2)) cos(eul(3))*cos(eul(1))+sin(eul(3))*sin(eul(2))*sin(eul(1)) -cos(eul(3))*sin(eul(1))+sin(eul(3))*sin(eul(2))*cos(eul(1));
       -sin(eul(2)) cos(eul(2))*sin(eul(1)) cos(eul(2))*cos(eul(1))]';
P(:,1) = P(:,1) + T(1);
P(:,2) = P(:,2) + T(2);
P(:,3) = P(:,3) + T(3);
P = P*R;
