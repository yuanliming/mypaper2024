tic
alpha=0.01;
m1= 0.128; % mass of the fine stage
m11=0.128*1.2;
m12=0.128*0.8;
m2= 1.977; % mass of the coarse stage
b1= 17.7; % counter electromotive force of VCM
b2= 19.6; % counter electromotive force of linear motor
k1= 6797; % stiffness of the flexure 用测力计测一下
kf1= 17.7;% VCM force constant N/A 
kv1= 0.2; % VCM driver gain A/V
kf2= 24;  % linear PM motor force constant N/A
kv2= 1; % linear PM motor driver gain A/V
f= 4.8;  % friction
A1=[0 1 0 0 0 0 0 0 0;
   0 0 1 0 0 0 0 0 0;
   0 -k1/m11 -b1/m11 0 k1/m11 b1/m11 0 0 0;
   0 0 0 0 1 0 0 0 0;
   0 0 0 0 0 1 0 0 0;
   0 k1/m2 b1/m2 0 -k1/m2 -b2/m2-b1/m2 0 0 b2/m2;
   0 0 0 0 0 0 0 1 0;
   0 0 0 0 0 0 0 0 1;
   0 0 0 0 0 0 -64 -48 -12];
A2=[0 1 0 0 0 0 0 0 0;
   0 0 1 0 0 0 0 0 0;
   0 -k1/m12 -b1/m12 0 k1/m12 b1/m12 0 0 0;
   0 0 0 0 1 0 0 0 0;
   0 0 0 0 0 1 0 0 0;
   0 k1/m2 b1/m2 0 -k1/m2 -b2/m2-b1/m2 0 0 b2/m2;
   0 0 0 0 0 0 0 1 0;
   0 0 0 0 0 0 0 0 1;
   0 0 0 0 0 0 -64 -48 -12];
B21=[0 0;
    0 0;
    -kf1*kv1/m11 0;
    0 0;
    0 0;
    kf1*kv1/m2 -kf2*kv2/m2;
    0 0;
    0 0;
    0 0];
B22=[0 0;
    0 0;
    -kf1*kv1/m12 0;
    0 0;
    0 0;
    kf1*kv1/m2 -kf2*kv2/m2;
    0 0;
    0 0;
    0 0];
B1=[0;0;0;0;0;f/m2;0;0;0];%

C=[100000 0 0 0 0 0 0 0 0;
   0 0 0 1000 0 0 0 0 0;
   0 0 0 0 0 0 0 0 0;
   0 0 0 0 0 0 0 0 0];

D=[0 0 ;
   0 0 ;
   1 0 ;
   0 1;];

P1 = sdpvar(9,9);
Z1 = sdpvar(1,1);
P2 = sdpvar(9,9);
Z2 = sdpvar(1,1);

X = [sdpvar(3,3,'full')     zeros(3,3)     sdpvar(3,3,'full');
    zeros(3,3)     sdpvar(3,3,'full')     sdpvar(3,3,'full');
    zeros(3,3)     zeros(3,3)     sdpvar(3,3,'full');];

R = [sdpvar(1,3),zeros(1,3),sdpvar(1,3);zeros(1,3),sdpvar(1,3),sdpvar(1,3)];

mat11 = -[zeros(9) P1 zeros(9,4);
         P1 zeros(9) zeros(9,4); 
         zeros(4,9) zeros(4,9) -eye(4)]-[A1*X-B21*R;-X;C*X-D*R]*[eye(9) alpha*eye(9) zeros(9,4)]-([A1*X-B21*R;-X;C*X-D*R]*[eye(9) alpha*eye(9) zeros(9,4)])';
mat12 = -[zeros(9) P2 zeros(9,4);
         P2 zeros(9) zeros(9,4); 
         zeros(4,9) zeros(4,9) -eye(4)]-[A2*X-B22*R;-X;C*X-D*R]*[eye(9) alpha*eye(9) zeros(9,4)]-([A2*X-B22*R;-X;C*X-D*R]*[eye(9) alpha*eye(9) zeros(9,4)])';

mat21 = [Z1 B1'; B1 P1];
mat22 = [Z2 B1'; B1 P2];

trz1 = trace(Z1);  
trz2 = trace(Z2);

gamma = sdpvar(1);
LMI = [mat21,mat22, mat11,mat12, P1, P2, Z1>=0, Z2>=0, trz1<=gamma, trz2<=gamma];
options=sdpsettings('solver','sdpt3'); % sdpt3
optimize(LMI, gamma,options); 

% result
K_sv = value(R)*inv(value(X)); % controller
H2norm_upper = sqrt(value(gamma)); % optimized H2 norm upper bound
toc
H2norm_upper
 vertex1=norm(ss(A1-B21*K_sv,B1,C-D*K_sv,zeros(4,1)),2)
 vertex1=norm(ss(A1-B22*K_sv,B1,C-D*K_sv,zeros(4,1)),2)