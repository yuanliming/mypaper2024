tic
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
f= 4.8; % friction
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
%  B1=eye(6);
%  B1=eye(6);
C=[100000 0 0 0 0 0 0 0 0;
   0 0 0 1000 0 0 0 0 0;
   0 0 0 0 0 0 0 0 0;
   0 0 0 0 0 0 0 0 0];
% C=[1 0 0 0 0 0 0 0 0;
%    0 0 0 1 0 0 0 0 0;
%    0 0 0 0 0 0 0 0 0;
%    0 0 0 0 0 0 0 0 0];
D=[0 0 ;
   0 0 ;
   1 0 ;
   0 1;];
%%%%system in the parameter space
F1 = [A1          -B21;
      zeros(2,9) zeros(2,2)];
F2 = [A2          -B22;
      zeros(2,9) zeros(2,2)];
Q = [B1*B1'     zeros(9,2);
      zeros(2,9) zeros(2,2)]; 
R = blkdiag(C'*C,D'*D); 

%%%%optimization in the parameter space
W1 = blkdiag(sdpvar(3,3),sdpvar(3,3),sdpvar(3,3));
W2 = [sdpvar(1,3),zeros(1,3),sdpvar(1,3);zeros(1,3),sdpvar(1,3),sdpvar(1,3)]';
W3 = sdpvar(2,2);
W = [W1  W2;
     W2' W3];
V=[eye(9)  zeros(9,2)];
z=trace(R*W);
MyLMI = [W >= 0, V*(F1*W + W*F1'+Q)*V' <= 0, V*(F2*W + W*F2'+Q)*V' <= 0];
options=sdpsettings('solver','sdpt3');%sdpt3
optimize(MyLMI, z,options);
K_qs=value(W2)'*value(W1)^-1; % controller
H2norm_upper=sqrt(value(z)); % optimized H2 norm upper bound
toc
% eig(A-B2*K)
H2norm_upper
vertex1=norm(ss(A1-B21*K_qs,B1,C-D*K_qs,zeros(4,1)),2)
vertex1=norm(ss(A1-B22*K_qs,B1,C-D*K_qs,zeros(4,1)),2)