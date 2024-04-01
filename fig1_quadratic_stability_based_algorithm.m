un=0.1; % uncertainty level
N=16; % number of extreme systems
up=1+un; % upper bound
lo=1-un; % lower bound

% the plant 
A=[0.19,0.25,0.18,0.70;
    0.31,0.89,0.21,0.55;
    0.30,0.70,0.07,0.31;
    0.21,0.55,0.91,0.16];
B2=[0.62,0.39,0.98;
    0.96,0.07,0.40;
    0.17,0.68,0.62;
    0.25,0.40,0.15];
% impose model uncertainties
A1=A;
A1(1,1)=A(1,1)*up;A1(3,4)=A(3,4)*up;A1(4,2)=A(4,2)*up;
A2=A;
A1(1,1)=A(1,1)*up;A1(3,4)=A(3,4)*up;A1(4,2)=A(4,2)*lo;
A3=A;
A1(1,1)=A(1,1)*up;A1(3,4)=A(3,4)*lo;A1(4,2)=A(4,2)*up;
A4=A;
A1(1,1)=A(1,1)*up;A1(3,4)=A(3,4)*lo;A1(4,2)=A(4,2)*lo;
A5=A;
A1(1,1)=A(1,1)*lo;A1(3,4)=A(3,4)*up;A1(4,2)=A(4,2)*up;
A6=A;
A1(1,1)=A(1,1)*lo;A1(3,4)=A(3,4)*up;A1(4,2)=A(4,2)*lo;
A7=A;
A1(1,1)=A(1,1)*lo;A1(3,4)=A(3,4)*lo;A1(4,2)=A(4,2)*up;
A8=A;
A1(1,1)=A(1,1)*lo;A1(3,4)=A(3,4)*lo;A1(4,2)=A(4,2)*lo;
B21=B2;
B21(2,1)=B21(2,1)*up;
B22=B2;
B22(2,1)=B22(2,1)*lo;
% construct the 16 extreme systems
A=cell(1,16);
A{1}=A1;A{2}=A2;A{3}=A3;A{4}=A4;A{5}=A5;A{6}=A6;A{7}=A7;A{8}=A8;
A{9}=A1;A{10}=A2;A{11}=A3;A{12}=A4;A{13}=A5;A{14}=A6;A{15}=A7;A{16}=A8;
B2=cell(1,16);
B2{1}=B21;B2{2}=B21;B2{3}=B21;B2{4}=B21;B2{5}=B21;B2{6}=B21;B2{7}=B21;B2{8}=B21;
B2{9}=B22;B2{10}=B22;B2{11}=B22;B2{12}=B22;B2{13}=B22;B2{14}=B22;B2{15}=B22;B2{16}=B22;

B1=eye(4);
C=[1 0 0 0;
    0 1 0 0;
    0 0 0 0;
    0 0 0 0;
    0 0 0 0];
D=[0 0 0;0 0 0;1 0 0;0 1 0;0 0 1];

V=[eye(4)  zeros(4,3)]; 
Q = [B1*B1'     zeros(4,3);
      zeros(3,4) zeros(3,3)]; 
% the diagonal Lyapunov matrix that is compatible to the desired structure  
W1 = blkdiag(sdpvar(1,1),sdpvar(1,1),sdpvar(1,1),sdpvar(1,1)); 
W2 = [sdpvar(1,1),sdpvar(1,1),sdpvar(1,1),zeros(1,1);
    zeros(1,1),sdpvar(1,1),sdpvar(1,1),sdpvar(1,1)
    sdpvar(1,1),zeros(1,1),sdpvar(1,1),zeros(1,1)]'; 
W3 = sdpvar(3,3);
W = [W1  W2;
     W2' W3];

 % use a single Lyapunov matrix to deal with all the uncertainties  
for i = 1:N; 
    F{i} = [A{i} -B2{i};
        zeros(3,4) zeros(3,3)];
    T{i}=-V*(F{i}*W + W*F{i}'+Q)*V'; 
end
R = blkdiag(C'*C,D'*D); 

% establish the optimization for constrained optimal robust controller synthesis
z=trace(R*W);
MyLMI = [W >= 0,W1>=0, T{1:N}];
options=sdpsettings('solver','sdpt3');
optimize(MyLMI, z,options)

% result
K_qs=value(W2)'*value(W1)^-1 % controller
H2norm_upper=sqrt(value(z)) % optimized H2 norm upper bound

% calculate the H2 costs at each vertex of the polytopic undertain domain
x=1:1:N;
h2norm_Kqs=[];
for i=1:N
    hnorm=norm(ss(A{i}-B2{i}*K_qs,B1,C-D*K_qs,zeros(5,4)),2);
    h2norm_Kqs=[h2norm_Kqs,hnorm];
end
figure
scatter(x,h2norm_Kqs,'d','r','filled')
set(gca,'FontSize',16,'Fontname', 'Times New Roman')
leg=legend('$K_{\rm QS}$');
xla=xlabel('Extreme systems');
yla=ylabel('$H_2$ cost');
set(leg,'interpreter','latex')
set(xla,'interpreter','latex')
set(yla,'interpreter','latex')
grid on