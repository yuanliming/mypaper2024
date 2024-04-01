un=0.1;
N=16;
revo_alpha=[];

up=1+un;
lo=1-un;
A=[0.19,0.25,0.18,0.70;
    0.31,0.89,0.21,0.55;
    0.30,0.70,0.07,0.31;
    0.21,0.55,0.91,0.16];
B2=[0.62,0.39,0.98;
    0.96,0.07,0.40;
    0.17,0.68,0.62;
    0.25,0.40,0.15];

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

for alpha=0.05:0.05:5 % optimize the H2norm upper bound with differnt alpha

for i = 1:N; 
    X{i} = sdpvar(4,4);
    Z{i} = sdpvar(4,4);
end
G = [sdpvar(1,1)     zeros(1,1)     sdpvar(1,1)     zeros(1,1);
     zeros(1,1)     sdpvar(1,1)     sdpvar(1,1)     zeros(1,1);
     zeros(1,1)     zeros(1,1)     sdpvar(1,1)     zeros(1,1);
     zeros(1,1)     sdpvar(1,1)    sdpvar(1,1)     sdpvar(1,1)];
W = [sdpvar(1,1) sdpvar(1,1) sdpvar(1,1) zeros(1,1);
  zeros(1,1) sdpvar(1,1) sdpvar(1,1) sdpvar(1,1) ;
    sdpvar(1,1) zeros(1,1) sdpvar(1,1) zeros(1,1)];
gamma = sdpvar(1);
for i = 1:N; 
    mat1{i} = -[zeros(4) X{i} zeros(4,5);
         X{i} zeros(4) zeros(4,5); 
         zeros(5,4) zeros(5,4) -eye(5)]-[A{i}*G-B2{i}*W;-G;C*G-D*W]*[eye(4) alpha*eye(4) zeros(4,5)]-([A{i}*G-B2{i}*W;-G;C*G-D*W]*[eye(4) alpha*eye(4) zeros(4,5)])';
end
for i = 1:N; mat2{i} = [Z{i} B1'; B1 X{i}];end  
for i = 1:N; trz{i} = trace(Z{i});end  
LMI = [mat2{1:N}, mat1{1:N}, X{1:N}, Z{i}>=0, trz{i}<=gamma];
options=sdpsettings('solver','sdpt3'); % sdpt3
optimize(LMI, gamma,options); 
%K_sv = value(W)*inv(value(G));
H2norm_upper = sqrt(value(gamma))
 revo_alpha=[revo_alpha H2norm_upper];
end
alpha=0.05:0.05:5;
plot(alpha,revo_alpha,'b','LineWidth',1.5)
hold on 
plot([0,5],[3.3479,3.3479],'r--','LineWidth',1.5)
box on
grid on
set(gca,'FontSize',16,'Fontname', 'Times New Roman')
xla=xlabel('$\alpha$');
yla=ylabel('$\beta$');%Upper bound to the $H_2$ norm 
set(xla,'interpreter','latex')
set(yla,'interpreter','latex')
legend('Slacking-variable method','Quadratic stability based method')