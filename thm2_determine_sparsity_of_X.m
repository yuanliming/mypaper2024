% the sparse structure of K
SK=[1 1 1 0 ;
    0 1 1 1 ;
    1 0 1 0;];
[m,n]=size(SK);
X=sym('X',[n,n]);
SM=cell(1,m);
M=ones(n,n);
% calculate the sparse structure of X by Theorem 2
for i=1:m
    SM{i}=SK(i,:)'*~SK(i,:);
    M=M.*~SM{i};
end
disp(M)