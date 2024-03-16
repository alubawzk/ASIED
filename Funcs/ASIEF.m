function [elp,info] = ASIEF(p_x, p_y)
%% normalize input
if isempty(p_x)
    info = 0;
    elp = [];
    return;
end
if size(p_x,1)>size(p_x,2)
    x=p_x(:);
else
    x=p_x(:)';
end
if size(p_y,1)>size(p_y,2)
    y=p_y(:);
else
    y=p_y(:)';
end

%% set intermediate variables
X  = [x.*x x.*y y.*y x y];
z_ = inv(X'*X)*X'*ones(size(x,1),1);
Q  = [z_(1), z_(2)/2;
      z_(2)/2, z_(3)];
q  = [z_(4)/2; z_(5)/2];

%% matrix decomposition
[gevec, geval] = eig(Q);
if( isempty(gevec) )
    info = 0;
    elp = [];
    return;
end
lambda = diag(geval);
b = 1 / (1+q'*inv(Q)*q);
M = [sqrt(b*lambda(1))*gevec(:,1), sqrt(b*lambda(2))*gevec(:,2)]';
t = b*inv(M')*q;

%% get ellipse parameters
A = M'*M;
c = -inv(M)*t;
A1 = A(1,1); A2 = A(1,2); A3 = A(2,1); A4 = A(2,2);
c1 = c(1,1); c2 = c(2,1);
p1 = A1;
p2 = (A3+A2);
p3 = A4;
p4 = -(A1*c1+A2*c2+A1*c1+A3*c2);
p5 = -(A3*c1+A4*c2+A2*c1+A4*c2);
p6 = A1*c1*c1+A3*c2*c1+A2*c1*c2+A4*c2*c2-1;
elp = coff2param([p1,p2,p3,p4,p5,p6]);
info = 1;
end