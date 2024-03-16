function [M,tau,info] = Coeff2Trans(coeffs)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
%% set ellipse parameters
% E = [10 -13 30 15 0.3*pi];
% th=0:pi/180:2*pi;
% Semi_major= E(3);
% Semi_minor= E(4);
% x0= E(1);
% y0= E(2);
% Phi= E(5);
% x=x0+Semi_major*cos(Phi)*cos(th)-Semi_minor*sin(Phi)*sin(th)+1.2*rand(1,361);
% y=y0+Semi_minor*cos(Phi)*sin(th)+Semi_major*sin(Phi)*cos(th)+1.2*rand(1,361);
% scatter(x,y);
%% F(a,x) = Ax^2+Bxy+Cy^2+Dx+Ey+F = 0
% coeffs = [-1/9 1/15 -1/5 -7 -39 -1];   % [-1/5 -1/15 -1/9 7 39 190];  [1 0 1 0 0 1]; 
A=coeffs(1); B=coeffs(2); C=coeffs(3); D=coeffs(4); E=coeffs(5); F=coeffs(6);
if B^2-4*A*C >= 0
    disp('The discriminant B^2-4AC should be less than zero');
    M=0; tau=0; info = 0;
    return
end
%% Matrix computing
Aq = [A B/2 D/2;...
      B/2 C E/2;...
      D/2 E/2 F;];
if sign(trace(Aq(1:2,1:2))) < 0
    Aq = -Aq; % 不会改变原二次曲线形式,去除复数值
end
Q = Aq(1:2,1:2);
[vec, lambda]=eig(Q);
lambda_sqrt=sqrt(lambda);
M = (vec*lambda_sqrt)'; 
tau = M' \ [Aq(1,3);Aq(2,3)];
if tau(1)^2+tau(2)^2-Aq(3,3) < 0
    mu = sqrt(-tau(1)^2-tau(2)^2+Aq(3,3));
else
    mu = sqrt(tau(1)^2+tau(2)^2-Aq(3,3));
end
M = M/mu;
tau = tau/mu;
info = 1;
%% solve 求解
% syms a b c d tx ty
% eq1 = a^2+c^2==A;
% eq2 = a*b+c*d==B/2;
% eq3 = b^2+d^2==C;
% eq4 = a*tx+c*ty==D/2;
% eq5 = b*tx+d*ty==E/2;
% eq6 = tx^2+ty^2-1==F;
% AA = solve([eq1,eq2,eq3,eq4,eq5,eq6],[a,b,c,d,tx,ty]);
%% Parameters of the ellipses
% ellipses_para = coff2param(coeffs/mu);
% % draw the ellipse
% th=0:pi/180:2*pi;
% Semi_major = ellipses_para(3);
% Semi_minor = ellipses_para(4);
% x0 = ellipses_para(1); c = -inv(M)*tau;
% y0 = ellipses_para(2);
% Phi= ellipses_para(5);
% x = x0+Semi_major*cos(Phi)*cos(th)-Semi_minor*sin(Phi)*sin(th);
% y = y0+Semi_minor*cos(Phi)*sin(th)+Semi_major*sin(Phi)*cos(th);

% subplot(121)
% plot(x,y,'color','r', 'LineWidth',2);
% axis equal
% draw the circle

% tform_matrix = [M', [0;0];tau' 1]; % 因为matlab中仿射变换矩阵变了所以这里用H'而不是用HW
% tform = affine2d(tform_matrix);
% [xc,yc] = transformPointsForward(tform,x,y);

% plot(xc,yc,'color','r', 'LineWidth',2);
% xlim([-1,1]);
% axis equal
end

