function [dist] = dist_two__Elps(elp1, elp2, coeff1, coeff2)
K = [0, 90, 180, 270];
K = (K .* pi) ./ 180; 
[M1,t1,info1] = Coeff2Trans(coeff1);
[M2,t2,info2] = Coeff2Trans(coeff2);
if ~info1 || ~info2
    dist = 1000;
    return
end

%% dist1
c2 = elp2(:,1:2);
axis2 = elp2(:,3:4);
theta2 = elp2(:,5);
px=c2(1)+axis2(1)*cos(theta2)*cos(K)-axis2(2)*sin(theta2)*sin(K);
py=c2(2)+axis2(2)*cos(theta2)*sin(K)+axis2(1)*sin(theta2)*cos(K);
p = [px'  py'];
tform_matrix = [M1', [0;0];t1' 1];
tform = affine2d(tform_matrix);
[ux,uy] = transformPointsForward(tform,p(:,1),p(:,2));
d1 = abs(sqrt((ux).^2+(uy).^2)-1);
dist1 = sum(d1) / size(d1,1);

%% dist2
c1 = elp1(:,1:2);
axis1 = elp1(:,3:4);
theta1 = elp1(:,5);
px = c1(1)+axis1(1)*cos(theta1)*cos(K)-axis1(2)*sin(theta1)*sin(K);
py = c1(2)+axis1(2)*cos(theta1)*sin(K)+axis1(1)*sin(theta1)*cos(K);
p = [px' py'];
tform_matrix = [M2', [0;0];t2' 1];
tform = affine2d(tform_matrix);
[ux,uy] = transformPointsForward(tform,p(:,1),p(:,2));
d2 = abs(sqrt((ux).^2+(uy).^2)-1);
dist2 = sum(d2) / size(d2,1);

%% output
dist = exp(-0.5*(dist1+dist2));

end

