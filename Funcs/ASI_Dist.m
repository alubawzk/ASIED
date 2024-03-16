function Maha_Dist = MahalaDist(x, c, axis, theta, r)
% 输入：x:[x1,y1;x2,y2;...xn,yn], c:椭圆中心[cx,cy], axis:[a,b], r: 马氏距离中圆的半径
    warning('off');
    U = [cos(theta), -sin(theta); sin(theta), cos(theta)];
    T = [axis(1)^2, 0; 0, axis(2)^2];
    S = U * T * U';
    points = x;
    inv_S = inv(S);
    x = points(1,:)-c(1);
    y = points(2,:)-c(2);
    dis = inv_S(1,1)*x.^2+inv_S(2,2)*y.^2+(inv_S(1,2)+inv_S(2,1)).*x.*y;
    Maha_Dist = abs(dis.^0.5 - r);
end