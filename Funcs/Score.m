function [s] = Score(e,p)
    c=e(1:2)';
    axis=e(3:4)';
    theta=e(5);
%     t = 3;
%     tol_p = c+[cos(theta)*(axis(1)+t); sin(theta)*(axis(1)+t)];
    % 计算椭圆e长轴+t像素上点的D4作为容忍度阈值。
%     tolerance = MahalaDist(tol_p,c,axis,theta,1);
    rho = 0.005;
    Maha_Dist = MahalaDist(p',c,axis,theta,1);
    s_index = find(Maha_Dist<rho);
    s=size(s_index,1);
end


