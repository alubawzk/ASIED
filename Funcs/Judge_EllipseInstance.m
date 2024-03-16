function flag = Judge_EllipseInstance(ellipse_i, ellipse_j, Mahala_tolerance)
    flag = false;
    cx = ellipse_j(1);
    cy = ellipse_j(2);
    Semi_major = ellipse_j(3);
    Semi_minor = ellipse_j(4);
    Phi = ellipse_j(5);
    th = linspace(0, 1.5*pi, 100); %0 90 180 270
    x=cx+Semi_major*cos(Phi)*cos(th)-Semi_minor*sin(Phi)*sin(th);
    y=cy+Semi_minor*cos(Phi)*sin(th)+Semi_major*sin(Phi)*cos(th);
    p_j=[x; y]';
    N = size(p_j,1);
    
    % 平均得分
%     s = Score(ellipse_i,p_j); 
%     H = s/N;
%     if H >= 0.9
%         flag = true;
%     end
    
    % 平均距离
    c=ellipse_i(1:2)';
    axis=ellipse_i(3:4)';
    theta=ellipse_i(5);
    Maha_Dist = MahalaDist(p_j',c,axis,theta,1);
    H = sum(Maha_Dist)/N;
    if H <= Mahala_tolerance
        flag = true;
    end
    
end

