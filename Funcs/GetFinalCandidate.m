function [elp] = GetFinalCandidate(elps_set,candidates, points, normals, Mahala_tolerance, Tmin)
%UNTITLED 返回当前同源集合的聚类结果
%   此处显示详细说明
    elp_num = size(elps_set,1);
    Homo_inliers = [];
    for i=1:elp_num
        E = candidates(elps_set(i),:);
%         ellipseCenter = E(1 : 2);
        ellipseAxes = E(3:4);
        tbins = min( [ 180, floor( pi * (1.5*sum(ellipseAxes)-sqrt(ellipseAxes(1)*ellipseAxes(2)) ) * Tmin ) ] );
        s_dx = find( points(:,1) >= (E(1)-E(3)-1) & points(:,1) <= (E(1)+E(3)+1) & points(:,2) >= (E(2)-E(3)-1) & points(:,2) <= (E(2)+E(3)+1));
        inliers_M = s_dx(ASI_Dist(points(s_dx,:)',E(1:2)',E(3:4)',E(5),1) <= 1.8*Mahala_tolerance);
        inliers = inliers_M(dRosin_square(E,points(inliers_M,:)) <= 1);
        ellipse_normals = computePointAngle(E,points(inliers,:));
        p_dot_temp = dot(normals(inliers,:), ellipse_normals, 2);
        p_cnt = sum(p_dot_temp>0);
        if(p_cnt > size(inliers,1)*0.5)
            inliers = inliers(p_dot_temp>0 & p_dot_temp >= 0.923879532511287 );
        else
            inliers = inliers(p_dot_temp<0 & (-p_dot_temp) >= 0.923879532511287 );
        end
        inliers = inliers(takeInliers(points(inliers, :), E(1:2), tbins));
        
        Homo_inliers = union(Homo_inliers,inliers);
%         Homo_inliers = intersect(Homo_inliers,inliers);
    end
    [new_ellipse,new_info] = fitEllipse(points(Homo_inliers,1),points(Homo_inliers,2));
    if new_info
        elp = new_ellipse;
    else
        elp = candidates(elps_set(1),:);
    end
end

