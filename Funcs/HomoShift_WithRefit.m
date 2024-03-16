function [MostReliableNum, HomoLabel, HomoVector] = HomoShift_WithRefit( candidates, reference, inliers_matrix, Homo_label, ...
                                                   candidates_cover, sigma_L, sigma_w, dis_tolerance )
%UNTITLED2 此处显示有关此函数的摘要
%   输入:
%   candidates: 候选椭圆
%   reference: 当前最佳椭圆索引
%   inliers_matrix: 各candidates的内点
%   Homo_label: candidates的同源使用情况
%   输出:
%   MostReliableNum: 最佳可靠候选索引
%   HomoLabel: 更新后的candidates的同源使用情况
    ellipseCenter_ref = candidates(reference,1:2);
    ellipseAxes_ref = candidates(reference,3:4);
    theta_ref = candidates(reference,5);
    num_Homologous = 1;
    temp_Homo_label = Homo_label;
    num_candidates = reference;
    % 找同源子集
    for i=1:size(candidates,1)
        if ~temp_Homo_label(i,:)
            continue;
        end
        ellipseCenter_que = candidates(i,1:2);
        ellipseAxes_que = candidates(i,3:4);
        theta_que = candidates(i,5);
        distance_center = norm(ellipseCenter_ref-ellipseCenter_que);
        distance_axis   = norm(ellipseAxes_ref-ellipseAxes_que);
        distance_theta  = norm(theta_ref-theta_que);
        if distance_center > dis_tolerance || distance_axis > dis_tolerance || distance_theta > dis_tolerance*pi/180% 如果距离大于3就将j加入到Homo_vector
            continue;
        end
        flag = Judge_Ellipse(candidates(reference,:), inliers_matrix{i}); % 判断jth candidate是否与ith同源
        if ~flag
            continue;
        end
        Homo_vector(num_Homologous,:) = i;  % Homo_vector中保存通过与i同源判断的candidates索引
        num_Homologous = num_Homologous + 1;
        temp_Homo_label(i,:) = false;  % 标记为false表示与其他同源
    end
    
    if num_Homologous == 1
        temp_Homo_label(reference,:) = false;
        MostReliableNum = reference;
        HomoLabel = temp_Homo_label;
        HomoVector(num_Homologous,:) = reference;
        return
    else
        num_candidates = Computr_Reliability( Homo_vector, inliers_matrix,...
                                          candidates, candidates_cover, sigma_L, sigma_w );
        if num_candidates ~= reference
            % 递归shift找到最佳代表
            reference = num_candidates;
            [MostReliableNum, HomoLabel,HomoVector] = HomoShift( candidates, reference, inliers_matrix, Homo_label, ...
                candidates_cover, sigma_L, sigma_w, dis_tolerance );
            return
        else
            MostReliableNum = reference;
            HomoLabel = temp_Homo_label;
            HomoVector = Homo_vector;
            return
        end
    end
    
end

% 判断两椭圆同源
function flag = Judge_Ellipse(ellipse_i, p_j)
    flag = false;
    N = size(p_j,1);
    s = Score(ellipse_i,p_j);
    H = s/N;
    if H >= 0.98
        flag = true;
    end
end

% score函数
function [s] = Score(e,p)
    c=e(1:2)';
    axis=e(3:4)';
    theta=e(5);
    t = 3;
    tol_p = c+[cos(theta)*(axis(1)+t); sin(theta)*(axis(1)+t)];
    % 计算椭圆e长轴+t像素上点的D4作为容忍度阈值。
    tolerance = MahalaDist(tol_p,c,axis,theta,1);
    Maha_Dist = MahalaDist(p',c,axis,theta,1);
    s_index = find(Maha_Dist<tolerance);
    s=size(s_index,1);
end

function Maha_Dist = MahalaDist(x,c,axis,theta,r)
% 输入：x:[x1,y1;x2,y2;...xn,yn], c:椭圆中心[cx,cy], axis:[a,b], r: 马氏距离中圆的半径
    U = [cos(theta), -sin(theta); sin(theta), cos(theta)];
    T = [axis(1)^2, 0; 0, axis(2)^2];
    S = U * T * U';
%     Maha_Dist = abs( sqrt((x-c)' * inv(S) * (x-c)) - r ); %D1
%     Maha_Dist = power( sqrt((x-c)' * inv(S) * (x-c)) - r,2 ); %D2
%     Maha_Dist = power( ((x(1,:)-c)' * inv(S) * (x(1,:)-c) - r),2 ); %D4
%     Maha_Dist=zeros(size(x,1),1);
%     for i=1:size(x,1)
%         Maha_Dist(i)=power( ((x(i,:)-c) * inv(S) * (x(i,:)-c)' - r),2 );
%     end
    Maha_Dist=diag(power( ((x-c)' * inv(S) * (x-c) - r),2 ),0);
end