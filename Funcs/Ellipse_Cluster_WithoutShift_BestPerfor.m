function [ candidates_ ] = Ellipse_Cluster_WithoutShift_BestPerfor( candidates, points, normals, candidates_cover, Tmin, sigma_L, sigma_w  )
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
    % [y,x]=find(edge);
    % points=[x,y];
    % get every candidates' inlier
    inliers_matrix=cell(1,size(candidates,1));
    for i=1:size(candidates,1)
        %ellipse circumference is approximate pi * (1.5*sum(ellipseAxes)-sqrt(ellipseAxes(1)*ellipseAxes(2))
        ellipseCenter = candidates(i, 1 : 2);
        ellipseAxes_i   = candidates(i, 3:4);
        % 取180与周长*Tr中的小者
        tbins = min( [ 180, floor( pi * (1.5*sum(ellipseAxes_i)-sqrt(ellipseAxes_i(1)*ellipseAxes_i(2)) ) * Tmin ) ] );%选分区 default Tmin=0.6
        %ellipse_normals = computePointAngle(candidates(i,:),points);
        %inliers = find( labels == 0 & dRosin_square(candidates(i,:),points) <= 1 );  % +-1个像素内找支持内点
        %加速计算，只挑出椭圆外接矩形内的边缘点(椭圆中的长轴a>b),s_dx存储的是相对points的索引
        %points中的是edge point，s_dx存的是符合条件的点的索引值
        s_dx = find( points(:,1) >= (ellipseCenter(1)-ellipseAxes_i(1)-2) & points(:,1) <= (ellipseCenter(1)+ellipseAxes_i(1)+2) & points(:,2) >= (ellipseCenter(2)-ellipseAxes_i(1)-2) & points(:,2) <= (ellipseCenter(2)+ellipseAxes_i(1)+2));
        inliers = s_dx(dRosin_square(candidates(i,:),points(s_dx,:)) <= 1);    % 距离限制
        ellipse_normals = computePointAngle(candidates(i,:),points(inliers,:));     % 计算points点在ellipse上的法线
        p_dot_temp = dot(normals(inliers,:), ellipse_normals, 2); % 点乘表示夹角 加速后ellipse_normals(inliers,:)改为加速后ellipse_normals
        p_cnt = sum(p_dot_temp>0);%无奈之举
        if(p_cnt > size(inliers,1)*0.5) 
            %极性相异,也就是内黑外白     挑选极性相同的inliers，并且满足角度限制22.5°之内
            %ellipse_polarity = -1;
            inliers = inliers(p_dot_temp>0 & p_dot_temp >= 0.923879532511287 );%cos(pi/8) = 0.923879532511287, 夹角小于22.5°  
        else
            %极性相同,也就是内白外黑 
            %ellipse_polarity = 1;
            inliers = inliers(p_dot_temp<0 & (-p_dot_temp) >= 0.923879532511287 );
        end
        inliers = inliers(takeInliers(points(inliers, :), ellipseCenter, tbins));%论文中的SI(Lj)>Length(Lj)
        inliers = inliers';
        inliers_matrix{i} = points(inliers,:);
    end
    %% 寻找是否有同源椭圆
    candidates_label = true(size(candidates,1),1); %
    Homo_label = true(size(candidates,1),1);
    dis_tolerance = 3;
    for i = 1:size(candidates,1)-1
        if ~Homo_label(i,:)
            continue;
        end
        reference = i;
        ellipseCenter_i = candidates(reference,1:2);
        ellipseAxes_i = candidates(reference,3:4);
        theta_i = candidates(reference,5);
        num_Homologous = 1;
        Homo_vector = zeros(size(candidates,1),1);
        % j从i后面开始遍历
        for j = i+1:size(candidates,1)
            ellipseCenter_j = candidates(j,1:2);
            ellipseAxes_j = candidates(j,3:4);
            theta_j = candidates(j,5);
%             distance_center = norm(ellipseCenter_i-ellipseCenter_j);
%             distance_axis   = norm(ellipseAxes_i-ellipseAxes_j);
%             distance_theta  = norm(theta_i-theta_j);
            distance_center = sqrt((ellipseCenter_i(1)-ellipseCenter_j(1))^2+(ellipseCenter_i(2)-ellipseCenter_j(2)));
            distance_axis   = sqrt((ellipseAxes_i(1)-ellipseAxes_j(1))^2+(ellipseAxes_i(2)-ellipseAxes_j(2)));
            distance_theta  = norm(theta_i-theta_j);
            
            if distance_center > dis_tolerance || distance_axis > dis_tolerance || distance_theta > dis_tolerance*pi/180% 如果距离大于3就将j加入到Homo_vector
                continue;
            end
            flag = Judge_Ellipse(candidates(reference,:), inliers_matrix{j}); % 判断jth candidate是否与ith同源
            if ~flag
                continue;
            end
            Homo_vector(num_Homologous,:) = j;  % Homo_vector中保存通过与i距离检测的candidates索引
            num_Homologous = num_Homologous + 1;
            candidates_label(j,:) = false;  % 标记为false表示不输出
        end
        if ~isempty(find(Homo_vector~=0, 1)) 
            Homo_vector(num_Homologous,:) = reference;  % Homo_vector中保存同源candidates的序号
            Homo_label(Homo_vector(1:num_Homologous,:)) = false;  % 将同源的几组candidates的label标记为false，不参与后续的i循环
            candidates_label(reference,:) = false;
            % 检测Homo_vector中谁更可靠，返回最可靠的candidate序号num_candidates
            num_candidates = Computr_Reliability( Homo_vector(1:num_Homologous,:), inliers_matrix,...
                                                  candidates, candidates_cover, sigma_L, sigma_w );
            candidates_label(Homo_vector(num_candidates,:),:) = true;  % 标记为true表示输出该candidate
        end
    end
    candidates_ = candidates( candidates_label,: )';
end

%%
function [dmin]= dRosin_square(param,points)
ae2 = param(3).*param(3);
be2 = param(4).*param(4);
x = points(:,1) - param(1);
y = points(:,2) - param(2);
xp = x*cos(-param(5))-y*sin(-param(5));
yp = x*sin(-param(5))+y*cos(-param(5));
fe2 = ae2-be2;
X = xp.*xp;
Y = yp.*yp;
delta = (X+Y+fe2).^2-4*fe2*X;
A = (X+Y+fe2-sqrt(delta))/2;
ah = sqrt(A);
bh2 = fe2-A;
term = A*be2+ae2*bh2;
xi = ah.*sqrt(ae2*(be2+bh2)./term);
yi = param(4)*sqrt(bh2.*(ae2-A)./term);
d = zeros(size(points,1),4);%n x 4
d(:,1) = (xp-xi).^2+(yp-yi).^2;
d(:,2) = (xp-xi).^2+(yp+yi).^2;
d(:,3) = (xp+xi).^2+(yp-yi).^2;
d(:,4) = (xp+xi).^2+(yp+yi).^2;
dmin = min(d,[],2); %返回距离的平方
%[dmin, ii] = min(d,[],2); %返回距离的平方
% for jj = 1:length(dmin)
%     if(ii(jj) == 1)
%         xi(jj) = xi(jj);
%         yi(jj) = yi(jj);
%     elseif (ii(jj) == 2)
%         xi(jj) = xi(jj);
%         yi(jj) = -yi(jj);
%     elseif (ii(jj) == 3)
%         xi(jj) = -xi(jj);
%         yi(jj) = yi(jj);
%     elseif(ii(jj) == 4)
%          xi(jj) = -xi(jj);
%         yi(jj) = -yi(jj);
%     end
% end
% 
% xi =  xi*cos(param(5))-yi*sin(param(5));
% yi =  xi*sin(param(5))+yi*cos(param(5));
% 
% testim = zeros(300,300);
% testim(sub2ind([300 300],uint16(yi+param(2)),uint16(xi+param(1)))) = 1;
% figure;imshow(uint8(testim).*255);
end
%% compute the points' normals belong to an ellipse, the normals have been already normalized. 
%param: [x0 y0 a b phi].
%points: [xi yi], n x 2
function [ellipse_normals] = computePointAngle(ellipse, points)
%convert [x0 y0 a b phi] to Ax^2+Bxy+Cy^2+Dx+Ey+F = 0
a_square = ellipse(3)^2;
b_square = ellipse(4)^2;
sin_phi = sin(ellipse(5));
cos_phi = cos(ellipse(5)); 
sin_square = sin_phi^2;
cos_square = cos_phi^2;
A = b_square*cos_square+a_square*sin_square;
B = (b_square-a_square)*sin_phi*cos_phi*2;
C = b_square*sin_square+a_square*cos_square;
D = -2*A*ellipse(1)-B*ellipse(2);
E = -2*C*ellipse(2)-B*ellipse(1);
% F = A*ellipse(1)^2+C*ellipse(2)^2+B*ellipse(1)*ellipse(2)-(ellipse(3)*ellipse(4)).^2;
% A = A/F;
% B = B/F;
% C = C/F;
% D = D/F;
% E = E/F;
% F = 1;
%calculate points' normals to ellipse
angles = atan2(C*points(:,2)+B/2*points(:,1)+E/2, A*points(:,1)+B/2*points(:,2)+D/2);
ellipse_normals = [cos(angles),sin(angles)];
end
%%
function idx = takeInliers(x, center, tbins)
    [theta, ~] = cart2pol(x(:, 1) - center(1), x(:, 2) - center(2));%得到[-pi,pi]的方位角，等价于 theta = atan2(x(:, 2) - center(2) , x(:, 1) - center(1)); 
    % theta中保存每个内点关于椭圆中心点的角度
    tmin = -pi; tmax = pi;
    % 将theta中每个点原本在[-180°,180°]内，现在将其映射到[1，tbins]内。(+0.5是什么意思)
    tt = round((theta - tmin) / (tmax - tmin) * tbins + 0.5);%将内点分区到[1 tbins]
    tt(tt < 1) = 1; 
    tt(tt > tbins) = tbins;
    h = histc(tt, 1 : tbins);%h为直方图[1 tbins]的统计结果
    mark = zeros(tbins, 1);
    compSize = zeros(tbins, 1);
    nComps = 0;
    queue = zeros(tbins, 1);
    du = [-1, 1];
    for i = 1 : tbins  %tbins个区域
        if (h(i) > 0 && mark(i) == 0)%如果落在第i个分区内的值大于0，且mark(i)为0
            nComps = nComps + 1;
            mark(i) = nComps;%mark(1)=1标记第nComps个连通区域
            front = 1; rear = 1;
            queue(front) = i;%queue(1) = 1将该分区加入队列，并以此开始任务
            while (front <= rear)  %11，23,
                u = queue(front);  % queue(2) = 180
                front = front + 1; %2，3
                for j = 1 : 2
                    v = u + du(j);  %0,2,179
                    if (v == 0)
                        v = tbins; %180
                    end
                    if (v > tbins)
                        v = 1;
                    end
                    if (mark(v) == 0 && h(v) > 0)
                        rear = rear + 1;%2,3,4
                        queue(rear) = v;%queue(2) = 180,queue(3) = 2,queue(4) = 179
                        mark(v) = nComps;%mark(179) = 1,mark(2) = 2标记第nComps个连通区域
                    end
                end
            end
            compSize(nComps) = sum(ismember(tt, find(mark == nComps)));%得到构成连通域为nComps的内点数量
        end
    end
    compSize(nComps + 1 : end) = [];
    maxCompSize = max(compSize);
    validComps = find(compSize >= maxCompSize * 0.1 & compSize > 10);%大于等于最大连通长度的0.1倍的连通区域是有效的
    validBins = find(ismember(mark, validComps));%有效的分区
    idx = ismember(tt, validBins);%有效的内点
end

%%
function num = Computr_Reliability( Homo_vector, inliers_matrix, candidates, candidates_cover,sigma_L,sigma_w)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
%     candidates = candidates';
    for i=1:size(Homo_vector,1)
        p_i = inliers_matrix{Homo_vector(i)};
        candidates_i = candidates(Homo_vector(i),:);
        c_i = candidates_i(:,1:2);
        semi_i = candidates_i(:,3:4);
        theta_i = candidates_i(:,5);
        f_i = sqrt(power(semi_i(1),2)-power(semi_i(2),2));
        f1_i_t = c_i + [cos(theta_i), -sin(theta_i)]*f_i;
        f2_i_t = c_i + [-cos(theta_i), sin(theta_i)]*f_i;
        f1_i = repmat(f1_i_t,size(p_i,1),1);
        f2_i = repmat(f2_i_t,size(p_i,1),1);
        % 计算距离Mi
        temp = p_i - f1_i;
        dis = sqrt(temp(:,1).^2 + temp(:,2).^2);
        temp = p_i - f2_i;
        dis = dis + sqrt(temp(:,1).^2 + temp(:,2).^2) - 2*f_i;
        Mi = sum(dis)/size(p_i,1);
        sigma_wj = 0;
        sigma_gj = 0;
        sigma_wgL = 0;
        MD_i = sum(MahalaDist(p_i',c_i',semi_i,theta_i,1))/size(p_i,1);
        
        for j=1:size(Homo_vector,1)
            if i==1
                continue;
            end
            p_j = inliers_matrix{Homo_vector(j)};
            candidates_j = candidates(Homo_vector(j),:);
            c_j = candidates_j(:,1:2);            
            semi_j = candidates_j(:,3:4);            
            theta_j = candidates_j(:,5);            
            f_j = sqrt(power(semi_j(1),2)-power(semi_j(2),2));            
            f1_j_t = c_j + [cos(theta_j), -sin(theta_j)]*f_j;
            f2_j_t = c_j - [-cos(theta_j), sin(theta_j)]*f_j;
            f1_j = repmat(f1_j_t,size(p_j,1),1);
            f2_j = repmat(f2_j_t,size(p_j,1),1);
            temp = p_j - f1_j;
            dis = sqrt(temp(:,1).^2 + temp(:,2).^2);
            temp = p_j - f2_j;
            dis = dis + sqrt(temp(:,1).^2 + temp(:,2).^2) - 2*semi_j(1);
            Mj = sum(dis)/size(p_j,1);
            MD_j = sum(MahalaDist(p_j',c_i',semi_i,theta_i,1))/size(p_j,1);
            MD_ij = abs(MD_j-MD_i);
%             tolerance = 
            Ellipse_error = abs(Mj-Mi);
            if Ellipse_error < 5
%                 L_j = exp( -(Mj-Mi)^2 / 2*(sigma_L^2) );
                L_j = exp( (-(MD_ij/MD_i)^2) / (2*(sigma_L^2)) );
            else
                L_j = 0;
            end
            n = size(Homo_vector,1);  % n是否正确
            w_j = 1/(sqrt(2*pi)*sigma_w) * exp( -((j-1)/n)^2 / 2*(sigma_L^2) );
            g_j = numel(p_j);
            sigma_wj = sigma_wj + w_j;
            sigma_gj = sigma_gj + g_j;
            sigma_wgL = sigma_wgL + (w_j*g_j*L_j);
        end
        temp_ = sigma_wgL / (sigma_wj*sigma_gj);
        temp_1 = sqrt(temp_^2 * (candidates_cover(i)/2*pi)^2);
        if i==1
            num = i;
            reliability_i = temp_1;
        else
            if reliability_i < temp_1
                num = i;
                reliability_i = temp_1;
            end
        end
    end
end

% 马氏距离
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