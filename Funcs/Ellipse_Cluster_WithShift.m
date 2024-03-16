function [ candidates_ ] = Ellipse_Cluster_WithShift( candidates, points, normals, candidates_cover, Tmin, sigma_L, sigma_w  )
%UNTITLED �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
    % [y,x]=find(edge);
    % points=[x,y];
    % get every candidates' inlier
    inliers_matrix=cell(1,size(candidates,1));
    for i=1:size(candidates,1)
        %ellipse circumference is approximate pi * (1.5*sum(ellipseAxes)-sqrt(ellipseAxes(1)*ellipseAxes(2))
        ellipseCenter = candidates(i, 1 : 2);
        ellipseAxes_ref   = candidates(i, 3:4);
        % ȡ180���ܳ�*Tr�е�С��
        tbins = min( [ 180, floor( pi * (1.5*sum(ellipseAxes_ref)-sqrt(ellipseAxes_ref(1)*ellipseAxes_ref(2)) ) * Tmin ) ] );%ѡ���� default Tmin=0.6
        %ellipse_normals = computePointAngle(candidates(i,:),points);
        %inliers = find( labels == 0 & dRosin_square(candidates(i,:),points) <= 1 );  % +-1����������֧���ڵ�
        %���ټ��㣬ֻ������Բ��Ӿ����ڵı�Ե��(��Բ�еĳ���a>b),s_dx�洢�������points������
        %points�е���edge point��s_dx����Ƿ��������ĵ������ֵ
        s_dx = find( points(:,1) >= (ellipseCenter(1)-ellipseAxes_ref(1)-2) & points(:,1) <= (ellipseCenter(1)+ellipseAxes_ref(1)+2) & points(:,2) >= (ellipseCenter(2)-ellipseAxes_ref(1)-2) & points(:,2) <= (ellipseCenter(2)+ellipseAxes_ref(1)+2));
        inliers = s_dx(dRosin_square(candidates(i,:),points(s_dx,:)) <= 1);    % ��������
        ellipse_normals = computePointAngle(candidates(i,:),points(inliers,:));     % ����points����ellipse�ϵķ���
        p_dot_temp = dot(normals(inliers,:), ellipse_normals, 2); % ��˱�ʾ�н� ���ٺ�ellipse_normals(inliers,:)��Ϊ���ٺ�ellipse_normals
        p_cnt = sum(p_dot_temp>0);%����֮��
        if(p_cnt > size(inliers,1)*0.5) 
            %��������,Ҳ�����ں����     ��ѡ������ͬ��inliers����������Ƕ�����22.5��֮��
            %ellipse_polarity = -1;
            inliers = inliers(p_dot_temp>0 & p_dot_temp >= 0.923879532511287 );%cos(pi/8) = 0.923879532511287, �н�С��22.5��  
        else
            %������ͬ,Ҳ�����ڰ���� 
            %ellipse_polarity = 1;
            inliers = inliers(p_dot_temp<0 & (-p_dot_temp) >= 0.923879532511287 );
        end
        inliers = inliers(takeInliers(points(inliers, :), ellipseCenter, tbins));%�����е�SI(Lj)>Length(Lj)
        inliers = inliers';
        inliers_matrix{i} = points(inliers,:);
    end
    %% Ѱ���Ƿ���ͬԴ��Բ
    candidates_label = false(size(candidates,1),1);  %
    Homo_label = true(size(candidates,1),1);
    Homo_vector = zeros(size(candidates,1),1);
    dis_tolerance = 3;
    for i = 1:size(candidates,1)-1
        if ~Homo_label(i,:)
            continue;
        end
        Homo_label(i,:) = false;  % �� ith candidate��Ϊ��ʹ�ã���Ϊ�����Ƿ���ͬԴ��Բ�������
        reference = i;
        num_Homologous = 1;
        Homo_vector(num_Homologous,:) = reference;  % reference����Homo_vector
        num_Homologous = num_Homologous + 1;
        ShiftFlag = true;
        % �ҵ�ͬԴ����
        while ShiftFlag
            Ellipse_ref = candidates(reference,:);  % ����reference center
            temp_Homo_label = Homo_label;  % ��ʱ��Ǻ�ѡ��ʹ�����
            for j = i+1:size(candidates,1)
                if ~temp_Homo_label(j,:)
                    continue;
                end
                Ellipse_que = candidates(j,:);
                distance_center = norm( Ellipse_ref(1:2) - Ellipse_que(1:2) );
                distance_axis   = norm( Ellipse_ref(3:4) - Ellipse_que(3:4) );
                distance_theta  = norm( Ellipse_ref(5) - Ellipse_que(5) );
                if distance_center > dis_tolerance || distance_axis > dis_tolerance || distance_theta > dis_tolerance*pi/180% ����������3�ͽ�j���뵽Homo_vector
                    continue;
                end
                % �ж�jth candidate�Ƿ���ithͬԴ
                flag = Judge_EllipseInstance(Ellipse_ref, Ellipse_que);
%                 flag = Judge_Ellipse(Ellipse_ref, inliers_matrix{j}); 
                if ~flag
                    continue;
                end
                Homo_vector(num_Homologous,:) = j;  % Homo_vector�б���ͨ����iͬԴ�жϵ�candidates����
                num_Homologous = num_Homologous + 1;
                temp_Homo_label(j,:) = false;  % ���Ϊfalse��ʾ������ͬԴ
            end
            if num_Homologous == 2  % ith candidate û��ͬԴ���ϣ����Լ�ֱ�����
                MostReliableNum = reference;
                ShiftFlag = false;
                Homo_label = temp_Homo_label;
                candidates_label(MostReliableNum, :) = true;
            else
                 % ���㵱ǰHomo_vector�е����
                 num_candidates = Computr_Reliability( Homo_vector(1:num_Homologous-1), inliers_matrix, candidates,...
                     candidates_cover, sigma_L, sigma_w );
                 if num_candidates == reference
                     MostReliableNum = reference;
                     ShiftFlag = false;
                     Homo_label = temp_Homo_label;  % ����Homo_label,�Ѿ�ʹ�ù���candidates�������������
                     candidates_label(MostReliableNum, :) = true;  % ֻ�����ɿ���
                 else
                     reference = num_candidates;
                 end
            end
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
dmin = min(d,[],2); %���ؾ����ƽ��
%[dmin, ii] = min(d,[],2); %���ؾ����ƽ��
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
    [theta, ~] = cart2pol(x(:, 1) - center(1), x(:, 2) - center(2));%�õ�[-pi,pi]�ķ�λ�ǣ��ȼ��� theta = atan2(x(:, 2) - center(2) , x(:, 1) - center(1)); 
    % theta�б���ÿ���ڵ������Բ���ĵ�ĽǶ�
    tmin = -pi; tmax = pi;
    % ��theta��ÿ����ԭ����[-180��,180��]�ڣ����ڽ���ӳ�䵽[1��tbins]�ڡ�(+0.5��ʲô��˼)
    tt = round((theta - tmin) / (tmax - tmin) * tbins + 0.5);%���ڵ������[1 tbins]
    tt(tt < 1) = 1; 
    tt(tt > tbins) = tbins;
    h = histc(tt, 1 : tbins);%hΪֱ��ͼ[1 tbins]��ͳ�ƽ��
    mark = zeros(tbins, 1);
    compSize = zeros(tbins, 1);
    nComps = 0;
    queue = zeros(tbins, 1);
    du = [-1, 1];
    for i = 1 : tbins  %tbins������
        if (h(i) > 0 && mark(i) == 0)%������ڵ�i�������ڵ�ֵ����0����mark(i)Ϊ0
            nComps = nComps + 1;
            mark(i) = nComps;%mark(1)=1��ǵ�nComps����ͨ����
            front = 1; rear = 1;
            queue(front) = i;%queue(1) = 1���÷���������У����Դ˿�ʼ����
            while (front <= rear)  %11��23,
                u = queue(front);  % queue(2) = 180
                front = front + 1; %2��3
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
                        mark(v) = nComps;%mark(179) = 1,mark(2) = 2��ǵ�nComps����ͨ����
                    end
                end
            end
            compSize(nComps) = sum(ismember(tt, find(mark == nComps)));%�õ�������ͨ��ΪnComps���ڵ�����
        end
    end
    compSize(nComps + 1 : end) = [];
    maxCompSize = max(compSize);
    validComps = find(compSize >= maxCompSize * 0.1 & compSize > 10);%���ڵ��������ͨ���ȵ�0.1������ͨ��������Ч��
    validBins = find(ismember(mark, validComps));%��Ч�ķ���
    idx = ismember(tt, validBins);%��Ч���ڵ�
end

%%
% function flag = Judge_Ellipse(ellipse_i, p_j)
%     flag = false;
%     N = size(p_j,1);
%     s = Score(ellipse_i,p_j);
%     H = s/N;
%     if H >= 0.98
%         flag = true;
%     end
% end

% function flag = Judge_EllipseInstance(ellipse_i, ellipse_j)
%     flag = false;
%     cx = ellipse_j(1);
%     cy = ellipse_j(2);
%     Semi_major = ellipse_j(3);
%     Semi_minor = ellipse_j(4);
%     Phi = ellipse_j(5);
%     
%     th = linspace(0, 2*pi, 6);
%     x=cx+Semi_major*cos(Phi)*cos(th)-Semi_minor*sin(Phi)*sin(th);
%     y=cy+Semi_minor*cos(Phi)*sin(th)+Semi_major*sin(Phi)*cos(th);
%     
%     p_j=[x; y]';
%     N = size(p_j,1);
%     s = Score(ellipse_i,p_j);
%     H = s/N;
%     if H >= 0.9
%         flag = true;
%     end
% end
% % score����
% function [s] = Score(e,p)
%     c=e(1:2)';
%     axis=e(3:4)';
%     theta=e(5);
%     t = 3;
%     tol_p = c+[cos(theta)*(axis(1)+t); sin(theta)*(axis(1)+t)];
%     % ������Բe����+t�����ϵ��D4��Ϊ���̶���ֵ��
%     tolerance = MahalaDist(tol_p,c,axis,theta,1);
%     rho = 0.005;
%     Maha_Dist = MahalaDist(p',c,axis,theta,1);
%     s_index = find(Maha_Dist<rho);
%     s=size(s_index,1);
% end

