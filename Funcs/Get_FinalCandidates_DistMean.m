function [final_candiates] = Get_FinalCandidates_MahalaDistMean(candidates, elps, points,normals,Homo_tolerance,Tmin)
%判断同源集合中的最大可靠elp
%   输入：
%   elps：候选椭圆
%   输出：
%   final_candidates：每一类的最佳候选
homosets_num = size(elps,1);
elp = zeros(homosets_num,5);
for i=1:homosets_num
    elps_set = elps{i};
    elps_num = size(elps_set,1);
    % 如果集合中只有一个elp则自成一类,如果集合中超过1个则计算可靠度
    if elps_num == 1
        elp(i,:) = candidates(elps{i},:);
        continue;
    end
    % get final ellpise candidate
%     final_candiates(i,1) = GetMostReliable(candidates, elps_set, sigma_L, sigma_w, Ellipse_error_tolerance, candidates_cover);
%     final_candiates(i,1) = elps_set(1); % used
    elp(i,:) = GetFinalCandidate(elps_set, candidates, points, normals, Homo_tolerance, Tmin);
end
final_candiates = elp;
end


