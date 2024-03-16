function p_dot_temp = point_ds(elp, p_ns, e_ns)
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
semi_major = elp(3);
semi_minor = elp(4);
ratio = semi_major / semi_minor;

p_ns_th = atan(p_ns(:,2) ./ p_ns(:,1));
e_ns_th = atan(e_ns(:,2) ./ e_ns(:,1));

p_ns_th_EN = atan((p_ns(:,2)/ratio) ./ p_ns(:,1));
e_ns_th_EN = atan((e_ns(:,2)/ratio) ./ e_ns(:,1));

p_ns_EN = [cos(p_ns_th_EN), sin(p_ns_th_EN)];
e_ns_EN = [cos(e_ns_th_EN), sin(e_ns_th_EN)];

p_dot_temp = dot(p_ns_EN, e_ns_EN, 2);
end