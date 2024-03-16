function [ellipse] = coff2param(a)
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
A=a(1); B=a(2); C=a(3); D=a(4); E=a(5); F=a(6);
theta = 0.5*atan2(B ,(A-C));
xc = (B*E-2*C*D) / (4*A*C-B^2);
yc = (B*D-2*A*E) / (4*A*C-B^2);
major = 2*(A*xc^2+C*yc^2+B*xc*yc-F) / (A+C+sqrt((A-C)^2+B^2));
minor = 2*(A*xc^2+C*yc^2+B*xc*yc-F) / (A+C-sqrt((A-C)^2+B^2));
major = sqrt(abs(major));
minor = sqrt(abs(minor));
if major<minor
    t = major;
    major = minor;
    minor = t;
    if theta<0
        theta = theta+1.570796326794897; %pi/2
    else
        theta = theta-1.570796326794897; %pi/2
    end
end
ellipse = [xc, yc, major, minor, theta];
end

