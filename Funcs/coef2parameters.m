function [ellipse] = coef2parameters(par)
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
thetarad = 0.5*atan2(par(2),par(1) - par(3));
cost = cos(thetarad);
sint = sin(thetarad);
sin_squared = sint.*sint;
cos_squared = cost.*cost;
cos_sin = sint .* cost;
Ao = par(6);
Au = par(4) .* cost + par(5) .* sint;
Av = - par(4) .* sint + par(5) .* cost;
Auu = par(1) .* cos_squared + par(3) .* sin_squared + par(2) .* cos_sin;
Avv = par(1) .* sin_squared + par(3) .* cos_squared - par(2) .* cos_sin;
% ROTATED = [Ao Au Av Auu Avv]
tuCentre = - Au./(2.*Auu);
tvCentre = - Av./(2.*Avv);
wCentre = Ao - Auu.*tuCentre.*tuCentre - Avv.*tvCentre.*tvCentre;
uCentre = tuCentre .* cost - tvCentre .* sint;
vCentre = tuCentre .* sint + tvCentre .* cost;
Ru = -wCentre./Auu;
Rv = -wCentre./Avv;
Ru = sqrt(abs(Ru));
Rv = sqrt(abs(Rv));
ellipse = [uCentre, vCentre, Ru, Rv, thetarad];
 %会出现Ru < Rv情况，对调一下
if(Ru < Rv )
   ellipse(3) = Rv;
   ellipse(4) = Ru;
   if(thetarad < 0)
     ellipse(5) = ellipse(5)+1.570796326794897; %pi/2
   else
     ellipse(5) = ellipse(5)-1.570796326794897;
   end
end
end

