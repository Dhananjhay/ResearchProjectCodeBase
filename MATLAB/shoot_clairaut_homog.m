function [f,xx,yy] = shoot_clairaut_homog(v,mm,rhom,rhof,xcmb)
% Shooting routine to be used in Clairaut's equation
% from x=xcmb to 1, (xcmb = cmb)

% to solve in mantle y1' = y2
%                    y2' = -6/x * z * y2 - 6/x^2 *(z-1) y1
%               where z = rhom / (rhom + (rhof-rhom)*(xcmb/x)^3)
%
% with b.c. y2(xcmb)=0, y2(1)+2y1(1)=2.5*mm  where mm = rot^2 * R^3 / G M


% Shoot from cmb to 1 (in mantle)
% initial condition
y0 = [v;0.0];

[xx,yy] = ode45('rhs_clairaut',[xcmb 1.0],y0,[],xcmb,rhof,rhom);
nx=length(xx);


% function to minimize: y2 at R = 1 
f = yy(nx,2) + 2.*yy(nx,1) - 2.5*mm;


end

