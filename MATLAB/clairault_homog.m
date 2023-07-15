function [ellf,ell] = clairault_homog(mm,rhom,rhof,xcmb)

% function to find the flatening at cmb and surface by solving 
% Clairaut's differential equation.  Case for homogeneous density core &
% mantle


% Numerical solution

v0=0.001; % initial guess for f=y1 at x=0 (and x=xcmb)

% 1) Find value of v which matches b.c. 
v = fzero(@(x) shoot_clairaut_homog(x,mm,rhom,rhof,xcmb), v0);
% 2) final solution: recompute shoot, with the determined of v
[f,xx,yy] = shoot_clairaut_homog(v,mm,rhom,rhof,xcmb);

ellf=v;

flat=yy(end,1);
ell = (rhof*xcmb^5*ellf + rhom * (flat - ellf*xcmb^5) ) ...
       / (rhof*xcmb^5 + rhom * (1 - xcmb^5) );

% Plot solution (Non-dimensional)
%figure(10)
%subplot(2,1,1)
%plot(xx,yy(:,1),'-r');  
%xlabel('radius'); ylabel('y1'); 
%subplot(2,1,2)
%plot(xx,yy(:,2),'-r');  
%xlabel('radius'); ylabel('y2'); 


end

