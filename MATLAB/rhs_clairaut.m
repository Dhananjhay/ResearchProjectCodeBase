hfunction dydx = rhs_clairaut(x,y,flag,xcmb,rhof,rhom)
% rhs_clairaut  Right-hand sides of coupled ODEs for Clairaut's
% differeential equation

% Input:    x      = the independent variable 
%           y      = vector (length 1) of dependent variables
%           flag   = dummy argument for compatibility with ode45
%           
%
% Output:   dydx = column vector of dy(i)/dx values

if x<xcmb  % in core
  dydx = [ y(2);
           -6./x * y(2)];
else       % in mantle
  zeta = rhom / (rhom + (rhof-rhom)*(xcmb/x)^3);  
  dydx = [ y(2);
           -6./x * zeta * y(2)   - 6./x^2 * (zeta -1.0) *y(1)];
end
