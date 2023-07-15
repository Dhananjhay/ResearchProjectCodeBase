function [J,f] = J_exo_epsm(x)
%  computes the Jacobian and function evaluation for Large angle system
% 
% input x = variable x (3 element vector)
% output f = function evaluation (2 function)
%        J = Jacobian matrix of derivatives

% x(1) = theta_m
% x(2) = theta_f
% x(3) = epsilon_m

global A Af ell ellf dsigma  ...
       pot3 inclI;

% initialize
n=length(x);
f = zeros(n,1); % f must be defined as a column vector
J = zeros(n,n);  

% compute the function f
 
s=-1.0-cos(inclI+x(3))*dsigma; % definition of sigma
 
 sIem=sin(inclI +x(3));
 cIem=cos(inclI +x(3));
 x12=x(1) + x(2);
 s12= sin(x12) - sin(x(1));
 c12= cos(x12) - cos(x(1));
 
 
 % compute f
 
 f1in = ( s - ell*cos(x(1)) ) *sin(x(1) ) ...
     + (Af/A)*( s*s12 + sin(x(2)) - ellf*sin(x(1))*c12 );
  % part that multiplies sigma, useful for J
 fs1 = sin(x(1))+ (Af/A)*s12 ;
 f(1) = f1in + pot3*ell*sin(x(3))*cos(x(3));  
     
 f(2) = ( sin(x(2)) + s*sin(x12) + ellf*cos(x12)*s12 );
 fs2 =  sin(x12)  ;
 
 
 f(3) = s*sin(inclI+x(3)) + sin(x(1)+inclI+x(3)) ;

% compute the Jacobian matrix
 J(1,1) = ( s*cos(x(1)) - ell*cos(x(1))^2 +ell*sin(x(1))^2  ) ...
        + Af/A*( s*c12 -ellf*cos(x(1))*c12 + ellf*sin(x(1))*s12 );
 J(1,2) = Af/A*(s*cos(x12) + cos(x(2))) ...
        + Af/A*(ellf*sin(x(1))*sin(x12) )  ;
 J(1,3) = pot3*ell*( cos(x(3))^2 - sin(x(3))^2 )  ...
        + sIem*dsigma*fs1;

    
 J(2,1) = ( s*cos(x12) + ellf*( cos(x12)*c12 -sin(x12)*s12) );
 J(2,2) = (cos(x(2)) + s*cos(x12) ...
              +ellf*(cos(x12)^2 -sin(x12)*s12) );
 J(2,3) = sIem*dsigma*fs2;
 
 
 
 J(3,1) = cos(x(1)+inclI+x(3));
 J(3,2) = 0.0;
 J(3,3) = s*cIem + cos(x(1)+inclI+x(3)) + sIem^2*dsigma;
 
 
end

