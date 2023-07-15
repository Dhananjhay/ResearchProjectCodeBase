% solves for theta_m, theta_f and epslon_m 
% for an exoplanet in synchronous rotation for a given 
%        -- precession rate (Omega_p),
%        -- inclination I 
%        -- rotation rate (Omega_o = n = mean motion)

% Assumes a simple mantle+fluid core
% assumes tidally locked (rot rate = orbital rate)
% Assumes an axisymmetric planet
% we sweep through a range of Poincare number

clear all

global A Af ell ellf dsigma  ...
       pot3 inclI;
%
npts=101; % number of loop points for Poincaré number
nicond=50; % number of initial conditions for large angle sol

%constants
bigG=6.67384e-11;
MassSun=1.989e30;
RadSun=6.597e8;
J2Sun=1.e-7;
MassEarth=5.974e24;
oneau=1.496e11;
RadiusEarth=6.371e6;

MassStar=1.0*MassSun;

% forcing parameters
ecc=0.0; % eccentricity
inclI=(0.1)*pi/180; %orbital inclinaison
%pot1 = 1.5*sin(inclI)*cos(inclI)/((1.-ecc*ecc)^1.5); 
%pot2 = 1.5*((cos(inclI))^2 - (sin(inclI))^2)/((1.-ecc*ecc)^1.5); 
%pot1 = 1.5/((1.-ecc*ecc)^1.5); 
%pot2 = 1.5/((1.-ecc*ecc)^1.5); 
pot3 = 1.5/((1.-ecc*ecc)^1.5);
%pot4 = 0.375*(1.- 2.5*ecc*ecc);
%pot5 = 0.375*(1. + 5.5*ecc*ecc);

% case 1, Earthlike
Mass=1.*MassEarth;
Radius=1.0*RadiusEarth;
Rcore = 0.5178*Radius;
rhom=4500.0;
rhof= (3.*Mass/4./pi - rhom*(Radius^3 -Rcore^3) ) / Rcore^3;

% case 2, 0.1*Earthlike
%rhom=3300;
%rhof=6000;
%Radius=0.1*RadiusEarth;
%Rcore = 0.5178*Radius;
%Mass = (4.*pi/3)*(rhom*(Radius^3-Rcore^3) + rhof*Rcore^3);

%------------

rhomean= (3.*Mass/4./pi ) / Radius^3;
Mcore=(4.*pi/3.)*rhof*Rcore^3;
Mm=(4.*pi/3.)*rhom*(Radius^3-Rcore^3);
  
Rfr5=(Rcore/Radius)^5;
Rfr3=(Rcore/Radius)^3;

% moments of inertia
Af = 8.*pi/15.*rhof*(Rcore^5);
A = Af + 8.*pi/15.*(rhom*(Radius^5-Rcore^5));
Am = A - Af;

      

% orbital frequency = rotation frequency
rorbit=0.1*oneau;
rk=rorbit/oneau;
mean_anomaly = sqrt(bigG*MassStar/(rorbit^3));
rot=mean_anomaly;

  
day=24.*3600.0;
rotearth=2.*pi/day;
   
% calculate dynamical ellipticites
  mm=rot^2*Radius^3/Mass/bigG;
  xcmb=Rcore/Radius;
%  [ellf,ell] = clairault(mm,rhom,rhof,xcmb);
  [ellf,ell] = clairault_homog(mm,rhom,rhof,xcmb);
  
  scale=pi/180.;
  
  kappa = 4e-4;
  
thetam=zeros(npts,nicond);
thetaf=zeros(npts,nicond);
epsm=zeros(npts,nicond);
pow=zeros(npts,1);
  
for k=1:npts  % loop over Poincaré number
    

  pow(k)=(-5.5+2*(k-1)/(npts-1));
  dsigma=10.^pow(k);
     
  % compute large angle solution with newtonSys
  
    xtol=1.e-12; % choose tolerance
    ftol=1.e-12; % choose tolerance
    maxit=50;    % maximum number of iterations
    xmin=-90.0*scale; % min angle
    xmax=90.0*scale;  % max angle
    x0=zeros(3,1);
    isol=0;
  % loop through a set of initial conditions
     for i=1:nicond
         
         x0(1,1)= 0.0;
         x0(2,1)=(xmin + (xmax-xmin)*rand(1));
         x0(3,1)=(xmin + (xmax-xmin)*rand(1)); % rand between [-90, 90]
         
       [x_root,suc] = newtonSys_mod('J_exo_epsm',x0,xtol,ftol,maxit,0);
    
    
       if (suc==0)
           x_root(:)=1.e3;
       end
       if (x_root(2)<xmin) || (x_root(2)>xmax)
           x_root(:)=1.e3;
       end
%       if  (x_root(3)<xmin)|| (x_root(3)>xmax)
%           x_root(:)=1.e3;
%       end
       if  (x_root(3)<-88.0*scale)|| (x_root(3)>88.0*scale)
           x_root(:)=1.e3;
       end
           
       thetam(k,i)=x_root(1)/scale;  % in degrees
       thetaf(k,i)=x_root(2)/scale; % in degrees
       epsm(k,i)=x_root(3)/scale; % in degrees
       
     end
     
     % solution we pick is the one which has lowest thetaf
     
       [dummy,imin]=min( abs(thetaf(k,:)));
       thetaf0(k)=thetaf(k,imin);  % in degrees
       thetam0(k)=thetam(k,imin);  % in degrees
       epsm0(k)=epsm(k,imin); % in degrees
     
       Qcmb(k)=kappa*4*pi^2/3*Af*(abs(thetaf0(k)*scale)*rot)^3; 
       
     % For initial condition x= 0
     
     %    x0(1,1)= 0.0;
     %    x0(2,1)= 0.0;
     %    x0(3,1)= 0.0;
     %    
     %  [x_root,suc] = newtonSys_mod('J_exo_epsm',x0,xtol,ftol,maxit,0);
    %
    %   thetam0(k)=x_root(1)/scale;  % in degrees
    %   thetaf0(k)=x_root(2)/scale; % in degrees
    %   epsm0(k)=x_root(3)/scale; % in degrees
     
  %end
end

figure(1)
hold off
plot(pow,epsm0,'r')
xlabel('log10 of Poicare number'); ylabel('epsilon_m')
hold on
for i=1:nicond
    plot(pow,epsm(:,i),'ro')
end
axis([pow(1) pow(npts) -90 90])
hold off

figure(2)
hold off
plot(pow,thetaf0,'r')
xlabel('log10 of Poicare number'); ylabel('theta_f')
hold on
for i=1:nicond
    plot(pow,thetaf(:,i),'ro')
end
axis([pow(1) pow(npts) -90 90])
hold off

figure(3)
hold off
plot(pow,thetam0,'r')
xlabel('log10 of Poicare number'); ylabel('theta_m')
hold on
for i=1:nicond
    plot(pow,thetam(:,i),'ro')
end
axis([pow(1) pow(npts) -2.e-3 2.e-3])
hold off

figure(4)
semilogy(pow,Qcmb,'r')
xlabel('log10 of Poicare number'); ylabel('log10 of Qcmb')
