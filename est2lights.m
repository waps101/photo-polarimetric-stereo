function [ s,t,residual ] = est2lights( theta,phi,mask,Iun1,Iun2 )
%EST2LIGHTS Estimate light source directions from pair of polarisation ims
%   Inputs:
%      theta     - zenith angle implied by degree of polarisation
%      phi       - phase angle
%      mask      - binary foreground mask
%      Iun1,Iun2 - unpolarised intensity images
%
%   Outputs:
%      s,t       - unit vector light source directions
%      residual  - final sum squared error
%
% William Smith
% University of York
% 2017

options = optimoptions('lsqnonlin');
options.Display = 'iter';

% One of the possible normals from the polarisation image
N1(:,1)=sin(phi(mask)).*sin(theta(mask));
N1(:,2)=cos(phi(mask)).*sin(theta(mask));
N1(:,3)=cos(theta(mask));

% Possible gradients from polarisation image
p1 = N1(:,1)./N1(:,3);
p2 = -N1(:,1)./N1(:,3);
q1 = N1(:,2)./N1(:,3);
q2 = -N1(:,2)./N1(:,3);

s_init = randn(3,1);
s_init(3) = abs(s_init(3));
s_init = s_init./norm(s_init);
[init(1),init(2),~]=cart2sph(s_init(1),s_init(2),s_init(3));
t_init = randn(3,1);
t_init(3) = abs(t_init(3));
t_init = t_init./norm(t_init);
[init(3),init(4),~]=cart2sph(t_init(1),t_init(2),t_init(3));

st_ang  = lsqnonlin(@(st_ang) objfun(st_ang,p1,q1,p2,q2,Iun1(mask),Iun2(mask)),init,[],[],options);

residual = sum(objfun(st_ang,p1,q1,p2,q2,Iun1(mask),Iun2(mask)).^2);

[s(1),s(2),s(3)]=sph2cart(st_ang(1),st_ang(2),1);
[t(1),t(2),t(3)]=sph2cart(st_ang(3),st_ang(4),1);
s=s';
t=t';

end

function residuals = objfun(st_ang,p1,q1,p2,q2,i1,i2)%

[s(1),s(2),s(3)]=sph2cart(st_ang(1),st_ang(2),1);
[t(1),t(2),t(3)]=sph2cart(st_ang(3),st_ang(4),1);

a = i2.*s(1)-i1.*t(1);
b = i2.*s(2)-i1.*t(2);
c = i2.*s(3)-i1.*t(3);

residuals=min(abs(a.*p1+b.*q1-c),abs(a.*p2+b.*q2-c));

end
