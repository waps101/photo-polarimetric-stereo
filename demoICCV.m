clear all
 
%% Parameters

bitdepth = 8;
n=1.5;
angles = (0:10:350).*(pi/180);
noise = 0;

%% Load Bunny

load bunnypol.mat

maskoriginal = mask;

% %% Synthetic lights
% maskoriginal = ones(201,201);
% phi2 = 0.1*(pi/180); 
% teta2 = 10.0*(pi/180);
% phi3 = 0.1*(pi/180);
% teta3 = 100.0*(pi/180);
% s= [sin(phi2)*cos(teta2); sin(phi2)*sin(teta2); cos(phi2)];
% t= [sin(phi3)*cos(teta3); sin(phi3)*sin(teta3); cos(phi3)];


%%

s = [0 1 3]';
s = s./norm(s);
t = [1 0 3]';
t = t./norm(t);

%%

Iun1 = render_diffuse(N,ones(size(maskoriginal)),s);
Iun2 = render_diffuse(N,ones(size(maskoriginal)),t);

mask = maskoriginal&Iun1>0&Iun2>0;

%% Generate synthetic data
im_diffuse = Iun1;
clear images
rho_d = ((n-1./n).^2.*sin(theta).^2)./(2+2.*n.^2-(n+1./n).^2.*sin(theta).^2+4.*cos(theta).*sqrt(n.^2-sin(theta).^2));
ImaxplusImin_d = 2.*im_diffuse;
ImaxminusImin_d = rho_d.*ImaxplusImin_d;
images = NaN(size(mask,1),size(mask,2),length(angles));
for i=1:length(angles)
    for row=1:size(mask,1)
        for col=1:size(mask,2)
            if mask(row,col)
                images(row,col,i)=im_diffuse(row,col)+(ImaxminusImin_d(row,col)./2).*cos(2.*angles(i)-2.*phi(row,col));
            end
        end
    end
end
disp([num2str(sum(mask(:))) ' pixels']);
clear ImaxplusImin_d ImaxminusImin_d ImaxplusImin_s ImaxminusImin_s rho_d rho_s phi_s
%% Corrupt synthetic images

% Add Gaussian image noise
images = images+randn(size(images)).*noise;
% Saturate
%images(images>1)=1;
%images(images<0)=0;
% Quantize
%images = round(images.*(2^bitdepth-1))./(2^bitdepth-1);

%% Estimate polarisation image from noisy images
[rho_est,phi_est,Iun1_est]=AMPolarisationImage(reshape(images,[size(images,1) size(images,2) 1 size(images,3)]),angles,mask);
theta_est = rho_diffuse(rho_est,1.5);

%% Estimate lighting
Iun2_est=Iun2+randn(size(mask)).*noise;
[ s_est,t_est,~ ] = est2lights( theta_est,phi_est,mask,Iun1_est,Iun2_est );

%% Reconstruct
T = [-1 0 0; 0 -1 0; 0 0 1];

options.usephase          = true;
options.usepolratio       = true;
options.useintensityratio = true;
options.usecolour         = false;
options.use2polratios     = true;

[ height1,mask ] = photoPolStereo( theta_est,Iun1_est,phi_est,s_est,mask,Iun2_est,t_est,options );
height1=height1-min(height1(mask));
vol1=sum(height1(mask));

[ height2,mask ] = photoPolStereo( theta_est,Iun1_est,phi_est,T*s_est,mask,Iun2_est,T*t_est,options );
height2=height2-min(height2(mask));
vol2=sum(height2(mask));
if vol1>vol2
    height=height1;
else
    height=height2;
end
%% Display

figure; surf(height','FaceColor', [1 1 1], 'EdgeColor', 'none', 'FaceLighting', 'phong'); axis equal; light
