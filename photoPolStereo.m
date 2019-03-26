function [ height,mask ] = photoPolStereo( theta,Iun,phi,l,mask,Iun2,l2,options )
%PHOTOPOLSTEREO Photo-polarimetric stereo
%   Inputs:
%      theta     - zenith angle implied by degree of polarisation
%      phi       - phase angle
%      l,l2      - light source directions in Iun and Iun2
%      mask      - binary foreground mask
%      Iun1,Iun2 - unpolarised intensity images

% Choose constraints (two must be true)
%usephase          = true;
%usepolratio       = true;
%useintensityratio = true;

%usecolour         = false; % Implies that Iun is rows by cols by 3
%use2polratios     = false;

usephase = options.usephase;
usepolratio = options.usepolratio;
useintensityratio = options.useintensityratio;
usecolour = options.usecolour;
use2polratios = options.use2polratios;

if usecolour
    IunR = Iun(:,:,1);
    IunG = Iun(:,:,2);
    IunB = Iun(:,:,3);
    if useintensityratio && usecolour
        Iun2R = Iun2(:,:,1);
        Iun2G = Iun2(:,:,2);
        Iun2B = Iun2(:,:,3);
    end
end

% Compute derivative matrices and update mask as appropriate
[ Dx,Dy,mask ] = gradMatrices( mask,'Backward' );

npix = sum(mask(:));

i = 1:npix;
j = 1:npix;

C = [];
d = [];

% Create linear system of equations for gradient terms
if usephase
    C = [C; sparse(i,j,-cos(phi(mask)),npix,npix) sparse(i,j,sin(phi(mask)),npix,npix)];
    d = [d; zeros(npix,1)];
end
if usepolratio
    if usecolour
        C = [C; sparse(i,j,ones(npix,1).*-l(1),npix,npix) sparse(i,j,ones(npix,1).*-l(2),npix,npix)];
        d = [d; IunR(mask)./cos(theta(mask))-l(3)];
        C = [C; sparse(i,j,ones(npix,1).*-l(1),npix,npix) sparse(i,j,ones(npix,1).*-l(2),npix,npix)];
        d = [d; IunG(mask)./cos(theta(mask))-l(3)];
        C = [C; sparse(i,j,ones(npix,1).*-l(1),npix,npix) sparse(i,j,ones(npix,1).*-l(2),npix,npix)];
        d = [d; IunB(mask)./cos(theta(mask))-l(3)];
    else
        C = [C; sparse(i,j,ones(npix,1).*-l(1),npix,npix) sparse(i,j,ones(npix,1).*-l(2),npix,npix)];
        if isfield(options,'albedo')
            d = [d; (Iun(mask)./options.albedo(mask))./cos(theta(mask))-l(3)];
        else
            d = [d; Iun(mask)./cos(theta(mask))-l(3)];
        end
    end
    
    if use2polratios
        % Iun2 and l2 must exist
        if usecolour
            C = [C; sparse(i,j,ones(npix,1).*-l2(1),npix,npix) sparse(i,j,ones(npix,1).*-l2(2),npix,npix)];
            d = [d; Iun2R(mask)./cos(theta(mask))-l2(3)];
            C = [C; sparse(i,j,ones(npix,1).*-l2(1),npix,npix) sparse(i,j,ones(npix,1).*-l2(2),npix,npix)];
            d = [d; Iun2G(mask)./cos(theta(mask))-l2(3)];
            C = [C; sparse(i,j,ones(npix,1).*-l2(1),npix,npix) sparse(i,j,ones(npix,1).*-l2(2),npix,npix)];
            d = [d; Iun2B(mask)./cos(theta(mask))-l2(3)];
        else
            C = [C; sparse(i,j,ones(npix,1).*-l2(1),npix,npix) sparse(i,j,ones(npix,1).*-l2(2),npix,npix)];
            if isfield(options,'albedo')
                d = [d; (Iun2(mask)./options.albedo(mask))./cos(theta(mask))-l2(3)];
            else
                d = [d; Iun2(mask)./cos(theta(mask))-l2(3)];
            end
        end
    end
end
if useintensityratio
    % Iun2 and l2 must exist
    if usecolour
        C = [C; sparse(i,j,Iun2R(mask).*l(1) - IunR(mask).*l2(1),npix,npix) sparse(i,j,Iun2R(mask).*l(2) - IunR(mask).*l2(2),npix,npix)];
        d = [d; Iun2R(mask).*l(3) - IunR(mask).*l2(3)];
        C = [C; sparse(i,j,Iun2G(mask).*l(1) - IunG(mask).*l2(1),npix,npix) sparse(i,j,Iun2G(mask).*l(2) - IunG(mask).*l2(2),npix,npix)];
        d = [d; Iun2G(mask).*l(3) - IunG(mask).*l2(3)];
        C = [C; sparse(i,j,Iun2B(mask).*l(1) - IunB(mask).*l2(1),npix,npix) sparse(i,j,Iun2B(mask).*l(2) - IunB(mask).*l2(2),npix,npix)];
        d = [d; Iun2B(mask).*l(3) - IunB(mask).*l2(3)];
    else
        C = [C; sparse(i,j,Iun2(mask).*l(1) - Iun(mask).*l2(1),npix,npix) sparse(i,j,Iun2(mask).*l(2) - Iun(mask).*l2(2),npix,npix)];
        d = [d; Iun2(mask).*l(3) - Iun(mask).*l2(3)];
    end
end
% Convexity constraint would go here since it uses surface gradients

% Post multiply by gradient matrices so that height finite differences are
% used
C = C * [Dx; Dy];

% Now add any equations that relate directly to depth

% Laplacian smoothness prior 
lambda=0;
smoothness=lambda*LaplacianMatrix(mask,mask);
C=[C;smoothness];
conRHS=zeros(npix,1);
d=[d;conRHS];

% Add row for resolving constant of integration
C(end+1,1)=1;
d = [d; 0];

% Solve linear system
z = C\d;
%disp(num2str(norm(C*z-d)));
% Copy height values back to height map
height = NaN(size(mask));
height(mask) = z;

end

