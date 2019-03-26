%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Alternative Optimization method for polarised images.
%   Test only with one light sources with multi-channels
%   If you have multi-light source test with
%   AMPolarisationImageMulti.m
%   
%   Input:
%   @ polImages: poloarised images with size [row,colum,nChannel,nImages]
%   @ polAng: polarised angles, the size should consistant with nImages.
%   @ mask: mask for all channels and images.
%
%   Output:
%   @ rho_est: Degree of polarisation.
%   @ phi_est: Phase angle.
%   @ Iun_est: Unpolarised intensity.
%
%   DIZHONG ZHU 2017/3/2. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [rho_est,phi_est,Iun_est]=AMPolarisationImage(polImages,polAng,mask)

[rows,cols,nChannel,nImages]=size(polImages);

noofPixels=sum(mask(:));

Ipol=zeros(noofPixels,nImages,nChannel);
for i=1:nChannel
    for j=1:nImages
        I=polImages(:,:,i,j);
        Ipol(:,j,i)=I(mask);
    end
end

%% 1. pick up one channel for initialization
I=Ipol(:,:,1)';
A = [ones(nImages,1) cos(2.*polAng') sin(2.*polAng')];
x=A\I;
x=x';
Imax = x(:,1)+sqrt(x(:,2).^2+x(:,3).^2);
Imin = x(:,1)-sqrt(x(:,2).^2+x(:,3).^2);
Iun = (Imin+Imax)./2;
rho = (Imax-Imin)./(Imax+Imin+eps);
phi = 0.5*atan2(x(:,3),x(:,2));
phi = mod(phi,pi);
rho=rho';
phi=phi';

Iun_o=Iun;
%% 2. Do Alternative optimization

% A is the matrix with unpolarised intensity and polarised angles.
% X is the column contain DoP and phi angles.
% d is the observed polarised intensity images with color channels.

% Build  matrix for AX=D, X contains information of DoP and phi angle
Af_row_idx=repmat([1:1:noofPixels*nChannel*nImages],[2 1]);
Af_row_idx=Af_row_idx(:);
Af_col_idx=repmat([1:1:2*noofPixels],[1 nChannel*nImages])';
Af_value=zeros(noofPixels*nChannel*nImages*2,1);
x=[rho.*cos(2*phi) ;rho.*sin(2*phi)];
X=x(:);
D=zeros(noofPixels*nImages*nChannel,1);

% Build matrix for TI=B
T_row_index=[1:1:noofPixels*nImages]';
T_col_index=[1:1:noofPixels]';
T_col_index=repmat(T_col_index,[nImages 1]);
T_value=zeros(noofPixels*nImages,1);
B=reshape(Ipol,noofPixels*nImages,nChannel);

% Permute and stack the observed images
Ipol=permute(Ipol,[1 3 2]);

nIter=0;
maxIter=20;
converged=false;

tolerance=1e-15;

while ~converged&&nIter~=maxIter
    % Compute matrix A with initialized value DoP and phase angles
    
    %     % test only show images
%     phi2 = zeros(rows,cols);
%     phi2(mask) = phi;
%     rho2 = zeros(rows,cols);
%     rho2(mask) = rho;
%     figure();imagesc(phi2);axis equal;
    %     figure();imagesc(rho2);axis equal;
    
    % a. Given DoP and phase anlge estimate Unpolarised Intensity
    for j=1:nImages
        t=1+rho.*cos(2*polAng(j)-2*phi);
        lb=noofPixels*(j-1)+1;ub=noofPixels*j;
        T_value(lb:ub)=t;
    end
    T=sparse(T_row_index,T_col_index,T_value,noofPixels*nImages,noofPixels);
    Iun=T\B;
    %     maskset=repmat(mask,[1 1 nChannel]);
    %     II=zeros(rows,cols,nChannel);
    %     II(maskset)=Iun;
    %     figure();imshow(II);
    % Prepare construct matrix A
    
    % b. Given unpolarised intensity, estimate DoP and phase angles.
    tIun=repmat(Iun(:),[1 2]);
    for j=1:nImages % same as length of polarization angles
        A_value=tIun.*repmat([cos(2*polAng(j)) sin(2*polAng(j))],[nChannel*noofPixels 1]);
        A_value=A_value';
        lb=(j-1)*noofPixels*nChannel*2+1;ub=j*noofPixels*nChannel*2;
        Af_value(lb:ub)=A_value(:);
        
%         lb=(j-1)*nChannel+1;ub=j*nChannel;
        Itemp=Ipol(:,:,j)-Iun;
%         Itemp=Ipol(:,lb:ub)-Iun;
        lb=(j-1)*noofPixels*nChannel+1;ub=j*noofPixels*nChannel;
        D(lb:ub,:)=Itemp(:);
    end
    
    A=sparse(Af_row_idx,Af_col_idx,Af_value,noofPixels*nChannel*nImages,2*noofPixels);
    Xnew=A\D;
   residual=norm(Xnew-X);
    disp(['Residual = ' num2str(residual)]);
    if residual<tolerance
        converged=true;
    end
    
    % retrive the DoP and phase angle from optimized result;
    x=zeros(noofPixels,2);
    x(:,1)=Xnew(1:2:2*noofPixels);
    x(:,2)=Xnew(2:2:2*noofPixels);
    
    X=Xnew;
    phi=0.5*atan2(x(:,2),x(:,1));
    phi=mod(phi,pi);
    rho=sqrt(x(:,1).^2+x(:,2).^2);
    phi=phi';
    rho=rho';
    nIter=nIter+1;
end

phi_est=mod(phi,pi);
rho_est=rho;

% mask=old_mask;

% We have a mask so need to reshape the estimated quantities to the
% masked pixels only
phi2 = zeros(rows,cols);
phi2(mask) = phi_est;
phi_est = phi2;
rho2 = zeros(rows,cols);
rho2(mask) = rho_est;
rho_est = rho2;
Iun2=zeros(rows,cols,nChannel);
maskset=repmat(mask,[1 1 nChannel]);
Iun2(maskset)=Iun;
Iun_est=Iun2;

end