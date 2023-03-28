function [PSF,PSFfine,ax,axfine]=MTF2PSF(MTF,profileorder,FOV)
% latest change: PSF output is kept complex now. (29/10-20)
Nx = numel(MTF);
if nargin < 3
    ax = -Nx/2:Nx/2-1;
else
    dx = FOV/Nx;
    ax = [-FOV/2+dx/2:dx:FOV/2-dx/2];
end

%figure,plot(abs(MTF)), hold on,
if nargin < 2
    profileorder = 'linear'; %default profile order is linear ; 
end
    
if strcmp(profileorder,'lowhigh') || strcmp(profileorder,'lh')
   if rem(Nx,2) 
    order = [Nx-1:-2:2, 1:2:Nx];
    MTF = MTF(order);
   else
    order = [Nx:-2:2, 1:2:Nx-1];
    MTF = MTF(order);  
   end
end
%plot(abs(MTF))
%real sampled PSF:
PSF = (sqrt(Nx)*fftshift(ifft(ifftshift(MTF)))); %IM(:,Nx/2+1);
%PSF1 = ifft2c(MTF); same as above

%figure,plot([-FOV/2:pix:FOV/2-1],PSF), xlabel('[mm]')
%figure,plot(ax,PSF), xlabel('pixels')
%xlim([-10 10])
%title('PSF')

%upsampled PSF: 
pad = 4096;

if nargin < 3
    axfine = linspace(-Nx/2,Nx/2-Nx/pad,pad);
else
    dxpad = FOV/pad;
    axfine = [-FOV/2+dxpad/2:dxpad:FOV/2-dxpad/2];
end
PSFfine = (sqrt(pad)*ifftshift(ifft(MTF,pad)));

%figure,plot(linspace(-FOV/2,FOV/2-FOV/pad,pad),PSFfine),xlabel('[mm]')
%figure,plot(axfine,PSFfine), xlabel('Pixels')
%xlim([-10 10])
%title('Upsampled PSF')
