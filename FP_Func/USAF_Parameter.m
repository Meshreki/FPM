%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All main.m, AtlerMin.m, GDUpdate_Multiplication_rank1.m, USAF_Parameter.m, 
% hela_Parameter.m files are modified and extended to:
%
% 1- allowing a FULL reconstruction of the input images for 
% any dimesnions (i.e. not only square crops of them!).
% 2- saving images (and figures) for multiple variables.
% 3- saving dirac peaks positions.
% 4- adding more descriptive comments (in addition to those added by Ivo Ihrke).
% 4- investigating the huge error of the Algorithm to be due to pixels w/ 
% very large (intensity) values in the input stacks (and not due to
% the Algorithm itself!).
%
% last modified on 27.05.2022
% by John Meshreki, john.meshreki@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last modofied on 4/22/2014 
% Lei Tian (lei_tian@alum.mit.edu)


F = @(x) fftshift(fft2(ifftshift(x)));
Ft = @(x) fftshift(ifft2(ifftshift(x)));
row = @(x) x(:).';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% wavelength of illumination, assume monochromatic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%lambda = 0.6292;
lambda = 0.519; %CSE nimrud

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% numerical aperture of the objective
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NA = 0.1;
NA = 0.3; %CSE nimrud

% maximum spatial frequency set by NA
um_m = NA/lambda;
% system resolution based on the NA
dx0 = 1/um_m/2;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% magnification of the system,
% need to calibrate with calibration slides
% on 4x objective, front port
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%mag = 8.1485;
mag = 10; %CSE nimrud

%dpix_c = 6.5; %6.5um pixel size on the sensor plane
dpix_c = 6.4; %6.5um pixel size on the sensor plane %CSE nimrud
% effective image pixel size on the object plane
dpix_m = dpix_c/mag; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% # of pixels at the output image patch
% each patch will assign a single k-vector, the image patch size cannot be
% too large to keep the single-k assumption holds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Np = 100;

% FoV in the object space in x-direction
FoV_idx = Np(1)*dpix_m;
% FoV in the object space in y-direction
FoV_idy = Np(2)*dpix_m;
% sampling size at Fourier plane set by the image size (FoV)
% sampling size at Fourier plane is always = 1/FoV
if mod(Np,2) == 1
    du = 1/dpix_m/(Np(1)-1);
    dv = 1/dpix_m/(Np(2)-1);
else
    du = 1/FoV_idx;
    dv = 1/FoV_idy;
end

% low-pass filter diameter set by the NA = bandwidth of a single measurment
% in index
% N_NA = round(2*um_m/du_m);
% generate cutoff window by NA
m_x = 1:Np(1);
m_y = 1:Np(2);
[mm,nn] = meshgrid(m_x-round((Np(1)+1)/2), m_y-round((Np(2)+1)/2));
mm = mm';
nn = nn';
ridx = sqrt(mm.^2+nn.^2);
um_idx = um_m/du;

% assume a circular pupil function, lpf due to finite NA
w_NA = double(ridx<um_idx);
% h = fspecial('gaussian',10,5);
% w_NA = imfilter(w_NA,h);

% support of OTF is 2x of ATF(NA)
%Ps_otf = double(ridx<2*um_idx);

phC = ones(Np(1), Np(2));
aberration = ones(Np(1), Np(2));
pupil = w_NA.*phC.*aberration;

clear m mm nn


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up image corrdinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% original image size: 2160x2560
% can calibrate the center of the illumination with respect to the image by
% looking at the data from the dark/bright field image transitions
%ncent = [1080,1280];
ncent = [188,282]; %CSE nimrud
% start pixel of the image patch
%nstart = [981,1181];
%nstart = [801,1001];
nstart = [1,1];
% center, start & end of the image patch
img_ncent = [nstart(1)-ncent(1)+Np(1)/2,  nstart(2)-ncent(2)+Np(2)/2];
img_center = [(nstart(1)-ncent(1)+Np(1)/2)*dpix_m, (nstart(2)-ncent(2)+Np(2)/2)*dpix_m];
%img_start = nstart*dpix_m; %unused variable
%img_end = (nstart+Np)*dpix_m; %unused variable


%% LED array geometries 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spacing between neighboring LEDs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ds_led = 4e3; %4mm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% distance from the LED to the object
% experientally determined by placing a grating object
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%z_led = 67.5e3;
z_led = 66.5e3; %CSE nimrud

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% diameter of # of LEDs used in the experiment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%dia_led = 19;
dia_led = 3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up LED coordinates
% h: horizontal, v: vertical
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lit_cenv = 13;
lit_cenh = 14;
vled = [0:31]-lit_cenv;
hled = [0:31]-lit_cenh;

[hhled,vvled] = meshgrid(hled,vled);
rrled = sqrt(hhled.^2+vvled.^2);
LitCoord = rrled<dia_led/2;
% total number of LEDs used in the experiment
Nled = sum(LitCoord(:));
% index of LEDs used in the experiment
Litidx = find(LitCoord);

% corresponding angles for each LEDs
dd = sqrt((-hhled*ds_led-img_center(1)).^2+(-vvled*ds_led-img_center(2)).^2+z_led.^2);
sin_thetav = (-hhled*ds_led-img_center(1))./dd;
sin_thetah = (-vvled*ds_led-img_center(2))./dd;


illumination_na = sqrt(sin_thetav.^2+sin_thetah.^2);

% corresponding spatial freq for each LEDs
%
vled = sin_thetav/lambda;
uled = sin_thetah/lambda;
% spatial freq index for each plane wave relative to the center
idx_u = round(uled/du);
idx_v = round(vled/dv);


illumination_na_used = illumination_na(LitCoord);

% number of brightfield image
NBF = sum(illumination_na_used<NA);


% maxium spatial frequency achievable based on the maximum illumination
% angle from the LED array and NA of the objective
um_p = max(illumination_na_used)/lambda+um_m;
% resolution achieved after freq post-processing
dx0_p = 1/um_p/2;

disp(['synthetic NA is ',num2str(um_p*lambda)]);

% assume the max spatial freq of the original object
% um_obj>um_p
% assume the # of pixels of the original object
N_obj = [round(2*um_p/du)*2, round(2*um_p/dv)*2] ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% need to enforce N_obj/Np = integer to ensure no FT artifacts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_obj = [ceil(N_obj(1)/Np(1))*Np(1), ceil(N_obj(2)/Np(2))*Np(2)];
% max spatial freq of the original object
um_obj = [du*N_obj(1)/2, dv*N_obj(2)/2];

% sampling size of the object (=pixel size of the test image)
dx_obj = [1/um_obj(1)/2, 1/um_obj(2)/2];

% end
[xp,yp] = meshgrid([-Np(1)/2:Np(1)/2-1]*dpix_m, [-Np(2)/2:Np(2)/2-1]*dpix_m);

%x0 = ([-N_obj(1)/2:N_obj(1)/2/2-1]*dx_obj(1), [-N_obj(2)/2:N_obj(2)/2/2-1]*dx_obj(2));
[xx0,yy0] = meshgrid([-N_obj(1)/2:N_obj(1)/2/2-1]*dx_obj(1), [-N_obj(2)/2:N_obj(2)/2/2-1]*dx_obj(2));

%% define propagation transfer function
[u,v] = meshgrid(-um_obj(1):du:um_obj(1)-du,-um_obj(2):dv:um_obj(2)-dv);

% Fresnel
% object defocus distance
z0=0;
H0 = exp(1i*2*pi/lambda*z0)*exp(-1i*pi*lambda*z0*(u.^2+v.^2));
% OR angular spectrum
% H0 = exp(1i*2*pi*sqrt((1/lambda^2-u.^2-v.^2).*double(sqrt(u.^2+v.^2)<1/lambda))*dz);

