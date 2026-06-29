%% Efficient Near-Field Reconstruction for Linear-Array MIMO Radar With Irregular 3-D Translational Motion
% Ding Changqing
clc;clear;close all

%% Enable vertical motion compensation?
use_height_mc = input('Enable vertical motion compensation? 1 = yes, 0 = no: ');
while ~ismember(use_height_mc, [0, 1])
    disp('Invalid input. Please enter 1 or 0.');
    use_height_mc = input('Use vertical motion compensation? Enter 1 for yes, 0 for no: ');
end

if use_height_mc == 1
    disp('Current setting: vertical motion compensation is enabled.');
else
    disp('Current setting: vertical motion compensation is disabled.');
end
%% load echo
load('sig_data.mat','sig_real_i16','sig_imag_i16','sig_scale','sig_size');
% Recover complex single signal
sig = complex(single(sig_real_i16), single(sig_imag_i16));
sig = sig * single(sig_scale / 32767);
% Ensure original size
sig = reshape(sig, sig_size);
clear sig_real_i16 sig_imag_i16 sig_scale sig_size
disp('Pre-M2S signal has been loaded and recovered as single complex data.');
[MU, MV, MZ] = size(sig);
%% load params
load('params.mat', 'P');
c = P.radar.c;
Br = P.radar.Br;
f0 = P.radar.f0;
fc = P.radar.fc;
Kr = P.radar.Kr;
Tr = P.radar.Tr;
Fr = P.radar.Fr;
Fax = P.radar.Fax;
Vax = P.radar.Vax;
fr = P.radar.fr;
fax = P.radar.fax;

utx = P.array.utx;
urx = P.array.urx;
T_motion = P.motion;

fr_new = P.stolt.fr_new;
fr_new_min = P.stolt.fr_new_min;
fr_new_max = P.stolt.fr_new_max;
Neg = P.stolt.Neg;

dz = P.grid.dz;
X = P.grid.X;
Y = P.grid.Y;
Z = P.grid.Z;
lamda_stolt = P.grid.lamda_stolt;
%% MIMO-to-SISO conversion
MZ_up=1600;
sig_up = cat(3,sig,zeros(MU,MV,(MZ_up-MZ))); 
S_1 = ifft(sig_up,[],3);
% --------------------------------------------------------------------
r_r0 = (0:MZ_up-1)*c/2/Br/MZ_up*MZ;
r_r0(1) = 1e-6;
r_r1 = r_r0+(((utx-urx).').^2+(2.635e-2)^2)/8./r_r0;
delta_r = r_r1 - r_r0;
delta_r = reshape(delta_r,[MU,1,MZ_up]);
S_1 = S_1 .*exp(1j*4*pi*delta_r*f0/c);
clear delta_r sig_up sig
% --------------------------------------------------------------------
S_2 = zeros(MU,MV,MZ);
for nn = 1:MU  % interp1
    % --------------------------------------------------------------------
    pn = (r_r1(nn,:)>=0);
    PN = find(pn~=0);
    tem = interp1(r_r0,squeeze(S_1(nn,:,:)).',r_r1(nn,:).','linear',0); 
    S_2(nn,:,PN) = tem(PN,:).';
    % --------------------------------------------------------------------
end
sigM2S = fft(S_2,[],3);
sigM2S = sigM2S(:,:,1:MZ);
clear sig S_2 S_1 r_r0 r_r1
%% w1 azimuth density-compensation
nchose = 1;
if nchose == 1
    [ysort,I] = sort(T_motion(2,:),'ascend');
    w1 = gradient(ysort);
    w1(I) = w1;
    w1 = reshape(w1,[1,MV,1]);
    sigM2S = sigM2S.*w1;
end
%% FFT on x
NFFT_a = 2*MU;
sigM2S_MC_padded = zeros(NFFT_a,MV,MZ);
s1 = floor((NFFT_a - MU)/2);
sigM2S_MC_padded(s1+1:s1+MU,:,:) = sigM2S;
S_rd = fftshift(fft(sigM2S_MC_padded,[],1),1);
clear sigM2S_MC_padded sigM2S
%% Stolt interpolation
S_rd(Neg)=0; %only positive Ky is acceptable
S_3 = zeros(size(S_rd));
for nn = 1:NFFT_a  % interp1
    % --------------------------------------------------------------------
    pn = (abs(fr_new(nn,:))<MZ/2/Fr*Kr);
    PN = find(pn~=0);
    tem = interp1(fr,squeeze(S_rd(nn,:,:)).',fr_new(nn,:).','spline',0); 
    S_3(nn,:,PN) = tem(PN,:).';
    % --------------------------------------------------------------------
end
clear S_rd fr_new tem Neg
%% vertical motion compensation
fax_mtx = reshape(fax,[NFFT_a,1,1]);
H_refX = exp((1j*2*pi/Vax*fax_mtx)*reshape(T_motion(1,:),[1,MV,1]));
if use_height_mc == 1
    S_3 = S_3.*H_refX;
    disp('Vertical motion compensation has been completed.');
else
    disp('Vertical motion compensation is skipped.');
end   
clear fax_mtx fax H_refX
%% IFFT on r
S_4_padded = zeros(NFFT_a,MV,MZ_up);
s3 = floor((MZ_up - MZ)/2);
S_4_padded(:,:,s3+1:s3+MZ) = S_3;
Srd5 = ifft(ifftshift(S_4_padded,3),[],3);
clear S_3 S_4_padded
%% IFFT on x
Sac_t = ifft(ifftshift(Srd5,1),[],1);
s21 = floor((NFFT_a - MU)/2);
Sac_t = Sac_t(s21+1:s21+MU,:,:);
Sac_t(:,:,MZ_up) = 0;
clear Srd5
%% BP accumulation
[Ygrid,Zgrid]=meshgrid(Y,Z);
Ygrid = Ygrid.';Zgrid = Zgrid.';
[Ny,Nz] = size(Ygrid);
f_back=zeros(Ny,Nz,MU);
h=waitbar(0,'waiting');
for ii = 1:MV 
	R_ijk = sqrt((Zgrid+T_motion(3,ii)).^2+(Ygrid+T_motion(2,ii)).^2);
    R_ijk = R_ijk.*(abs(Ygrid+T_motion(2,ii))./abs(Zgrid+T_motion(3,ii))<=0.4570);
	t_ijk = R_ijk/dz+1;                     
    it_ijk=(t_ijk>1&t_ijk<=MZ_up);
 	t_ijk=t_ijk.*it_ijk+(MZ_up)*(1-it_ijk); 
    t_ijk1=ceil(t_ijk); 
    t_ijk2=floor(t_ijk);
    for mu = 1:MU
	    sig_rdta=Sac_t(mu,ii,:);
        sig_rdtachazhi=((t_ijk-t_ijk2)./(t_ijk1-t_ijk2+1e-9).*(sig_rdta(t_ijk1)-sig_rdta(t_ijk2))+sig_rdta(t_ijk2));
        f_back(:,:,mu)=f_back(:,:,mu)+sig_rdtachazhi.*exp(1j*4*pi*R_ijk/lamda_stolt); 
    end
 waitbar(ii/MV);
end
close(h);
clear Sac_t Ygrid Zgrid sig_rdta sig_rdtachazhi t_ijk t_ijk1 t_ijk2
f_back = flip(f_back,3); 
%% 2D maximum-intensity projection
f_2D = squeeze(max(abs(f_back),[],2));
[Ypic,Xpic]=meshgrid(Y,X);

figure;
    pcolor(Ypic',Xpic', 20*log10(rescale(abs((f_2D))))); 
    shading interp;
    colorbar;colormap(jet);
    xlabel('Y/m');
    ylabel('X/m');
    axis image
    clim([-40 0]);

%% figure: Measured 3-D translation vectors
figure('Color', 'w', 'Position', [250, 250, 560, 380]);
    frames = 1:MV; 
    h_x = plot(frames, T_motion(1,:), 'Color', [0.0, 0.45, 0.74], 'LineStyle', '-',  'LineWidth', 2); hold on;   
    h_y = plot(frames, T_motion(2,:), 'Color', [0.85, 0.33, 0.1], 'LineStyle', '-.', 'LineWidth', 2);            
    h_z = plot(frames, T_motion(3,:), 'Color', [0.47, 0.67, 0.19], 'LineStyle', '--', 'LineWidth', 2);            
    grid on; box on;
    xlabel('Frame Index (n_{\rho})', 'FontName', 'Times New Roman', 'FontSize', 12);
    ylabel('Displacement (m)', 'FontName', 'Times New Roman', 'FontSize', 12);
    legend([h_x, h_y, h_z], {'X component (\rho_x)', 'Y component (\rho_y)', 'Z component (\rho_z)'}, ...
           'FontName', 'Times New Roman', 'FontSize', 10, 'Location', 'best');
    ylim([-0.8 1.1])   
    ax4 = gca;
    set(ax4, 'FontName', 'Times New Roman', 'FontSize', 12, 'LineWidth', 1.1);
    ax4.GridColor = [0.85, 0.85, 0.85]; ax4.GridAlpha = 0.6;