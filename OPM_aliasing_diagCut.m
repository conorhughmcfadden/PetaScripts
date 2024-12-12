clear all;
%%
mip = @(x,dim) squeeze(max(rescale(x),[],dim))';
pwr = @(x,dim) squeeze(log10(max(abs(fftshift(x)),[],dim)))';
p_rl = @(x,dim) squeeze(log10(max(real(fftshift(x)),[],dim)))';
p_im = @(x,dim) squeeze(log10(max(imag(fftshift(x)),[],dim)))';
apply = @(mask,x) fftshift(mask) .* x;
% dataPath = '/archive/bioinformatics/Danuser_lab/Fiolka/MicroscopeDevelopment/OPM/ZoomSystem/SiliconeOil/100nmBeads/green/241105';

%%
% dataPath = '/endosome/archive/bioinformatics/Danuser_lab/Fiolka/MicroscopeDevelopment/omniOPM/PSFparade/Beads100nm/YG/240926';
% imPath = fullfile(dataPath, 'Cell1', '1_CH00_000000.tif');

% beads
% dataPath = '/archive/bioinformatics/Danuser_lab/Fiolka/MicroscopeDevelopment/omniOPM/Reto/Beads/YG100/241125/';
% imPath = fullfile(dataPath, 'Cell11', '1_CH00_000000.tif');

% big beads
% dataPath = '/archive/bioinformatics/Danuser_lab/Fiolka/MicroscopeDevelopment/omniOPM/Reto/Beads/YG100nm/241203';
% imPath = fullfile(dataPath, 'Cell1', '1_CH00_000000.tif');

% big beads
dataPath = '/archive/bioinformatics/Danuser_lab/Fiolka/MicroscopeDevelopment/omniOPM/Retp/Beads/100nmYG/241209';
imPath = fullfile(dataPath, 'Cell6', '1_CH00_000000.tif');

% mito
% dataPath = '/archive/bioinformatics/Danuser_lab/Fiolka/MicroscopeDevelopment/omniOPM/U2OS/TOMM20/241118';
% imPath = fullfile(dataPath, 'Cell14', '1_CH00_000000.tif');

% dataPath = '/endosome/archive/bioinformatics/Danuser_lab/Fiolka/LabMembers/Conor/OPMv3/Aliasing/100nm_YG';
% imPathFull = fullfile(dataPath, '2024-09-29', 'dz_200nm_001', 'CH00_000000.tiff');
% imPath = fullfile(dataPath, '2024-09-29', 'dz_1000nm_001', 'CH00_000000.tiff');

%% deskewing

% beads
ovFactor = 2;
dsFactor = 3;
xyPixelSize = 0.147;
dz = 0.210;
skewAngle = 45.0;

% mito
% dsFactor = 3;
% xyPixelSize = 0.147;
% dz = 0.190;
% skewAngle = 50.0;

%% load, deskew and downsample image

im = permute(readtiff(imPath), [2,1,3]);
im = im(:, :, 1:ovFactor:end);

im_full = im;
im_dsp = im(:, :, 1:dsFactor:end);
% im_ds = im;

fillVal = median(im_dsp(:));

% deskew
im_dsk = deskewFrame3D( ...
    im_dsp, ...
    skewAngle, ...
    dz * dsFactor * ovFactor, ...
    xyPixelSize, ...
    'reverse', true ...
    );

im_full = deskewFrame3D( ...
    im_full, ...
    skewAngle, ...
    dz * ovFactor, ...
    xyPixelSize, ...
    'reverse', true ...
    );

% median fill:
% im_dsk(im_dsk(:) == 0) = fillVal;

%%
x1 = find(squeeze(im_full(1,:,1)) > 0);
x1 = x1(1);

x2 = find(squeeze(im_full(1,:,end)) > 0);
x2 = x2(end);

im_full = im_full(:, (x1+1):(x2-1), :);
im_dsk = im_dsk(:, (x1+1):(x2-1), :);

figure(1); clf;
set(gcf, 'color', [1,1,1]);
imagesc(mip(im_full, 1));
% imagesc(rescale(squeeze(im_full(256, :, :))));

figure(2); clf;
set(gcf, 'color', [1,1,1]);
imagesc(mip(im_dsk, 1));
axis image;
% title('stack (deskewed, downsampled)');

%% "upsampling" downsampled stack: G

G = fftn(im_dsk);
G = repmat(G, [1, 1, dsFactor]);

% handle even downsample case;
if ~mod(dsFactor, 2)
    z_ds = uint16(size(G,3)/dsFactor) + mod(size(G,3),2);
    G = cat(3, G(:, :, (z_ds+1):end), G(:, :, 1:z_ds));
end

figure(3); clf;
set(gcf, 'color', [1,1,1]);

% subplot(1,2,1);
imagesc(pwr(G, 1));
colormap hot;
axis image;
% title("fft(stack) repeated: G");
pbaspect([1,1,1]);

%%
[sy, sx, sz] = size(G);

[x, y, z] = meshgrid(1:sx, 1:sy, 1:sz);
x = x - mean(x(:));
y = y - mean(y(:));
z = z - mean(z(:));

blurSize = 0.0;

th = (cosd(skewAngle)*sz + sind(skewAngle)*sx)/2 - blurSize/2;
mask = (z > -cosd(skewAngle).*(x + th/dsFactor));
% mask = mask | (z > 0);
mask = mask & mask & flip(flip(mask, 3), 2);

if blurSize > 0
    mask = imgaussfilt3(double(mask), blurSize/dsFactor);
end

G_mask = fftshift(mask) .* G;

figure(4); clf; set(gcf, 'color', [1,1,1]);
imagesc(pwr(G_mask, 1));
colormap hot;
axis image;

%%
g_recon = ifftn(G_mask);
g_recon = real(g_recon);

figure(5); clf; set(gcf, 'color', [1,1,1]);
imagesc(mip(g_recon, 1));
% imagesc(rescale(squeeze(g_recon(256, :, :))));
colormap parula;

%%
rz = size(im_full,3) / size(im_dsk,3);

S = [1 0 0 0
    0 1 0 0
    0 0 rz 0
    0 0 0 1];

im_interp = imwarp(im_dsk, affine3d(S), "cubic");

figure(6); clf;
imagesc(mip(im_interp, 1));

%%
outSize = [512,512,128];

im_full_rot = rotateFrame3D( ...
    im_full, ...
    skewAngle, ...
    1.0, ...
    'reverse', true, ...
    'Crop', true, ...
    'outSize', outSize ...
    );
im_full_rot = norm_u16(im_full_rot);

figure(7); clf;
imagesc(mip(im_full_rot, 1));

writetiff(im_full_rot, fullfile(dataPath, "im_full_rot.tif"));

%%
im_recon_rot = rotateFrame3D( ...
    g_recon, ...
    skewAngle, ...
    1.0, ...
    'reverse', true, ...
    'Crop', true, ...
    'outSize', outSize ...
    );
im_recon_rot = norm_u16(im_recon_rot);

figure(8); clf;
imagesc(mip(im_recon_rot, 1));

writetiff(im_recon_rot, fullfile(dataPath, "im_recon_rot.tif"));

%%
im_interp_rot = rotateFrame3D( ...
    im_interp, ...
    skewAngle, ...
    1.0, ...
    'reverse', true, ...
    'Crop', true, ...
    'outSize', outSize ...
    );
im_interp_rot = norm_u16(im_interp_rot);

figure(9); clf;
imagesc(mip(im_interp_rot, 1));

writetiff(im_interp_rot, fullfile(dataPath, "im_interp_rot.tif"));

%% functions

function [ out ] = norm_u16( in )
    out = double(in);
    out = out - median(out(:));
    out(out < 0) = 0;
    out = out ./ max(out(:));
    out = uint16(65535 .* out);
end