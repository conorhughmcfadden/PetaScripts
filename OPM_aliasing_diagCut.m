clear all;
%%
mip = @(x,dim) squeeze(max(rescale(x),[],dim))';
% mip = @(x,dim) squeeze(x(256,:,:))';
pwr = @(x,dim) squeeze(log10(max(abs(fftshift(x)),[],dim)))';
p_rl = @(x,dim) squeeze(log10(max(real(fftshift(x)),[],dim)))';
p_im = @(x,dim) squeeze(log10(max(imag(fftshift(x)),[],dim)))';
apply = @(mask,x) fftshift(mask) .* x;

%%
% big beads
% dataPath = '/archive/bioinformatics/Danuser_lab/Fiolka/MicroscopeDevelopment/omniOPM/Calibration60X/mito/OMP/241211';
% dataPath = '/archive/bioinformatics/Danuser_lab/Fiolka/MicroscopeDevelopment/Aliasing_OPM/U2OS/241211';
dataPath = '/archive/bioinformatics/Danuser_lab/Fiolka/Manuscripts/OPM-ALIAS/DataToShare';

imPath = fullfile(dataPath, 'highOPM', 'Cell1', '1_CH00_000000.tif');
% imFullPath = fullfile(dataPath, 'Cell20', '1_CH00_000000.tif');

%% deskewing

% omniOPM
dsFactor = 3;
xyPixelSize = 0.147;
dz = 0.210;
skewAngle = 45.0;

% mesoOPM
% dsFactor = 4;
% xyPixelSize = 1.15;
% dz = 1.60;
% skewAngle = 45.0;

%% load, deskew and downsample image

im = permute(readtiff(imPath), [2,1,3]);
% im_full = permute(readtiff(imFullPath), [2,1,3]);
im_full = im;

im_dsp = im(:, :, 1:dsFactor:end);
% im_dsp = im;

fillVal = median(im_dsp(:));

% deskew
im_dsk = deskewFrame3D( ...
    im_dsp, ...
    skewAngle, ...
    dz * dsFactor, ...
    xyPixelSize, ...
    'reverse', true ...
    );

im_full = deskewFrame3D( ...
    im_full, ...
    skewAngle, ...
    dz, ...
    xyPixelSize, ...
    'reverse', true ...
    );

% median fill:
% im_dsk(im_dsk(:) == 0) = fillVal;

x1 = find(squeeze(im_dsk(1,:,1)) > 0);
x1 = x1(1);

x2 = find(squeeze(im_dsk(1,:,end)) > 0);
x2 = x2(end);

im_dsk = im_dsk(:, (x1+1):(x2-1), :);
im_full = im_full(:, (x1+1):(x2-1), :);

figure(1); clf;
set(gcf, 'color', [1,1,1]);
imagesc(mip(im_full, 1));
axis image;
title('Fully sampled (deskewed)');

figure(2);
set(gcf, 'color', [1,1,1]);
imagesc(mip(im_dsk, 1));
% imagesc(squeeze(im_dsk(512,:,:)));
axis image;

%% "upsampling" downsampled stack: G

nz_ds = size(im_dsk,3);
G = fftn(im_dsk);
G_rep = repmat(G, [1, 1, dsFactor]);

% handle even downsample case;
if ~mod(dsFactor, 2)
    z_ds = uint16(size(G_rep,3)/dsFactor) + mod(size(G_rep,3),2);
    G_rep = cat(3, G_rep(:, :, (z_ds+1):end), G_rep(:, :, 1:z_ds));
end

figure(3); clf;
set(gcf, 'color', [1,1,1]);

t = tiledlayout(2, 3, 'TileSpacing', 'none', 'Padding', 'compact');

nexttile(t);
imagesc(mip(im_full, 1));
colormap(gca, 'parula');
axis image;

nexttile(t);
imagesc(mip(im_dsk, 1));
colormap(gca, 'parula');
axis image;
set(gca, 'YTick', []);

im_rep = ifftn(G_rep);
im_rep = real(im_rep);

nexttile(t);
imagesc(mip(im_rep, 1));
colormap(gca, 'parula');
axis image;
set(gca, 'YAxisLocation', 'right');

G_full = fftn(im_full);

nexttile(t);
imagesc(pwr(G_full, 1));
colormap(gca, 'hot');
axis image;

nexttile(t);
imagesc(pwr(G, 1));
colormap(gca, 'hot');
axis image;
set(gca, 'YTick', []);

nexttile(t);
imagesc(pwr(G_rep, 1));
colormap(gca, 'hot');
axis image;

set(gca, 'YAxisLocation', 'right');
[sy, sx, sz] = size(G_rep);

[x, y, z] = meshgrid(1:sx, 1:sy, 1:sz);
x = x - mean(x(:));
y = y - mean(y(:));
z = z - mean(z(:));

blurSize = 0.0;

alpha = skewAngle + 3.0;
th = (cosd(alpha)*sz + sind(alpha)*sx)/2 - blurSize/2;
mask = (z > -cosd(alpha).*(x + th/dsFactor));
% mask = mask | (z > 0);
mask = mask & mask & flip(flip(mask, 3), 2);

hold on;

notMask = mip(~mask,1);
tri1 = bwtraceboundary(notMask, [1, 1], "E");
tri2 = bwtraceboundary(notMask, [size(notMask,1), size(notMask,2)], "W");

fill(tri1(:,2), tri1(:,1), [1,1,1], 'FaceAlpha', 0.2, 'EdgeColor', 'w');
fill(tri2(:,2), tri2(:,1), [1,1,1], 'FaceAlpha', 0.2, 'EdgeColor', 'w');

for j = 1:(dsFactor - mod(dsFactor,2))
    yline((j - ~mod(dsFactor,2) * 0.5) * nz_ds, '--', 'color', 'cyan');
end

hold off;

if blurSize > 0
    mask = imgaussfilt3(double(mask), blurSize/dsFactor);
end

%%
G_mask = fftshift(mask) .* G_rep;

figure(4); clf; set(gcf, 'color', [1,1,1]);
imagesc(pwr(G_mask, 1));
colormap hot;
axis image;
title('G masked');

%%
g_recon = ifftn(G_mask);
g_recon = real(g_recon);

figure(5); clf; set(gcf, 'color', [1,1,1]);
imagesc(mip(g_recon, 1));
% imagesc(rescale(squeeze(g_recon(256, :, :))));
colormap parula;
axis image;
title('Reconstructed (deskewed)');

%%
rz = size(im_full,3) / size(im_dsk,3);

S = [1 0 0 0
    0 1 0 0
    0 0 rz 0
    0 0 0 1];

im_interp = imwarp(im_dsk, affine3d(S), "cubic");

figure(6); clf;
imagesc(mip(im_interp, 1));
axis image;
title('Interpolated (deskewed)');

%%
outSize = [16,256,64];

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
set(gcf, 'color', [1,1,1]);
subplot(3,1,1);
imagesc(mip(im_full_rot, 1));
set(gca, 'XTick', []);
ylabel('fully sampled');
axis image; 

% 
% writetiff(im_full_rot, fullfile(dataPath, "im_full_rot.tif"));
% 

im_recon_rot = rotateFrame3D( ...
    g_recon, ...
    skewAngle, ...
    1.0, ...
    'reverse', true, ...
    'Crop', true, ...
    'outSize', outSize ...
    );
im_recon_rot = norm_u16(im_recon_rot);

subplot(3,1,2);
imagesc(mip(im_recon_rot, 1));
set(gca, 'XTick', []);
ylabel('reconstructed');
axis image;

% 
% writetiff(im_recon_rot, fullfile(dataPath, "im_recon_rot.tif"));
% 

im_interp_rot = rotateFrame3D( ...
    im_interp, ...
    skewAngle, ...
    1.0, ...
    'reverse', true, ...
    'Crop', true, ...
    'outSize', outSize ...
    );
im_interp_rot = norm_u16(im_interp_rot);

subplot(3,1,3);
imagesc(mip(im_interp_rot, 1));
ylabel('interpolated');
axis image;
% 
% writetiff(im_interp_rot, fullfile(dataPath, "im_interp_rot.tif"));

%% functions

function [ out ] = norm_u16( in )
    out = double(in);
    out = out - median(out(:));
    out(out < 0) = 0;
    out = out ./ max(out(:));
    out = uint16(65535 .* out);
end