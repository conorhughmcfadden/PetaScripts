clear all;
%%
mip = @(x,dim) squeeze(max(rescale(x),[],dim))';
pwr = @(x,dim) squeeze(log10(max(abs(fftshift(x)),[],dim)))';
p_rl = @(x,dim) squeeze(log10(max(real(fftshift(x)),[],dim)))';
p_im = @(x,dim) squeeze(log10(max(imag(fftshift(x)),[],dim)))';
apply = @(mask,x) fftshift(mask) .* x;

%%
% dataPath = '/archive/bioinformatics/Danuser_lab/Fiolka/MicroscopeDevelopment/omniOPM/Calibration60X/mito/OMP/241211';
% dataPath = '/archive/bioinformatics/Danuser_lab/Fiolka/MicroscopeDevelopment/Aliasing_OPM/U2OS/241211';
dataPath = '/archive/bioinformatics/Danuser_lab/Fiolka/LabMembers/Conor/omniOPM/Reto/Aliasing_OPM/U2OS/241211';

cellStr = 'Cell19';

chanToken = 'CH00';

tiffList = dir(fullfile(dataPath, cellStr, '*.tif*'));

%% microscope params

% omniOPM
dsFactor = 4;
xyPixelSize = 0.147;
dz = 0.207;
skewAngle = 45.0;

%% reconstruction params

zCrop = 64;

%% set up save folder

condition = sprintf("recon_dz=%dnm", uint16((dz*dsFactor)*1000));
saveDir = fullfile(dataPath, cellStr, condition);

if ~exist(saveDir, "dir")
    mkdir(saveDir);
end

%% loop through timepoints

maxNumTimePoints = -1;

if maxNumTimePoints <= 0
    maxNumTimePoints = length(tiffList);
end

for t = 1:maxNumTimePoints

    timePoint = tiffList(t);

    fprintf("Processing t = %d : %s...\n", t, timePoint.name);

    display = (t == 1) || (t == maxNumTimePoints);

    if ~contains(timePoint.name, chanToken)
        continue;
    end

    imPath = fullfile(timePoint.folder, timePoint.name);

    %% load, deskew and downsample image
    
    im = permute(readtiff(imPath), [2,1,3]);
    
    % deskew
    im_dsk = deskewFrame3D( ...
        im, ...
        skewAngle, ...
        dz * dsFactor, ...
        xyPixelSize, ...
        'reverse', true ...
        );
    
    x1 = find(squeeze(im_dsk(1,:,1)) > 0);
    x1 = x1(1);
    
    x2 = find(squeeze(im_dsk(1,:,end)) > 0);
    x2 = x2(end);
    
    im_dsk = im_dsk(:, (x1+1):(x2-1), :);
    
    outSize = [size(im_dsk, [1,2]), zCrop];
    
    im_dsk_rot = rotateFrame3D( ...
        im_dsk, ...
        skewAngle, ...
        dsFactor, ...
        'reverse', true, ...
        'Crop', true, ...
        'outSize', outSize ...
        );
    im_dsk_rot = norm_u16(im_dsk_rot);

    writetiff(im_dsk_rot, fullfile(saveDir, timePoint.name));

    if display
        figure(t); clf;
        subplot(1,4,1);
        set(gcf, 'color', [1,1,1]);
        imagesc(mip(im_dsk_rot, 1));
        colormap(gca, "parula");
        axis image;
    end

    %% "upsampling" downsampled stack: G
    
    G = fftn(im_dsk);
    G = repmat(G, [1, 1, dsFactor]);
    
    % handle even downsample case;
    if ~mod(dsFactor, 2)
        z_ds = uint16(size(G,3)/dsFactor) + mod(size(G,3),2);
        G = cat(3, G(:, :, (z_ds+1):end), G(:, :, 1:z_ds));
    end
    
    if display
        subplot(1,4,2);
        imagesc(pwr(G, 1));
        colormap(gca, "hot");
        axis image;
    end

    %% create mask and apply to G
    [sy, sx, sz] = size(G);
    
    [x, y, z] = meshgrid(1:sx, 1:sy, 1:sz);
    x = x - mean(x(:));
    y = y - mean(y(:));
    z = z - mean(z(:));
    
    blurSize = 0.0;
    
    th = (cosd(skewAngle)*sz + sind(skewAngle)*sx)/2 - blurSize/2;
    mask = (z > -cosd(skewAngle).*(x + th/dsFactor));
    mask = mask & mask & flip(flip(mask, 3), 2);
    
    if blurSize > 0
        mask = imgaussfilt3(double(mask), blurSize/dsFactor);
    end
    
    G_mask = fftshift(mask) .* G;
     
    if display
        subplot(1,4,3);
        imagesc(pwr(G_mask, 1));
        colormap(gca, "hot");
        axis image;
    end

    %% ifft on G
    g_recon = ifftn(G_mask);
    g_recon = real(g_recon);

    %% rotate recon
    im_recon_rot = rotateFrame3D( ...
        g_recon, ...
        skewAngle, ...
        1.0, ...
        'reverse', true, ...
        'Crop', true, ...
        'outSize', outSize ...
        );
    im_recon_rot = norm_u16(im_recon_rot);

    if display
        subplot(1,4,4);
        imagesc(mip(im_recon_rot, 1));
        colormap(gca, "parula");
        axis image;
    end

    %% save
    writetiff(im_recon_rot, fullfile(saveDir, ['recon_', timePoint.name]));
end

%% functions

function [ out ] = norm_u16( in )
    out = double(in);
    out = out - median(out(:));
    out(out < 0) = 0;
    out = out ./ max(out(:));
    out = uint16(65535 .* out);
end