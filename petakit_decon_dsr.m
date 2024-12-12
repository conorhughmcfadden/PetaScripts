%   Use your specific path to PetaKit5D here:
run("/archive/bioinformatics/Danuser_lab/Fiolka/LabMembers/Conor/MATLAB/PetaKit5D/setup.m");

% default dir when you first start script:
defaultDir = "/archive/bioinformatics/Danuser_lab/Fiolka/MicroscopeDevelopment";
% defaultDir = ".";

%% Microscope parameters:

xyPixelSize = 0.147;
dz = 0.240;
dzPSF = 0.212;
skewAngle = 45.0;
reverse = true; % {OmniOPM: true, OPMv2: false}

%% Set up data paths:
clear dir dirs;

if ~exist("rt", "var")
    rt = defaultDir;
end

channelPatterns = {};
dataPaths = {};
psfFullpaths = {};

while true
    [fName, fDir] = uigetfile(fullfile(rt, "*.tif*"), "Choose any stack in folder... (CANCEL to stop)");
    if ~fName
        break;
    end
    
    rt = fDir;

    if ~any(strcmp(dataPaths, rt))
        dataPaths{end+1} = rt;
    else
        fprintf("Folder %s \nhas already been added to dataPaths...\nNavigate to a new data folder to add, or CANCEL if done.\n\n", rt);
        continue;
    end

    fileList = dir(fullfile(rt, "*.tif*"));
    
    for i = 1:length(fileList)
        ch = strsplit(fileList(i).name, '_');
        ch = ch{end-1};
        
        if ~any(strcmp(channelPatterns, ch))
            channelPatterns{end+1} = ch;
        end
    end
end

if isempty(dataPaths)
    return;
end

%% Set up PSF paths:

[psfName, psfDir] = uigetfile(fullfile(rt, "*.tif*"), "Choose PSF file...");
psfPath = fullfile(psfDir, psfName);

% use same psf for each channel (for now):
[psfFullpaths{1:length(channelPatterns)}] = deal(psfPath);

doDecon = false;
if psfName
    doDecon = true;
end

%% -------- Deconvolution parameters --------

showOTFMasks = true;           % 'true' to show OTF masks

numberOfIterations = 2;         % 2-3 if 'omw', 20-30 otherwise
method = 'omw';                 % 'original'|'simplified'|'omw'
weinerAlpha = 0.005;            
OTFCumThresh = 0.6;             % OTF percentile threshold (< 1.0)
hannWindowBounds = [0.8, 1];

cpusPerTask = 72;
GPUJob = true;

if showOTFMasks && doDecon
    for i = 1:length(psfFullpaths)
        XR_visualize_OTF_mask_segmentation(psfFullpaths{i}, OTFCumThresh, false);
    end
end

%% ------------- Deconvolution --------------

deconResultDirName = '';

if doDecon
    deconResultDirName = 'decon';

    XR_decon_data_wrapper( ...
        dataPaths, ...
        'resultDirName', deconResultDirName, ...
        'overwrite', true, ...
        'channelPatterns', channelPatterns, ...
        'skewAngle', skewAngle, ...
        'dz', dz, ...
        'xyPixelSize', xyPixelSize, ...
        'background', 100, ...
        'dzPSF', dzPSF, ...
        'psfFullPaths', psfFullpaths, ...
        'deconIter', numberOfIterations, ...
        'RLMethod', method, ...
        'wienerAlpha', weinerAlpha, ...
        'OTFCumThresh', OTFCumThresh, ...
        'hannWinBounds', hannWindowBounds, ...
        'skewed', true, ...
        'Save16bit', true, ...
        'GPUJob', GPUJob, ...
        'cpusPerTask', cpusPerTask, ...
        'parseCluster', false, ...
        'largeFile', true ...
        );


    if GPUJob && gpuDeviceCount('available') > 0
        reset(gpuDevice);
    end
end

%% ------------- Deskew-Rotation -------------

dsrDataPaths = {};
for i = 1:length(dataPaths)
    dsrDataPaths{i} = fullfile(dataPaths{i}, deconResultDirName);
end

dsrResultDirName = 'dsr';

XR_deskew_rotate_data_wrapper( ...
    dsrDataPaths, ...
    'DSRDirName', dsrResultDirName, ...
    'overwrite', true, ...
    'channelPatterns', channelPatterns, ...
    'dz', dz, ...
    'xyPixelSize', xyPixelSize, ...
    'skewAngle', skewAngle, ...
    'reverse', reverse, ...
    'inputAxisOrder', 'xyz', ...
    'outputAxisOrder', 'xyz', ...
    'DSRCombined', true, ...
    'cpusPerTask', cpusPerTask, ...
    'parseCluster', false, ...
    'crop', true ...
    );