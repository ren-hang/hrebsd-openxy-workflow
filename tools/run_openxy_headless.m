function run_openxy_headless(openxyCode, dataDir, stem, outName, opts)
% Headless OpenXY run mirroring MainGUI/RunButton_Callback, for local copies only.
%   run_openxy_headless(codeDir, dataDir, stem, outName)
%   dataDir contains <stem>.ang, <stem>.txt and Corrected_EBSPs\Platinum_r0c0.tif ...
% Settings follow Brady App. E / lab protocol. Pattern center from Pathfinder
% (X,Y,L) via Brady E.1-E.3. Logs go to the MATLAB diary (set by the caller).
if nargin < 5, opts = struct; end
if ~isfield(opts,'thermoPC'),  opts.thermoPC = [16.31 6.86 19.09]; end % protocol practice values
if ~isfield(opts,'mperpix'),   opts.mperpix  = 80.64; end             % Brady E: "Screen size (um/pixel)"
if ~isfield(opts,'camElev'),   opts.camElev  = 0; end                  % Brady E: Camera Elevation 0
if ~isfield(opts,'nproc'),     opts.nproc    = max(1,str2double(getenv('OPENXY_NPROC'))); end % protocol: 3
if isnan(opts.nproc),          opts.nproc    = 1; end
if ~isfield(opts,'firstImage'),opts.firstImage = 'Platinum_r0c0.tif'; end

% Environment overrides (used by the synthetic validation runs; defaults keep the practice-data behaviour)
if ~isfield(opts,'pcDirect') && ~isempty(getenv('OPENXY_PC')), opts.pcDirect = sscanf(getenv('OPENXY_PC'),'%f')'; end % "x* y* z*" in OpenXY convention
if ~isempty(getenv('OPENXY_FIRSTIMAGE')), opts.firstImage = getenv('OPENXY_FIRSTIMAGE'); end
if ~isempty(getenv('OPENXY_CAMELEV')),    opts.camElev = str2double(getenv('OPENXY_CAMELEV')); end
if ~isempty(getenv('OPENXY_MPERPIX')),    opts.mperpix = str2double(getenv('OPENXY_MPERPIX')); end
if ~isfield(opts,'gnd'),                  opts.gnd = ~strcmp(getenv('OPENXY_GND'),'0'); end

cd(openxyCode);
addpath(genpath(openxyCode));
addpath(fullfile(fileparts(mfilename('fullpath')),'headless_shims')); % GUI dialogs -> errors (must precede built-ins)
if ~exist('temp','dir'), mkdir('temp'); end   % MainGUI error handler saves temp/ErrorSettings.mat

Settings = GetHROIMDefaultSettings();
Settings.ScanFilePath  = fullfile(dataDir,[stem '.ang']);
Settings.FirstImagePath = fullfile(dataDir,'Corrected_EBSPs',opts.firstImage);
Settings.OutputPath    = fullfile(dataDir,'out',outName);
if ~exist(fullfile(dataDir,'out'),'dir'), mkdir(fullfile(dataDir,'out')); end
Settings.ScanType = 'Square';
Settings.Material = 'Scan File';
Settings.DoParallel = opts.nproc;
Settings.DisplayGUI = false;
Settings.DoShowPlot = false;
Settings.Exit = false;
% Microscope
Settings.AccelVoltage = 20;
Settings.SampleTilt = 70*pi/180;
Settings.SampleAzimuthal = 0;
Settings.CameraElevation = opts.camElev*pi/180;
Settings.CameraAzimuthal = 0;
Settings.mperpix = opts.mperpix;
% Advanced
Settings.DoStrain = 1;
Settings.HROIMMethod = 'Real';
Settings.StandardDeviation = 2;
Settings.MisoTol = 5;
Settings.GrainRefImageType = 'IQ > Fit > CI';
Settings.GrainMethod = 'Grain File';
Settings.CalcDerivatives = opts.gnd;
Settings.GNDMethod = 'Full';
Settings.NumSkipPts = 0;
Settings.IQCutoff = 0;
Settings.DoDDS = false;
% ROI
Settings.ROISizePercent = 25;
Settings.NumROIs = 48;
Settings.ROIStyle = 'Annular';
Settings.ROIFilter = [2 50 1 1];
Settings.ImageFilterType = 'standard';
Settings.ImageFilter = [9 90 0 0];
% preprocessing / ROI overrides for the sweeps (env vars; defaults above = Brady App. E)
if ~isempty(getenv('OPENXY_FILTER')),    Settings.ImageFilter = sscanf(getenv('OPENXY_FILTER'),'%f')'; end   % "lc uc 0 0"
if ~isempty(getenv('OPENXY_ROIFILTER')), Settings.ROIFilter   = sscanf(getenv('OPENXY_ROIFILTER'),'%f')'; end % "lc uc 1 1"
if ~isempty(getenv('OPENXY_ROIPCT')),    Settings.ROISizePercent = str2double(getenv('OPENXY_ROIPCT')); end
if ~isempty(getenv('OPENXY_NROI')),      Settings.NumROIs = str2double(getenv('OPENXY_NROI')); end

% --- Scan file (same calls as MainGUI SetScanFields)
Settings = ImportScanInfo(Settings, Settings.ScanFilePath);

% --- Image fields (same as MainGUI SetImageFields)
pats = patterns.makePatternProvider(Settings);
x = pats.imSize(1);
Settings.PixelSize = x;
Settings.ROISize = round(Settings.ROISizePercent*.01*Settings.PixelSize);
Settings.PhosphorSize = Settings.PixelSize*Settings.mperpix;
Settings.imsize = [x, pats.imSize(2)];
Settings.patterns = pats;
% The image filter actually used by the pattern-provider route is pats.filter (defaults 9 / 90 px), NOT Settings.ImageFilter
% (which only the legacy ReadEBSDImage route reads). Apply the requested filter here.
f = Settings.ImageFilter;
pats.filter.lowerRadius = f(1); pats.filter.upperRadius = f(2);
if numel(f) >= 4, pats.filter.lowerSmoothing = logical(f(3)); pats.filter.upperSmoothing = logical(f(4)); end
fprintf('pattern-provider filter: lower %g upper %g (doFilter=%d)\n', pats.filter.lowerRadius, pats.filter.upperRadius, pats.filter.doFilter);
Settings.patterns = pats;   % value object: re-assign after changing the filter
fprintf('PixelSize=%d  PhosphorSize=%g um  ScanLength=%d  Nx=%d Ny=%d\n', ...
    Settings.PixelSize, Settings.PhosphorSize, Settings.ScanLength, Settings.Nx, Settings.Ny);

% --- Pattern center: Brady E.1-E.3 applied to Pathfinder (X,Y,L)
if isfield(opts,'pcDirect')
    xs = opts.pcDirect(1); ys = opts.pcDirect(2); zs = opts.pcDirect(3);   % already OpenXY (x*,y*,z*)
else
    X = opts.thermoPC(1); Y = opts.thermoPC(2); L = opts.thermoPC(3);
    xs = 0.5 + 0.5*((X-15.8)/20);
    ys = 1 - Y/20;
    zs = L/20;
end
fprintf('PC (x*,y*,z*) = %.5f %.5f %.5f\n', xs, ys, zs);
Settings.ScanParams.xstar = xs; Settings.ScanParams.ystar = ys; Settings.ScanParams.zstar = zs;
Settings.XStar = xs*ones(Settings.ScanLength,1);
Settings.YStar = ys*ones(Settings.ScanLength,1);
Settings.ZStar = zs*ones(Settings.ScanLength,1);

Settings = HREBSDMain(Settings);
disp('run_openxy_headless: HREBSDMain returned');
end
