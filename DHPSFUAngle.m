% DHPSFU - script to find double-dots and extract molecular coordinates from DHPSF microscopy data, pre-processed by extracting single localisations.

% Inputs - files:
% - Data for analysis in the form of a text file containing a list of single localisations, with coordinates in pixels. The recommended fitting software is GDSC SMLM PeakFit. To enhance batch processing, the code loops through all files with the correct extension in the supplied folder.
% - Calibration file, similarly formatted. The calibration must contain a single fiducial, imaged at equally-spaced z-coordinates in the desired range. E.g. it could have been imaged every 50 nm from -2 um to 2 um in z. NB! There must be one frame per z-coordinate, i.e. there must be a one-step movement between each frame.

% Inputs - parameters:
% - Pixel size
% - Precision cutoff, above which the Peak Fit localisations will be filtered out
% - The step size in z in the calibration sequence
% - The range of z or angle in which to fit the data. Used to remove parts of calibration going near 0 or pi
% - Filtering parameters: the initial values for filtering out pairs of localisations that are too far and too close, and the allowed deviation in the distance compared to the calibration

% Output will be written into the specified folder in the ViSP .3d format: x y z intensity frame, tab-separated.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

%% User-defined variables %%

folder = ''; %folder with data
%fT = 'Pos1_slice2_PF70nm'; %filename of the data file, without extension
ext = '.xls'; %extension of the data file
cols = [10 11 9 1 14]; % the column numbers, corresponding to the following info: [x y intensity frame precision]
skip = 9; %how many lines to skip when reading in data (i.e. num lines in header)
separator = ','; %separator in the data file, typically either '\t' or ','

calibFolder = '';
calib = 'calib.xls';
sepCalib = '\t'; %separator in the calibration file, typically either '\t' or ','

outputFolder = '';

pixelSize = 271; %nm
precision = 70; %cutoff in nm
calibStep = 40; %nm
rangeToFit = [600 3800]; %range for fitting (closed interval), nm for z and radians for angle
ZorAngleorFrame = 'Z'; %is rangeToFit wrt 'Angle' or 'Z' coordinate, or 'Frame'?

initDistFilt = [8 3]; %top and bottom limits for distance between dots in pixels; best left lenient if subsequent filters enabled
filterOut = true; %enable/disable all filters
filterCalibRange = true; % filter localisations out of the angular range of the calibration? If false, polynomial fit is extrapolated beyond the range. Recommended for use especially if calib range is not full
filterDist = true; %filter based on the expected distance between the dots
distDev = 0.2; %allowed relative deviation of the distance between the dots, compared to calibration
filterRat = true; %filter based on the ratio of intensities between the dots
intDev = 1; %allowed deviation of the ratio of the intensities, compared to calibration

%% Calibration %%

[dx dy dz dd dr aRange] = calibAngle([calibFolder calib], sepCalib, calibStep, 'rangeToFit', rangeToFit, 'ZorAngleorFrame', ZorAngleorFrame);

%% Looping %%
files = dir(fullfile(folder, ['*' ext]));
fTs = {files.name};

for i=1:numel(fTs)
fT = regexprep(fTs{i}, ext, '');

%% Read in data, initial precision filter %%

x1 = zeros(3000000,1);
x2 = x1;
y1 = x1;
y2 = x1;
aAll = x1;
zAll = x1;
iAll = x1;
fAll = x1;
dAll = x1;
rAll = x1;

try
    data = dlmread([folder fT ext],separator,skip,0);
catch
    disp(['Error reading file ' folder fT ext ', maybe empty. Continuing...']);
    continue
end
data(find(data(:,cols(5))>precision),:)=[]; %precision cutoff

%% For each frame, find all pairs of dots that are at an appropriate distance from each other to be lobes of DHPSF, determinae and record their xy coordinates, the angle, the distance between them, their intensity ratio and average intensity %%

count = 1;

poses = accumarray(data(:,cols(4)),(1:size(data,1))',[],@(x){data(x,:)}); %separate the data into chunks corresponding to each frame

for i = 1:max(data(:,cols(4))) %iteration through frames
    if(i/100==floor(i/100)) %display progress
        disp([num2str(i) '/' num2str(max(data(:,1)))]);
    end
    
    thisPos = poses{i};
    if(isempty(thisPos))
        continue
    end
    xs = thisPos(:,cols(1));
    is = thisPos(:,cols(3));
    fs = thisPos(:,cols(4));
    if isempty(xs)
        continue
    end
    ys = thisPos(:,cols(2));

    D = squareform(pdist([xs ys],'euclidean')); %find distances between all localisations in frame i
    D(D>initDistFilt(1))=0; %if localisations are too far or too close, discard this pair (i.e. set distance to 0)
    D(D<initDistFilt(2))=0;
    D=triu(D);

    [I,J] = ind2sub(size(D),find(D)); %indices of all non-zero distances

    x1s = xs(I);
    x2s = xs(J);
    y1s = ys(I);
    y2s = ys(J);
    dists = D(sub2ind(size(D),I,J));
    ratio = (max([is(I) is(J)],[],2)./min([is(I) is(J)],[],2)); %ratio always >=1
    intensity = mean([is(I) is(J)],2);
    angles = atan2(ys(I)-ys(J),xs(I)-xs(J));
    angles(angles<0)=angles(angles<0)+pi;

    ll = length(x1s);
    ranges = count:(count+ll-1);
    x1(ranges) = x1s;
    x2(ranges) = x2s;
    y1(ranges) = y1s;
    y2(ranges) = y2s;
    aAll(ranges) = angles;
    iAll(ranges) = intensity;
    fAll(ranges) = fs(1)*ones(ll,1);
    dAll(ranges) = dists;
    rAll(ranges) = ratio;
    count = count+ll;

    %PLOTTING
    % plot(data(:,11),data(:,12),'.')
    % hold on
    
%     plot(xs,ys,'r.')
%     plot([x1s(1) x2s(1)],[y1s(1) y2s(1)],'.')
% % plot(xs,ys,'ro')
% %     hold on
% %     plot(x1s-3*cos(angles),y1s-3*sin(angles),'x')
% %     plot3(mean([x1s x2s],2),mean([y1s y2s],2),polyval(dz,angles),'.')
% %     hold off
% %     
% % %     hold on
% % %     plot(x1s+3*cos(angles),y1s+3*sin(angles),'x')
% %     pause(2)
    
end

x1(x1==0)=[];
x2(x2==0)=[];
y1(y1==0)=[];
y2(y2==0)=[];
aAll(aAll==0)=[];
iAll(iAll==0)=[];
fAll(fAll==0)=[];
dAll(dAll==0)=[];
rAll(rAll==0)=[];
xAll = mean([x1 x2],2);
yAll = mean([y1 y2],2);

%% Using the calibration data, find the full coordinates for each molecule %%

zmin = polyval(dz, aRange);
zN = polyval(dz,aAll); %find z, using angle
xN = (xAll-polyval(dx,zN)+polyval(dx, zmin))*pixelSize; %correct x and y for DHPSF-induced displacement, convert to nm
yN = (yAll-polyval(dy,zN)+polyval(dy, zmin))*pixelSize;

%% Quality filtering %%

if(filterOut)
    marker = ones(length(fAll),1);
    if(filterCalibRange)
        marker(aAll<aRange(1))=-1001; %Angle below calibration range
        marker(aAll>aRange(2))=-1002; %Angle above calibration range
    end
    if(filterDist)
        diffDist = abs(dAll-polyval(dd,aAll))./dAll; %relative difference between the actual distance between the dots and that expected from calibration
        marker(diffDist>distDev)=-1003; %deviation in distance larger than allowed
    end
    if(filterRat)
        diffRat = abs(rAll-polyval(dr,aAll)); %difference between the ratio of the dots' intensities and the expected ratio from calibration
        marker(diffRat>intDev)=-1004; %ratio off by more than allowed
        marker(rAll>4)=-1005; %ratio too large
    end

    rCutOff = find(marker>0);
    xN = xN(rCutOff);
    yN = yN(rCutOff);
    zN = zN(rCutOff);
    iAll = iAll(rCutOff);
    fAll = fAll(rCutOff);
end

%Plot
% figure
% colormap jet
% 
% scatter(xN,yN,0.2,zN)
% axis equal

%% Write the result %%
if (length(xN) > 0)
    disp(['Locs = ' num2str(length(xN))])

    mkdir(outputFolder);
    dlmwrite([outputFolder fT '.3d'],[xN yN zN iAll fAll],'delimiter','\t','precision',7);
else
    disp(['Locs = ' num2str(length(xN)) ', output not written'])
end

end
