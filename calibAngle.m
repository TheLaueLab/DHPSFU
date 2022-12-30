% Function to calibrate the xyz coordinates to the pitch of the double helix, using the calibration data series.
% It takes the calibration Peak Fit result file as the first argument, the step in z between consecutive frames as its second obligatory argument. Optional arguments include the column numbers in the file, corresponding to [x y intensity frame]; and the range in which to fit the data, either wrt angle/pitch of the double helix or the z coordinate. The optional arguments have defaults. Default column numbers correspond to the current formatting of Peak Fit results.
% The outputs are polynomial functions, describing the dependency of the z coordinate, the distance between the dots and the ratio of their intensities on the angle, as well as the dependency of the xy coordinates on the position in z. The final output value is the angular range of this particular calibration.

function [dx dy dz dd dr aRange model] = calib( fileN,sepCalib,step,varargin)
params = struct('cols', [10 11 9 1], 'rangeToFit', [600 3800], 'ZorAngleorFrame', 'Z');
paramNames = fieldnames(params);

nArgs = length(varargin);
if round(nArgs/2)~=nArgs/2
   error('calibAngle needs propertyName/propertyValue pairs')
end

for pair = reshape(varargin,2,[])

   if any(strcmp(pair{1},paramNames))
      params.(pair{1}) = pair{2};
   else
      error('%s is not a recognized parameter name',pair{1})
   end
end

cols = params.cols;
rangeToFit = params.rangeToFit;
ZorAngleorFrame = params.ZorAngleorFrame;

beadID = [];
x1 = [];
x2 = [];
y1 = [];
y2 = [];
aAll = [];
zAll = [];
iAll = [];
fAll = [];
ratAll = [];
distAll = [];

%Iterate over multiple calibration files
for j = 1:numel(fileN)
    data = dlmread(fileN{j},sepCalib,9,0);
    %% For each frame, find the two lobes, extract their xy coordinates, the distance and the angle between them, their intensities and the ratio of the intensities %%
    for i = 1:max(data(:,cols(4))) %iteration through frames
        if(i/100==floor(i/100))
            i %display progress
        end
    
        poses = find(data(:,cols(4))==i);
        xs = data(poses,cols(1));
        is = data(poses,cols(3));
        fs = data(poses,cols(4));
        if isempty(xs)
            disp(['Warning: no localisations found in frame ', ...
                num2str(i), ', skipping'])
            continue
        end
        if length(xs) == 1
            disp(['Warning: only 1 localisation found in frame ', ...
                num2str(i), ', skipping'])
            continue
        end
        ys = data(poses,cols(2));
    
        D = squareform(pdist([xs ys],'euclidean')); %find distances     between all localisations in frame i
        D(D>8)=0; %if localisations are too far or too close,     discard this pair (i.e. set distance to 0)
        D(D<3)=0;
        D=triu(D);
    
        [I,J] = ind2sub(size(D),find(D)); %indices of all non-zero distances, should be just a single distance if there are no spurious localisations
        if length(I) > 1 %if >1 distance, skip this frame with a warning
            disp(['Warning: >2 localisations found in frame ', ...
                num2str(i), ', skipping'])
            continue
        end
        if length(I) == 0 %if distance too large/small, skip this frame with a warning
            disp(['Warning: localisations found in frame ', ...
                num2str(i), ' are too far/close, skipping'])
            continue
        end
    
        x1s = xs(I);
        x2s = xs(J);
        y1s = ys(I);
        y2s = ys(J);
        dists = D(I,J);
        intensity = mean([is(I) is(J)],2);
        ratio = max([is(I) is(J)])./min([is(I) is(J)]);
        angles = atan2(ys(I)-ys(J),xs(I)-xs(J));
        angles(angles<0)=angles(angles<0)+pi;
    
        x1 = [x1;x1s];
        x2 = [x2;x2s];
        y1 = [y1;y1s];
        y2 = [y2;y2s];
        aAll = [aAll;angles];
        iAll = [iAll;intensity];
        fAll = [fAll;fs(1)*ones(length(x1s),1)];
        ratAll = [ratAll;ratio];
        distAll = [distAll;dists];

        beadID = [beadID;j*ones(length(x1s),1)];
    end
end

xAll = mean([x1 x2],2);
yAll = mean([y1 y2],2);

%% Correcting the tilt %%
%For each bead at each frame, find the relative shift of angle from the mean angle of all beads
aRel = [];
for i = 1:max(fAll)
    poses = find(fAll==i);
    aRel(poses) = aAll(poses)-mean(aAll(poses));
end
%For each bead, find the average (over frames) relative shift wrt mean
aAvg = [];
xAvg = [];
yAvg = [];
for j = 1:max(beadID)
    poses = find(beadID==j);
    aAvg = [aAvg; mean(aRel(poses))];
    xAvg = [xAvg; mean(xAll(poses))];
    yAvg = [yAvg; mean(yAll(poses))];
end

%Fit linear model, return it later to apply to data
model = fitlm([xAvg yAvg], aAvg);
%Correct the angles
aAll = aAll - predict(model, [xAll yAll]);

%% Exclude the frames outside the specified range %%
if ZorAngleorFrame == 'Z'
    idx = find((fAll*step >= rangeToFit(1)) & (fAll*step <= rangeToFit(2)));
elseif ZorAngleorFrame == 'Angle'
    idx = find((aAll >= rangeToFit(1)) & (aAll <= rangeToFit(2)));
elseif ZorAngleorFrame == 'Frame'
    idx = find((fAll >= rangeToFit(1)) & (fAll <= rangeToFit(2)));
end
xAll = xAll(idx);
yAll = yAll(idx);
fAll = fAll(idx);
aAll = aAll(idx);
distAll = distAll(idx);
ratAll = ratAll(idx);
beadID = beadID(idx);


disp(aAll)
disp(length(aAll))

aRange = [min(aAll) max(aAll)];

%% Fit polynomial functions to... %%
%First, obtain normalised x and y coordinates of each bead, which has >95% observations
xAllNorm = [];
yAllNorm = [];
fAllNorm = [];
numObs = accumarray(beadID, beadID, [], @length);
for j = 1:max(beadID)
    poses = find(beadID==j);
    aAvg = [aAvg; mean(aRel(poses))];
    xAvg = [xAvg; mean(xAll(poses))];
    yAvg = [yAvg; mean(yAll(poses))];
    if length(poses) > 0.95*max(numObs)
        xAllNorm = [xAllNorm; xAll(poses)-mean(xAll(poses))];
        yAllNorm = [yAllNorm; yAll(poses)-mean(yAll(poses))];
        fAllNorm = [fAllNorm; fAll(poses)];
    end
end
dx = polyfit(fAllNorm*step,xAllNorm,5); % ... x vs z, defined as frame*z-step
dy = polyfit(fAllNorm*step,yAllNorm,5); % ... y vs z
dz = polyfit(aAll,fAll*step,5); % ... z vs angle

dd = polyfit(aAll,distAll,8); % ... distance vs angle
dr = polyfit(aAll,ratAll,8); % ... ratio of intensities vs angle

%% Plot %%
subplot(3,2,1)
plot(fAll*step,xAll,'o',fAll*step,polyval(dx,fAll*step),'-')
ylabel('xPos (px)')
xlabel('zPos (nm)')
subplot(3,2,2)
plot(fAll*step,yAll,'o',fAll*step,polyval(dy,fAll*step),'-')
ylabel('yPos (px)')
xlabel('zPos (nm)')
subplot(3,2,3)
plot(fAll*step,aAll,'o',polyval(dz,aAll),aAll,'-')
ylabel('angle (rad)')
xlabel('zPos (nm)')
subplot(3,2,4)
plot(fAll*step,distAll,'o',fAll*step,polyval(dd,aAll),'-')
ylabel('px-dist')
xlabel('zPos (nm)')
subplot(3,2,5)
plot(fAll*step,ratAll,'o',fAll*step,polyval(dr,aAll))
ylabel('ratio')
xlabel('zPos (nm)')

end
