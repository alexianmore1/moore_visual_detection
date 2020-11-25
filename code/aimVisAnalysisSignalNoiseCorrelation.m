

%%
% Align sample frequencies

frameSamples = (frameDelta:frameDelta:frameDelta*size(somaticF,2))*1000;
imTime  = ((frameSamples)/1000);
hertz = 1/frameDelta*1000;

% baseline data with df/f

blCutOffs = computeQuantileCutoffs(somaticF);
somaticF_BLs=slidingBaseline(somaticF,500,blCutOffs);
deltaFDS = (somaticF - somaticF_BLs)./somaticF_BLs;


% upsample data
tseries = timeseries(deltaFDS',imTime);
tseries_rs = resample(tseries,bData.sessionTime);
deltaF = tseries_rs.Data;


% pre/post are ms before/after stimulus onset
pre = 1000;
post = 3000;

% Considering whisker stimulation

[stimOnSamps,stimOnVect]=getStateSamps(bData.teensyStates,2,1);
[catchOnSamps,catchOnVect]=getStateSamps(bData.teensyStates,3,1);
stimSamps = [stimOnSamps, catchOnSamps];
stimSamps = sort(stimSamps);
whiskerSample = stimSamps; 
[rewardOnSamps, rewardOnVect] = getStateSamps(bData.teensyStates,4,1);
[licks, licksVect] = getStateSamps(bData.thresholdedLicks, 1, 1);
lickWindow = 1500;
zStackCellMat=[];

% Get the event triggered average for fluorescence
for cellID = 1:size(somaticF,1)
    cell = deltaF(:,cellID); % this is the cell you want
  
    for stimid = 2:100
        whiskfluormat(stimid,:) = cell(whiskerSample(stimid) - pre: whiskerSample(stimid)+ post);
    end
    zStackCellMat(:,:,cellID) = whiskfluormat;
end


% Make subplots for each cell
    figure
    for cellID = 1:size(somaticF, 1)
        subplot (3,4,cellID)
        %plot([-1*pre:post],(zstackcellmat(:,:,cellID)'))
        plot([-1*pre:post], mean((zStackCellMat(:,:,cellID)'),2)-mean(mean(zStackCellMat(:,1:900,cellID))))
        title(['Cell', num2str(cellID)])
        ylim([-0.05,0.1])
        xlim([-1000,3000])
    end

%%



% capture evoked scores for every trial of every cell
zStackCellMatBaselined = (zStackCellMat(:,:,:)  -mean(zStackCellMat(:,1:900,:),2));

preScore = zStackCellMatBaselined(:,1:1000,:);
postScore = zStackCellMatBaselined(:,1000:2000,:);

preScore2 = trapz(preScore,2);
preScore2 = squeeze(preScore2);

postScore2 = trapz(postScore,2);
postScore2 = squeeze(postScore2);

evoked = (postScore2-preScore2)';


% calculate correlation between cells
corrMatrix = corr(evoked');
corrVector = corrMatrix(:);
figure
title('Noise Correlation Across Cells')
imagesc(corrMatrix)
xlabel('Cell Number')
ylabel('Cell Number')
colorbar

% convenience downsample for now to 100
orientation = bData.orientation;
orientation = orientation(1:100);

contrast = bData.contrast;
contrast = contrast(1:100);

s{1} = orientation == 90;
s{2} = orientation == 0;
s{3} = orientation == 270;

sig{1} = mean(evoked(:,s1),2);
sig{2} = mean(evoked(:,s2),2);
sig{3} = mean(evoked(:,s3),2);

sig = [sig{1}' ; sig{2}' ; sig{3}' ];

sigcorr = corr(sig);

% 
for i = 1:3
   ncs{i} = corr(evoked(:,s{i})');
end

ncs = cell2mat(ncs);
ncs = reshape(ncs, 11,11,3);
ncs = mean(ncs, 3);

plot(get_upper(sigcorr), get_upper(ncs), 'o')

for i = 1:11
   sigmat(:,i) = sigmat(:,i)./norm(sigmat(:,i));
end

scatter()
