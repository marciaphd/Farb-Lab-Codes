function pfs=findPF(posFR,meanFR,contiguousPixels)

% Determine placefields from firing rate map.
%
% Arguments:
% posFR -- 2d array of firing rate map
% meanFR -- mean of posFR
% contiguousPixels -- minimum number contiguous pixels to identify
% placefield.
%
% Returns:
% pfs -- placefields and associated placefield characteristics
% pfs.posFRpf - the FR of all place fields (meeting the criterion)
% pfs.numberPFfound - the number of place fields
% pfs.posFRpfLargest - the largest place field
% pfs.pfLargestSize - the size of the largest place field
% pfs.pfLargestmeanFR - the mean FR of the largest place field
%
% To test this script, use findPF();
% which should return numberPFfound=5.


debugMode = 0;

if nargin == 0
    debugMode = 1;
    posFRs{1}=[2.1 3.1 4.1 5.1 6.1 1 2 .5 2 0;.5 .5 0 .8 7.1 0 3 .5 3 0;0 2 3 4 .2 0 4 0 4 5;0 7 6 5 0 0 5 0 6 0;...
        0 0 0 0 0 0 6 0 7 0;0 0 0 0 0 0 1 0 0 0;0 4 3 2 0 0 0 0 .6 0;0 5 0 5 0 9 8 7 7 8;...
        0 6 0 0 0 0 0 0 0 0;0 7 0 9 3 4 5 6 7 0];
    posFRs{2}=[12 3.1 4.1 5.1 4.1 1 2 .5 2 0;.5 .5 0 .8 7.1 0 3 .5 3 0;0 2 3 4 .2 0 4 0 4 5;0 7 6 5 0 0 5 0 6 0;...
        0 0 0 0 0 0 6 0 7 0;0 0 0 0 0 0 7 0 0 0;0 4 3 2 0 0 0 0 .6 0;0 5 0 5 0 9 8 7 7 8;...
        0 6 0 0 0 0 0 0 0 0;0 7 0 9 3 4 5 6 7 0];
    posFRs{3}=[0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0;...
        0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0;...
        0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0];
    
    contiguousPixels = 6;
    meanFR = 0.5;
    fprintf('...Testing findPFs.m with minimum bin firing rate of 1 Hz and at least %i contiguous pixels and firing rates:\n',...
        contiguousPixels)
    
    pfs=findPF(posFRs{1},meanFR,contiguousPixels);
    figure(200);clf
    subplot(2,2,1);imagesc(pfs.posFR);title('Original firerate map');axis xy, axis image; colorbar
    subplot(2,2,2);imagesc(pfs.allPFs);title('Place fields');axis xy, axis image; colorbar
    subplot(2,2,3);imagesc(pfs.largestMeanFRPF);title('Place field with largest meanFR');axis xy, axis image; colorbar
    subplot(2,2,4);imagesc(pfs.largestSizePF);title('Place field with largest number of pixels');axis xy, axis image; colorbar
    return
end

pfs = struct([]);  % initialize struct to hold results
pfs(1).posFR = posFR;

% Identify bins with firing rate greater than minimum
FF_FRcriterion = max(0.5, 2 * meanFR);
binaryImage = posFR > FF_FRcriterion;

% Identify each blob of contiguous firing pixels so can do calcs on it
% use 4 nearest neighbors, north, east, south, west to determine blobs
labeledImage = bwlabel(binaryImage, 4);  % 4 neighbors (diagonals not included)
% a=bwconncomp(posFR > FF_FRcriterion)  % newer method than bwlabel

% Get all properties for all blobs
blobMeasurements = regionprops(labeledImage, 'all');

% Determine which blobs have enough bins to be classified as a place field
idx2=find([blobMeasurements.Area] >= contiguousPixels);  % index into blobMeasurements struct

% Find largest blob ***actually only returns one if more than one are largest
posFRpf = zeros(size(posFR));
posFRpfLargest=zeros(size(posFR));
posFRpfLargestSize=zeros(size(posFR));


%% Keep only those blobs that qualify as a place field
if ~isempty(idx2)
    posFRpfmean = [];
    posFRpfpeak = [];
    posFRpfSize = [];
    for k=1:length(idx2)
        % copy bins (firing rates) from posFR to posFRpf that meet PF criteria
        % to make FR map
        posFRpf(blobMeasurements(idx2(k)).PixelIdxList)=posFR(blobMeasurements(idx2(k)).PixelIdxList);
        posFRpfmean(k)=mean(posFR(blobMeasurements(idx2(k)).PixelIdxList));  % mean FR of placefield blob
        posFRpfpeak(k)=max(posFR(blobMeasurements(idx2(k)).PixelIdxList));  % peak FR of placefield blob
        posFRpfSize(k) = length(blobMeasurements(idx2(k)).PixelIdxList);  % number pixels in blob
    end
    
    pfs.allPFs = posFRpf;
    
    % Of the valid placefield blobs, find the one with largest mean FR and
    % largest peak FR
    [~,pos]=max(posFRpfmean);  % pos is index into idx2 (=k) of blob
    %[pfLargestpeakFR,pos_pk]=max(posFRpfpeak);
    if debugMode;fprintf('Blob %i had highest mean FR of %5.3f\n',idx2(pos),pfLargestmeanFR);end
    if debugMode;fprintf('Blob %i had highest peak FR of %5.3f\n',idx2(pos_pk),pfLargestpeakFR);end
    
    % Make blob FR map with largest mean FR
    posFRpfLargest(blobMeasurements(idx2(pos)).PixelIdxList)=posFR(blobMeasurements(idx2(pos)).PixelIdxList);
    pfs.largestMeanFRPF = posFRpfLargest;
    
    % Make blob FR map with largest PF
    [~, pos_size] = max(posFRpfSize);
    posFRpfLargestSize(blobMeasurements(idx2(pos_size)).PixelIdxList)=posFR(blobMeasurements(idx2(pos_size)).PixelIdxList);
    pfs.largestSizePF = posFRpfLargestSize;
    pfs.largestSizePFinBins = length(find(posFRpfLargestSize));
    if debugMode; fprintf('Blob %i had largest PF size of %i pixels\n',idx2(pos_size),maxNumPixels);end
    
else
    pfs.allPFs = zeros(size(posFR));
    pfs.largestMeanFRPF = zeros(size(posFR));
    pfs.largestSizePF = zeros(size(posFR));
    pfs.largestSizePFinBins = 0;
end

% Display stats of PF blobs
if debugMode
    fprintf('\nSummary of PF blobs\n')
    fprintf('k     idx2  # pixels    meanFR   peakFR\n')
    if ~isempty(idx2)
        for k=1:length(idx2)
            fprintf('%i\t%i\t%i\t%5.2f\t%5.2f',k,idx2(k),blobMeasurements(idx2(k)).Area,posFRpfmean(k),posFRpfpeak(k))
            fprintf('\n')
        end
    end
end
pfs.numberPFfound=length(idx2);


