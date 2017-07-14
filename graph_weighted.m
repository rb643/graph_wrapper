function [wResult p] = graph_weighted(Matrix,varargin)
%
% For example:
% [wResult p] = graph_weighted('Matrix',adjacencymatrix);
% [wResult p] = graph_weighted('Matrix',adjacencymatrix,'regionLabels',regions,'PlotLocal',1,'PlotGlobal', 1,...
%            `                       'subjectmask',include, 'groups', groups, 'nRand', 100);
%
% ---------------------------- INPUT ----------------------------
%
% Matrix            -       3D connectivity matrix(node*node*subject)
%
% --------------------- OPTIONAL ARGUMENTS ----------------------
%
% regionLabels      -       A cell structure containing names of all nodes in your matrix (at this point its only used for plotting)
% filelist          -       list of filenames
% nRand             -       number of randomizations used in random networks for normalization (default = 10)
% groups            -       1D vector with grouplabels  (default = random)
% percentage        -       percentage at which to threshold the matrix
% subjectmask       -       Logical to set which subjects to include (default = all a.k.a. a row of ones)
% PlotLocal         -       Logical to set if we want to plot local metrics (default = 0)
% PlotGlobal        -       Logical to set if we want to plot global metrics (default = 0)
% PlotMatrices      -       Logical to set if we want to plot group average matrices (default = 0)
%
% ---------------------------- OUTPUT ----------------------------
% wResult           -       Graph metrics from weighted matrices based on NOS
% 
% p                 -       Also return all the input settings to check

%% check all the inputs and if they do not exist then revert to default settings
% input parsing settings
p = inputParser;
p.CaseSensitive = true;
p.Parameters;
p.Results;
p.KeepUnmatched = true;
% set the desired and optional input arguments
addRequired(p,'Matrix',@isnumeric);
addOptional(p,'regionLabels',[],@iscell);
addOptional(p,'filelist',[],@iscell);
addOptional(p,'nRand',100,@isnumeric);
addOptional(p,'PlotLocal',0,@isnumeric);
addOptional(p,'PlotGlobal',0,@isnumeric);
addOptional(p,'PlotMatrices',0,@isnumeric);
addOptional(p,'groups',[], @isnumeric);
addOptional(p,'subjectmask', [], @isnumeric);
addOptional(p,'percentage',0.1, @isnumeric);
% parse the input
parse(p,varargin{:});
% then set/get all the inputs out of this structure
nRand = p.Results.nRand; regionLabels = p.Results.regionLabels; PlotLocal = p.Results.PlotLocal; PlotGlobal = p.Results.PlotGlobal; groups = p.Results.groups; Adj = p.Results.Matrix;
PlotMatrices = p.Results.PlotMatrices; subjectmask = p.Results.subjectmask; percentage = p.Results.percentage; filenames = p.Results.filelist;

% set the larger defaults in case they are not specified
if isempty(groups); groups = round(rand(1,size(Adj,3))*3)'; end % generate some random numbers
if isempty(regionLabels); regionsLabels = num2cell([1:size(Adj,1)])'; end % generate a set of numbered labels
if isempty(subjectmask); subjectmask = ones(1,size(Adj,3)); end
if isempty(filelist); for i = 1:size(Adj,3); filelist{i} = strcat('file_',num2str(i));end; end

wResult.group = groups;
wResult.regionLabels = regionLabels;
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Weighted Analyses %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Weighted networks
nSubjects = size(Adj,3);
h = waitbar(0,'Running analysis on weighted networks');
for i = 1:nSubjects
    waitbar(i/nSubjects);
    A = squeeze(Adj(:,:,i));
    A = (A+A')./2.;
    A = threshold_proportional(A,percentage);
    
    % create non-normalized output
    wResult.deg(i,:) = degrees_und(A); %degree
    wResult.dens(i,:) = density_und(A); %density
    wResult.cpl(i,:) = charpath(distance_wei(A)); %characteristic path length
    wResult.trans(i,:) = transitivity_wu(A); %transitivity
    wResult.strength(i,:) = strengths_und(A); %strength
    
    [M Q] = modularityConsensusFun(A,1,nRand); %optimized clustering coefficient
    wResult.M(i,:) = M;
    wResult.clustcoeff(i,:) = mean(clustering_coef_bu(A));
    wResult.part(i,:) = participation_coef(A,M);
    
    for iR = 1:nRand
        R = randmio_und_connected(A, 10);% create a random matrix from original
        cpl_random(iR,:) = charpath(distance_wei(R));
        clust_random(iR,:) = mean(clustering_coef_wu(R));
        trans_random(iR,:) = transitivity_bu(R);
    end
    CplRand = mean(cpl_random(find(not(isinf(cpl_random)))));
    ClustRand = mean(clust_random,1);
    TransRand = mean(trans_random,1);
    
    wResult.Norm.cpl(i,:) = wResult.cpl(i,:)/CplRand;
    wResult.Norm.clustcoeff(i,:) = wResult.clustcoeff(i,:)/ClustRand;
    wResult.Norm.trans(i,:) = wResult.trans(i,:)/TransRand;
    
    [wResult.mask(i,:), wResult.net(:,:,i)] = getTop(wResult.deg(i,:),A,percentage);
end
close(h)


%% generate groupwise tables
for g = 1:length(unique(groups));
    varname = strcat('T',num2str(g-1));
    wResult.(varname) = table(squeeze(wResult.deg((groups == g-1),:)),squeeze(wResult.strength((groups == g-1),:)),...
        'VariableNames',{'degree','strength'},...
        'RowNames', filelist(groups == g-1)');   
end

% get some dimensions for subplots based on groupsize
dim1 = length(unique(groups));
dim2 = ceil(sqrt(dim1));

if PlotGlobal == 1
    figure;
    subplot(2,2,1); boxplot((squeeze(Result.strength(:,plotCost,:))'),groups,'colorgroup',1:length(unique(groups)));title('strength');legend;
    subplot(2,2,2); boxplot((squeeze(Result.cpl(:,plotCost,:))'),groups,'colorgroup',1:length(unique(groups)));title('characteristic path length');
    subplot(2,2,3); boxplot((squeeze(Result.deg(:,plotCost,:))'),groups,'colorgroup',1:length(unique(groups)));title('degree');
    subplot(2,2,4); boxplot((squeeze(Result.trans(:,plotCost,:))'),groups,'colorgroup',1:length(unique(groups)));title('transitivity');
    
    ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
    text(0.5, 1,'\bf Weighted Networks','HorizontalAlignment','center','VerticalAlignment', 'top');
    
end

if PlotLocal == 1
    figure; hold on;
    
    for p = 1:dim1
        subplot(dim1,1,p); 
            bar(mean(wResult.strength(groups == (p-1),:)),...
            'FaceColor',[0 .5 .5],'EdgeColor',[0 .9 .9],'LineWidth',1); set(gca,'XTick',1:1:(length(regionLabels)),...
            'XLim',[0 (length(regionLabels)+1)],'XTickLabel',regionLabels, 'XTickLabelRotation',90, 'Fontsize',...
            10); ylabel(num2str(p-1)); title('Strength');
    end
    hold off
end

if PlotMatrices == 1
    figure; hold on;
    for p = 1:dim1
        mat = zscore(squeeze(mean(Adj(:,:,(groups == p-1)),3)));
        
        subplot(dim1,ceil(dim1/dim2),p);
        imagesc(mat);            % Create a colored plot of the matrix values
        title(strcat('Mean connection probability for Group: ',num2str(p-1)));
        colormap(parula);  % Change the colormap to gray (so higher values are
        %#   black and lower values are white)
        c = colorbar; ylabel(c,'connection probability z-score ')
        
        %         textStrings = num2str(mat(:),'%0.2f');  %# Create strings from the matrix values
        %         textStrings = strtrim(cellstr(textStrings));  %# Remove any space padding
        %         [x,y] = meshgrid(1:size(mat,1));   %# Create x and y coordinates for the strings
        %         hStrings = text(x(:),y(:),textStrings(:),...      %# Plot the strings
        %             'HorizontalAlignment','center');
        %         midValue = mean(get(gca,'CLim'));  %# Get the middle value of the color range
        %         textColors = repmat(mat(:) > midValue,1,3);  %# Choose white or black for the
        %         %#   text color of the strings so
        %         %#   they can be easily seen over
        %         %#   the background color
        %         set(hStrings,{'Color'},num2cell(textColors,2));  %# Change the text colors
        %
        %         set(gca,'XTick',1:size(mat,1),...                         %# Change the axes tick marks
        %             'XTickLabel',regionLabels,...  %#   and tick labels
        %             'YTick',1:size(mat,1),...
        %             'YTickLabel',regionLabels,...
        %             'TickLength',[0 0]);
        
        
    end
    hold off
    
end

