function [Result p] = graph_binary(Matrix,varargin)
%
% For example:
% [Result p] = graph_binary('Matrix',adjacencymatrix);
% [Result p] = graph_binary('Matrix',adjacencymatrix,'regionLabels',regions,...
%                                   1'nRand',100, 'subjectmask',include, 'groups',...
%                                   groups, 'cost', [0:0.01:0.3]);
%
% ---------------------------- INPUT ----------------------------
%
% Matrix            -       3D connectivity matrix(node*node*subject)
%
% --------------------- OPTIONAL ARGUMENTS ----------------------
%
% regionLabels      -       A cell structure containing names of all nodes in your matrix (at this point its not really used)
% nRand             -       number of randomizations used in random networks for normalization (default = 100)
% groups            -       1D vector with grouplabels  (default = random)
% cost              -       range of costs (default = [0:0.02:0.3])
% percentage        -       threshold for getting the top degree nodes (default = 0.1)
% PlotMatrices      -       Logical to set if we want to plot group average matrices (default = 0)
% PlotLocal         -       Logical to set if we want to plot local metrics (default = 0)
% PlotGlobal        -       Logical to set if we want to plot global metrics (default = 0)
% PlotCost          -       Costpoint at which to plot (default = 5)
%                           NOTE: this will fail if there are less than 5 costpoints
%
% ---------------------------- OUTPUT ----------------------------
% Result            -       Graph metrics from binary matrices
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
addOptional(p,'nRand',10,@isnumeric);
addOptional(p,'PlotLocal',0,@isnumeric);
addOptional(p,'PlotGlobal',0,@isnumeric);
addOptional(p,'PlotMatrices',0,@isnumeric);
addOptional(p,'groups',[], @isnumeric);
addOptional(p,'PlotCost',5,@isnumeric);
addOptional(p,'cost',[0:0.05:0.3], @isnumeric);

% parse the input
parse(p,varargin{:});
% then set/get all the inputs out of this structure
nRand = p.Results.nRand; regionLabels = p.Results.regionLabels; PlotLocal = p.Results.PlotLocal; PlotGlobal = p.Results.PlotGlobal; groups = p.Results.groups; Adj = p.Results.Matrix;
cost = p.Results.cost; PlotMatrices = p.Results.PlotMatrices; ; plotCost = p.Results.PlotCost;

% set the larger defaults in case they are not specified
if isempty(groups); groups = round(rand(1,size(Adj,3))*3)'; end % generate some random numbers
if isempty(regionLabels); regionLabels = num2cell([1:size(Adj,1)])'; end % generate a set of numbered labels

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Unweighted/Binary %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate degree, density, length, and clustering for each subject
% Everything here is done on a single-level, only use prevalence-thresholds
% when looking at group average connectomes
h = waitbar(0,'Running analysis on binary graphs across cost range');
nSubjects = size(Adj,3);
Result.cost = cost;

for i = 1:nSubjects
    fprintf('Working on subject %d \n', i);
    
    A = squeeze(Adj(:,:,i));
    A = (A+A')./2.;
    
    % threshold
    for c = 1:length(cost)
        waitbar(c/length(cost));
        A = mst_threshold(A,cost(c),'T');
        % get some random matrices for normalization later
        for iR = 1:nRand
            R = randmio_und_connected(A, 10);
            RichRand(iR,:) = rich_club_bu(R,50);
            cpl_random(iR,:) = charpath(distance_bin(R));
            clust_random(iR,:) = mean(clustering_coef_bu(R));
            trans_random(iR,:) = transitivity_bu(R);
        end
        
        RichRand = mean(RichRand,1);
        CplRand = mean(cpl_random(find(not(isinf(cpl_random)))));
        ClustCoeffRand = mean(clust_random,1);
        TransRand = mean(trans_random,1);
        
        % create non-normalized output
        Result.deg(i,c,:) = degrees_und(A); % degree
        Result.strength(i,c,:) = strengths_und(A); % strength
        Result.dens(i,c,:) = density_und(A); %density
        Result.cpl(i,c,:) = charpath(distance_bin(A)); %characteristic path length
        [M Q] = modularityConsensusFun(A,1,nRand); %optimized clustering coefficient
        Result.M(i,c,:) = M;
        Result.clustcoeff(i,c,:) = mean(clustering_coef_bu(A));
        Result.rich(i,c,:) = rich_club_bu(A,50); %rich-club coefficient
        Result.trans(i,c,:) = transitivity_bu(A); %transitivity
        Result.assor(i,c,:) = assortativity_bin(A,0); %assortativity
        Result.effN(i,c,:) = efficiency_nodal(A); %efficiency
        Result.part(i,c,:) = participation_coef(A,M);
        
        % create a random matrix for normalization for small-worldness
        R = randmio_und(A,nRand);
        Crand = mean(clustering_coef_bu(R)); % get random clustering
        Lrand = charpath(distance_bin(R)); % get path length from random matrix
        Result.Sigma(i,c,:)=(Result.clustcoeff(i,c,:)./Crand)./(Result.cpl(i,c,:)./Lrand); % get small world coefficient
        
        
        Result.Norm.rich(i,c,:) = squeeze(Result.rich(i,c,:))./RichRand';
        Result.Norm.cpl(i,c,:) = Result.cpl(i,c,:)/CplRand;
        Result.Norm.clustcoeff(i,c,:) = Result.clustcoeff(i,c,:)/ClustCoeffRand;
        Result.Norm.trans(i,c,:) = Result.trans(i,c,:)/TransRand;
        
    end
end
close(h);

% get some dimensions for subplots based on groupsize
dim1 = length(unique(groups));
dim2 = ceil(sqrt(dim1));

if PlotGlobal == 1
figure;
    subplot(2,2,1); boxplot((squeeze(Result.strength(:,plotCost,:))'),groups,'colorgroup',1:length(unique(groups)));title('strength');legend;
    subplot(2,2,2); boxplot((squeeze(Result.cpl(:,plotCost,:))'),groups,'colorgroup',1:length(unique(groups)));title('characteristic path length');
    subplot(2,2,3); boxplot((squeeze(Result.clustcoeff(:,plotCost,:))'),groups,'colorgroup',1:length(unique(groups)));title('clustering');
    subplot(2,2,4); boxplot((squeeze(Result.trans(:,plotCost,:))'),groups,'colorgroup',1:length(unique(groups)));title('transitivity');
    
    ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
    text(0.5, 1,'\bf Weighted Networks','HorizontalAlignment','center','VerticalAlignment', 'top');
       
end


if PlotLocal == 1
    figure; hold on;
    
    for p = 1:dim1
        subplot(dim1,1,p); bar(mean(Result.strength(groups == (p-1),:)),...
            'FaceColor',[0 .5 .5],'EdgeColor',[0 .9 .9],'LineWidth',1); ...
            set(gca,'XTick',1:1:(length(regionLabels)),'XLim',[0 (length(regionLabels)+1)],...
            'XTickLabel',regionLabels, 'XTickLabelRotation',90, 'Fontsize', 10); ...
            ylabel(num2str(p-1)); title('Strength');
    end
    hold off
end


if PlotMatrices == 1
    dim1 = length(unique(groups));
    dim2 = ceil(sqrt(dim1));
    
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