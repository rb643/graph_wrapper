function [Result p] = graph_binary(Matrix,varargin)
%
% For example:
% [Result p] = graph_binary('Matrix',adjacencymatrix);
% [Result p] = graph_binary('Matrix',adjacencymatrix,'regionLabels',regions,'prev',0.75,'PlotLocal',1,'PlotGlobal',...
%                                   1,'nos',10, 'subjectmask',include, 'groups', groups, 'cost', [0:0.01:0.3]);
%
% ---------------------------- INPUT ----------------------------
%
% Matrix            -       3D connectivity matrix(node*node*subject)
%
% --------------------- OPTIONAL ARGUMENTS ----------------------
%
% regionLabels      -       A cell structure containing names of all nodes in your matrix (at this point its not really used)
% nRand             -       number of randomizations used in random networks for normalization (default = 10)
% prev              -       prevalence for prevalence weighted matrices  (default = 0.75)
% groups            -       1D vector with grouplabels  (default = random)
% cost              -       range of costs (default = [0:0.05:0.3])
% percentage        -       threshold for getting the top degree nodes (default = 0.1)
%
% ---------------------------- OUTPUT ----------------------------
% Result            -       Graph metrics from binary matrices 
% p                 -       Also return all the input settings to check
 
%% check all the inputs and if they do not exist then revert to default settings
% set the larger defaults in case they are not specified
groupdefault = round(rand(1,size(Matrix,3))*3)'; % generate some random numbers
regionsdefault = num2cell([1:size(Matrix,1)])'; % generate a set of numbered labels
subjectsdefault = ones(1,size(Matrix,3));
% input parsing settings
p = inputParser;
p.CaseSensitive = true;
p.Parameters;
p.Results;
p.KeepUnmatched = true;
% set the desired and optional input arguments
addRequired(p,'Matrix',@isnumeric);
addOptional(p,'regionLabels',regionsdefault,@iscell);
addOptional(p,'nRand',10,@isnumeric);
addOptional(p,'PlotLocal',0,@isnumeric);
addOptional(p,'PlotGlobal',0,@isnumeric);
addOptional(p,'PlotMatrices',0,@isnumeric);
addOptional(p,'prev',0.75,@isnumeric);
addOptional(p,'groups',groupdefault, @isnumeric);
addOptional(p,'cost',[0:0.05:0.3], @isnumeric);
addOptional(p,'percentage',0.1, @isnumeric);

% parse the input
parse(p,varargin{:});
% then set/get all the inputs out of this structure
nRand = p.Results.nRand; regionLabels = p.Results.regionLabels; PlotLocal = p.Results.PlotLocal; PlotGlobal = p.Results.PlotGlobal; prev = p.Results.prev; groups = p.Results.groups; Adj = p.Results.Matrix;
cost = p.Results.cost; PlotMatrices = p.Results.PlotMatrices; percentage = p.Results.percentage;
 
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
    th=prctile((A(:)),90); 
    % threshold
    for c = 1:length(cost)
    waitbar(c/length(cost));
    th=prctile((A(:)),100-(cost(c)*100));
    A = (A>th);
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
