function [Result p] = runStats(ResultIn,varargin)
%
% For example: [Results p] = runStats('ResultsIn', nResults);
%              [Results p] = runStats('ResultsIn', Result, 'cost', 6, 'nPerm', 100, 'PHocT', 0.1, 'Local', 1, 'localmetric','strength');
%              [Results p] = runStats('ResultsIn', Result, 'cost', 6, 'nPerm', 100, 'PHocT', 0.1, 'Global', 1, 'globalmetric', 'Sigma',);
%
% ---------------------------- INPUT ----------------------------
%
% ResultsIn             -       Results structure obtained from running graph analyses
%
% --------------------- OPTIONAL ARGUMENTS ----------------------
% groups                -       Grouping vector
% globalmetric          -       Which global metric you want to test (should correspond to a field in the result structure and can only have one vector per subject)
% localmetric           -       Which local metric you want to test (should correspond to a field in the result structure)
% covariates            -       Optional set of covariates to include (NB: HAS NOT BEEN TESTED YET)
% nPerm                 -       Number of permutation for the post-hoc test
% Global                -       Switch on global tests (default is on)
% Local                 -       Switch on local tests (defeault is on)
% PHoCT                 -       Threshold for conducting post-hoc tests
% ---------------------------- OUTPUT ----------------------------
%
% Result                -       Result structure containing fields corresponding to the requested localmetric and globalmetric inputs
%                               Should also contain a 'pairs' field that lists the pairs of groupwise comparisons for post-hoc tests
% p                     -       Inputparser structure listing all the inputs

groupdefault = round(rand(1,100)*3)'; %note that this now assumes a group size of 100
% input parsing settings
p = inputParser;
p.CaseSensitive = true;
p.Parameters;
p.Results;
p.KeepUnmatched = true;
% set the desired and optional input arguments
addRequired(p,'ResultsIn',@isstruct);
addOptional(p,'groups',groupdefault,@isnumeric);
addOptional(p,'cost',6,@isnumeric);
addOptional(p,'globalmetric','dens',@isstr);
addOptional(p,'localmetric','deg',@isstr);
addOptional(p,'covariates',[],@isnumeric);
addOptional(p,'nPerm',500,@isnumeric);
addOptional(p,'Global',0,@isnumeric);
addOptional(p,'Local',0,@isnumeric);
addOptional(p,'PHocT',0.1,@isnumeric);
% parse the input
parse(p,varargin{:});
% then set/get all the inputs out of this structure
ResultIn = p.Results.ResultsIn; groups = p.Results.groups; globalmetric = p.Results.globalmetric; covariates = p.Results.covariates;
nPerm = p.Results.nPerm; Global = p.Results.Global; Local = p.Results.Local;
localmetric = p.Results.localmetric; PHocT = p.Results.PHocT; cost = p.Results.cost;

%% get the number of unique group numbers
%% NB this assumes the groups are numbered sequentially
numgroups = length(unique(groups));
[a,b] = meshgrid([1:numgroups], [1:numgroups]);
Result.pairs = [a(:) b(:)];



%% Compute stats on global metrics
if Global == 1
    
    %% get the relevant cost point for that metrix
    ResultIn = squeeze(ResultIn.(globalmetric)(:,cost,:));
    
    for i = 1:numgroups
        mask = (groups == i-1);
        group{i} = ResultIn(mask,:);
    end
    
    if isempty(covariates)
        %%%%%% ANOVA %%%%%
        data = ResultIn;
        [Result.Global.(globalmetric).p,Result.Global.(globalmetric).tbl,Result.Global.(globalmetric).Astats] = anova1(data,groups,'off');
        
        fprintf('P-Value for %s is %d \n', globalmetric, Result.Global.(globalmetric).p);
        %%%%%% posthoc tests %%%%%
        if Result.Global.(globalmetric).p < PHocT
            for i = 1:length(Result.pairs)
                Result.Global.(globalmetric).permuted(i,:) = permutation_2tailed(group{Result.pairs(i,1)},group{Result.pairs(i,2)},nPerm);
                Result.Global.(globalmetric).meandiff(i,:) = mean(group{Result.pairs(i,1)}) - mean(group{Result.pairs(i,2)});
            end
        else
            Result.Local.(globalmetric).permuted(i,:) = NaN;
            Result.Local.(globalmetric).meandiff(i,:) = NaN;
        end
    else
        %%%%%% ANCOVA %%%%%
        data = ResultIn.(globalmetric);
        [Result.Global.(globalmetric).h,Result.Global.(globalmetric).atab,Result.Global.(globalmetric).ctab,Result.Global.(globalmetric).Cstats] = aoctool(ResultIn.(globalmetric),covariates,groups);
    end
end

if Local == 1
    
    %% get the relevant cost point for that metric
    ResultIn = squeeze(ResultIn.(localmetric)(:,cost,:));
    
    for i = 1:numgroups
        mask = (groups == i-1);
        group{i} = ResultIn(mask,:);
    end
    
    for inode = 1:size(ResultIn,2);
        %%%%%% ANOVA %%%%%
        data = ResultIn(:,inode);
        [Result.Local.(localmetric).p(inode),Result.Local.(localmetric).tbl{inode},Result.Local.(localmetric).Astats{inode}] = anova1(data,groups,'off');
        %%%%%% posthoc tests %%%%%
        if Result.Local.(localmetric).p(inode) < PHocT
            for i = 1:length(Result.pairs)
                Result.Local.(localmetric).permuted(i,inode,:) = permutation_2tailed(group{Result.pairs(i,1)}(:,inode),group{Result.pairs(i,2)}(:,inode),nPerm);
                Result.Local.(localmetric).meandiff(i,inode,:) = mean(group{Result.pairs(i,1)}(:,inode)) - mean(group{Result.pairs(i,2)}(:,inode));
            end
        else
            Result.Local.(localmetric).permuted(i,inode,:) = NaN;
            Result.Local.(localmetric).meandiff(i,inode,:) = NaN;
        end
    end
    
    plot(Result.Local.(localmetric).p);xlabel('nodes'); ylabel('p-value'); title(strcat('nodal permutation p-values for:',' ',localmetric));
    hline = refline(0,0.05);hline.Color = 'r';hline.LineStyle = '--';
end
end
