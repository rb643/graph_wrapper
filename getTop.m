function [mask, net] = getTop(deg_data,net,percentage)
% get the top x% of nodes from a degree distribution and return the degree mask and matrix

Ms = sort(deg_data(:),'descend');                       % Sort Descending
Result = Ms(1:ceil(length(Ms)*percentage));                % Desired Output
minimalV = min(Result);

deg_data(deg_data<minimalV)=0;
deg_data(deg_data>=minimalV)=1;                             % Threshold

mask = logical(deg_data);

net(~mask,~mask)=0; %Remove all the connections between no hubs
net(mask,~mask)=0;    %Remove all connections between hubs and no hubs
net(~mask,mask)=0;    %The same

mask = mask';

end


