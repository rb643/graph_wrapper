function nodeff = efficiency_nodal(A,bin)

% input
% A = binary or weighted adjacency matrix
% bin = "binary" = optional flag: T if binary (default), F if weighted
% 
% output
% nodeff = nodal efficiency
% 
% author: Frantisek Vasa

if nargin < 2
    bin = true;
end

%Co(1:n+1:n*n)=1;    % set diagonal to ones

if bin
    D = distance_bin(A);
else
    % prepare "connection-length" matrix ('topological distance' = 1/'edge weight')
    CL = A;
    CL(CL>0) = 1./CL(CL>0);
    % compute topological distance matrix
    D = distance_wei(CL);
end

% Nodal efficiency
invD = 1./D;
invD(isinf(invD)) = NaN; % ignores disconnected nodes; not an issue with MST
nodeff = nansum(invD) / (length(A)-1);

end