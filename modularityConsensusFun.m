function [CiConsensus consensus]=modularityConsensusFun(net,gamma,numRand)
n=size(net,1);
consensus=zeros(n);
ind_mod=zeros(n);
m=[];
    for ir=1:numRand
        ir;
        [Ci Q]=community_louvain(net,gamma);
        nmod=max(Ci);
        ind_mod=zeros(n);
        for im=1:nmod
            posMod=find(Ci==im);
            ind_mod(posMod,posMod)=1;
        end
        consensus=consensus+ind_mod;
    end
    consensus_norm=consensus./numRand;
    [CiConsensus Q]=community_louvain(consensus_norm,gamma);
end
