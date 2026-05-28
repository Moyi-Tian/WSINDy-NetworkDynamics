%{
Name:
    O2O_DesignNetwork_NodesAveDegK


Version:
    wessler
    2024 July 11
    update to: 2024 July 10
    changes:
        *now that running, making this and other network option functions


Description:
    *This algorithm makes a network of nodes with each node having an
    average degree=AveDeg_input (TotalEdges/TotalNodes=AveDeg_input)
    *It can be thought of as creating one of the Erdős–Rényi model graphs
    G(NumNodes,NumNodes*AveDeg_input), with the condition that the output
    network has NO self-loops and NO multi-edges


Inputs:
    *NumNodes
    *AveDeg_input
    
    


Outputs/does:
    *MatNetwork
    


Used by:
    


Uses:
    *NOTHING


NOTES:




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%==========================================================================
%--------------------------------------------------------------------------
%__________________________________________________________________________
%}

function MatNetwork=O2O_DesignNetwork_NodesAveDegK(NumNodes,AveDeg_input)



%==========================================================================
% make matrix for connections in network (all weights = 1)
%==========================================================================
%at index (i,j): 1 means node ID i connected to node ID j; 0 means i not connected to j
%the number of 1's is designed to make the ave degree of nodes to be input_AveDeg
NumElementsUpperMatrix=(NumNodes*NumNodes-NumNodes)/2;
NumEdges=round(NumNodes*AveDeg_input/2);
tempVec_Connections=zeros(NumElementsUpperMatrix,1);
tempVec_Connections(1:NumEdges)=1;
tempVec_Connections=tempVec_Connections(randperm(NumElementsUpperMatrix));
MatNetwork=zeros(NumNodes);
Index_tempVec=0;
for ii=1:NumNodes
    for jj=ii+1:NumNodes
        Index_tempVec=Index_tempVec+1;
        MatNetwork(ii,jj)=tempVec_Connections(Index_tempVec);
    end
end
MatNetwork=MatNetwork+MatNetwork';








end
