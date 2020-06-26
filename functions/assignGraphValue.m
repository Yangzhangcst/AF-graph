% In this function, the global graph value in the set named global_nodes 
% W: fused graph, initialized by the adjacency graph;
% W_given: the given graph needed to be fused 
% gloabl_set: store the index of W_given could be fused;
% local_set: store the index of initilized W,i.e. adjacency graph;

function W=assignGraphValue(W,W_give,global_set)
[nonzero_row nonzero_col]=find(W_give>0);
[local_row local_col]=find(W>0);
for i=1:length(global_set)
    % find the node index denoted as the col number
    index_col_global=find(nonzero_col==global_set(i));
    for j=1:length(index_col_global)
        if length(nonzero_row)<nonzero_col(index_col_global(j))...
                ||length(nonzero_col)<index_col_global(j)
            index_col_global = []; break; 
        end
        % check all its connected neighbors in the global graph
         global_neighbors=nonzero_row(nonzero_col(index_col_global(j)));
         % check this node's neighbors in the local graph
         local_neighbors=local_row(nonzero_col(index_col_global(j)));
         % fuse the neighborhood sets
         new_row=unique([global_neighbors local_neighbors]);
         new_col=nonzero_col(index_col_global(j));
    end
          %[intersect_value,i_global, i_local]=intersect(global_neighbors,local_neighbors);
         % remove all  new_row(i_remove)=[];
          %new_col=nonzero_col(index_col);%repmat(local_nodes(i),length(index_col),1);
   if isempty(index_col_global), continue; end
   clear index_col_global;
   W(new_row,new_col)=W_give(new_row,new_col);
        
end
end