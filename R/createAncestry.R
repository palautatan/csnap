#' createAncestry()
#'
#' Outputs a list[[]] that corresponds with the phylogenetic tree structure via library(ape).
#' @param tre_file A path to your .tre output from RevBayes
#' @param outputlist If you want the function to output a list
#' @param df If you want the function to output a dataframe
#' @keywords ancestry
#' @export

#  *                        #
#         ANCESTRY          #
#                        *  #


createAncestry = function(tre_file, outputlist=TRUE, df=FALSE) {
  
  # * DEFINITIONS
  parsed_tree = phylotate::read_annotated(tre_file, format="nexus")
  num_tips    = length(parsed_tree$tip.label)
  edges       = parsed_tree$edge
  
  # * GET SPECIES AND LAST BRANCHES ON TREE
  species           = parsed_tree$tip.label
  last_branches     = edges[which(edges[,2] <= num_tips),]
  unique_branches_1 = unique(last_branches[,1])
  
  # * START MATRIX FOR CONVENIENCE
  merge_mat = cbind(last_branches, species)
  colnames(merge_mat)[1:2] = c("branches_1", "tip_index")
  
  
  # * GET NEXT BRANCHES
  next_branches_connected = intersect(edges[,2], last_branches[,1])
  next_branches_rows      = edges[which(edges[,2] %in% next_branches_connected),]
  
  colnames(next_branches_rows) = c("branches_2", "branches_1")
  merge_mat = merge(next_branches_rows, merge_mat, by="branches_1", all.y=TRUE)
  
  
  # * 1. BUILDING ANCESTRY DF ("MERGE MAT")
  next_branches_connected = intersect(edges[,2], merge_mat[,2])
  starting_index          = 3
  
  while (length(next_branches_connected)!=0) {
    # * CHOOSE THE BRANCHES' ROWS
    next_branches_rows    = edges[which(edges[,2] %in% next_branches_connected),]
    if (length(next_branches_rows)==2) {
      next_branches_rows    = t(data.frame(next_branches_rows))
    }
    
    
    # * MERGE BASED ON COLUMN NAME
    name_1 = paste0("branches_", starting_index)
    name_2 = paste0("branches_", (starting_index-1))
    colnames(next_branches_rows) = c(name_1, name_2)
    merge_mat = merge(next_branches_rows, merge_mat, by=name_2, all.y=TRUE)
    
    # * UPDATE VALUES
    starting_index = starting_index + 1
    next_branches_connected = intersect(edges[,2], merge_mat[,2])
  }
  
  merge_mat = merge_mat[c(2,1,3:ncol(merge_mat))]
  
  
  # * 2. BUILDING ANCESTRY LIST
  num_interior_nodes = (num_tips+1):max(edges)
  ancestry = list()
  
  for (this_node in num_interior_nodes) {
    rows = which(merge_mat == this_node, arr.ind=T)
    if (nrow(rows)>1) {
      rows = rows[,1]
    }
    rows = merge_mat[rows,]
    remove_species = rows[,-ncol(rows)]
    
    values_under = unlist(lapply(1:nrow(remove_species), function(ix) {
      this_row = remove_species[ix,]
      location = which(this_row == this_node)
      collect  = this_row[(location+1):length(this_row)]
      as.numeric(as.character(as.matrix(collect)))
    }))
    # values_under = as.vector(unlist(rows[,-ncol(rows)]))
    values_under = unique(values_under[!is.na(values_under)])
    
    rows = lapply(values_under, function(this_node) which(merge_mat == this_node, arr.ind=T))
    if (class(rows)=="list") {
      these_rows = do.call(rbind, rows)
      these_rows = these_rows[,1]
    } else {
      if (class(rows)=="matrix") {
        these_rows = rows[,1]
      }
    }
    
    all_these_species = unique(as.vector(merge_mat[these_rows,]$species))
    
    info_for_node = list(nodes_under=values_under, species_under=all_these_species)
    
    ancestry[[this_node]] = info_for_node
  }
  
  # * RETURN VALUES
  if (df==TRUE & outputlist==TRUE) {
    return(list(ancestry_df=merge_mat, ancestry_list=ancestry))
  }
  
  if (df==TRUE) {
    return(merge_mat)
  }
  
  if (outputlist==TRUE) {
    return(ancestry)
  }
}
