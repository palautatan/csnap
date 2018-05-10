#' addComments()
#'
#' Appends comments to an existing ancestry list from createAncestry()
#' @param tre_file A path to your .tre output from RevBayes
#' @param ancestry_list A list object from createAncestry()
#' @keywords ancestry
#' @export



# *                                    #
#       ADD COMMENTS TO ANCESTRY       #
#                                    * #

addComments = function(tre_file, ancestry_list) {
  
  # * DEFINITIONS
  parsed_tree = phylotate::read_annotated(tre_file, format="nexus")
  num_tips    = length(parsed_tree$tip.label)
  
  # * CLEAN UP NODE COMMENTS
  node_info = parsed_tree$node.comment
  node_info = gsub("\\{|\\}", "", node_info)
  node_info = strsplit(node_info, ",[a-z]")
  node_info = do.call(rbind, node_info)
  node_info = sapply(1:ncol(node_info), function(x) sapply(1:length(node_info[,x]), function(z) strsplit(node_info[z,x], "\\=")[[1]][2]))
  colnames(node_info) = c("index", "posterior", "age_hpd")
  
  # * MATCH TO EACH OTHER
  copy_ancestry = ancestry_list
  
  for (ix in ((num_tips+1):length(copy_ancestry))) {
    for (jx in 1:ncol(node_info)) {
      copy_ancestry[[ix]][colnames(node_info)[jx]] = node_info[ix,jx]
    }
  }
  
  return(copy_ancestry)
}
