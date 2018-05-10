#' matchNodes()
#'
#' Outputs a matrix[,] that shows which nodes match on species
#' @param list_1 A list object from createAncestry()
#' @param list_2 Another list object from createAncestry()
#' @keywords ancestry
#' @export

# *                          #
#       MATCH ANCESTRY       #
#                          * #

matchNodes = function(list_1, list_2) {
  list_1_nodes = c()
  list_2_nodes = c()
  
  for (ix in 1:length(list_1)) {
    this_list = list_1[[ix]]
    
    for (jx in 1:length(list_2)) {
      other_list = list_2[[jx]]
      
      if (identical(this_list, other_list)) {
        list_1_nodes = c(list_1_nodes, ix)
        list_2_nodes = c(list_2_nodes, jx)
      }
    }
  }
  return(cbind(list_1=list_1_nodes, list_2=list_2_nodes))
}