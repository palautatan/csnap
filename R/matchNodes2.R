#' matchNodes2()
#'
#' Outputs a matrix[,] that shows which nodes match on species
#' @param list_1 A list object from createAncestry()
#' @param list_2 Another list object from createAncestry()
#' @keywords ancestry
#' @export

# *                          #
#       MATCH ANCESTRY       #
#                          * #

matchNodes2 = function(ancestry_1, ancestry_2, intersection, num_tips) {
  ix_1 = which(ancestry_1 %in% intersection)
  ix_1 = ix_1[ix_1>num_tips]

  ix_2 = which(ancestry_2 %in% intersection)
  ix_2 = ix_2[ix_2>num_tips]

  species_1 = ancestry_1[ix_1]
  species_2 = ancestry_2[ix_2]

  match_df = matchNodes(species_1, species_2)

  order_1 = ix_1[match_df[,1]]
  order_2 = ix_2[match_df[,2]]
  cbind(order_1, order_2)
}
