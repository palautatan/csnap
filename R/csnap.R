#' csnap()
#'
#' Create a Compare Shared Node Ages Plot (CSNAP)
#' @param tre_files_path A path to where your .tre files (assumes you have AIC, Saturated, BIC, Gene, and Uniform). Future edits will make this generalized.
#' @keywords ancestry
#' @export

csnap = function(tre_files_path) {
  # * READ IN FILES
  tre_files = list.files(tre_files_path, pattern=".tre$", full.names=TRUE)
  aic_anc   = createAncestry(tre_files[1])
  sat_anc   = createAncestry(tre_files[2])
  bic_anc   = createAncestry(tre_files[3])
  gen_anc   = createAncestry(tre_files[4])
  uni_anc   = createAncestry(tre_files[5])
  
  # * CREATE SUBLISTS
  aic_species = sapply(1:length(aic_anc), function(x) aic_anc[[x]]$species_under)
  sat_species = sapply(1:length(sat_anc), function(x) sat_anc[[x]]$species_under)
  bic_species = sapply(1:length(bic_anc), function(x) bic_anc[[x]]$species_under)
  gen_species = sapply(1:length(gen_anc), function(x) gen_anc[[x]]$species_under)
  uni_species = sapply(1:length(uni_anc), function(x) uni_anc[[x]]$species_under)
  
  
  # * SORT SUBLISTS
  aic_species = lapply(aic_species, sort)
  sat_species = lapply(sat_species, sort)
  bic_species = lapply(bic_species, sort)
  gen_species = lapply(gen_species, sort)
  uni_species = lapply(uni_species, sort)
  
  
  # * INTERSECT TO GET SHARED NODES
  aic_sat = intersect(aic_species, sat_species)
  
  
  shared_nodes = intersect(uni_species, intersect(gen_species, intersect(bic_species, intersect(aic_species, sat_species))))
  shared_nodes[sapply(shared_nodes, is.null)] = NULL
  
  
  # * COLLECT NODE NUMBERS
  aic_num = which(aic_species %in% shared_nodes)
  sat_num = which(sat_species %in% shared_nodes)
  bic_num = which(bic_species %in% shared_nodes)
  gen_num = which(gen_species %in% shared_nodes)
  uni_num = which(uni_species %in% shared_nodes)
  
  
  # * GET THE MATCHES
  aic_match = aic_species[aic_num]
  sat_match = sat_species[sat_num]
  bic_match = bic_species[bic_num]
  gen_match = gen_species[gen_num]
  uni_match = uni_species[uni_num]
  
  
  # * MATCH AND MERGE ON AIC
  table_1  = matchNodes(aic_match, sat_match)
  table_2  = matchNodes(aic_match, bic_match)
  table_3  = matchNodes(aic_match, gen_match)
  table_4  = matchNodes(aic_match, uni_match)
  match_df = data.frame(table_1, table_2[,2], table_3[,2], table_4[,2])
  colnames(match_df) = c("aic", "sat", "bic", "gen", "uni")
  
  
  
  # * ATTACH ACTUAL NODE VALUES TO NEW DF
  aic_col = aic_num[match_df$aic]
  sat_col = sat_num[match_df$sat]
  bic_col = bic_num[match_df$bic]
  gen_col = gen_num[match_df$gen]
  uni_col = uni_num[match_df$uni]
  node_df = data.frame(aic_col, sat_col, bic_col, gen_col, uni_col)
  colnames(node_df) = colnames(match_df)
  
  aic_anc   = addComments(tre_files[1], aic_anc)
  sat_anc   = addComments(tre_files[2], sat_anc)
  bic_anc   = addComments(tre_files[3], bic_anc)
  gen_anc   = addComments(tre_files[4], gen_anc)
  uni_anc   = addComments(tre_files[5], uni_anc)
  
  
  # * GET AGE HPD'S FOR THE NODES
  aic_age = sapply(node_df$aic, function(x) aic_anc[[x]]$age_hpd)
  sat_age = sapply(node_df$sat, function(x) sat_anc[[x]]$age_hpd)
  bic_age = sapply(node_df$bic, function(x) bic_anc[[x]]$age_hpd)
  gen_age = sapply(node_df$gen, function(x) gen_anc[[x]]$age_hpd)
  uni_age = sapply(node_df$uni, function(x) uni_anc[[x]]$age_hpd)
  age_df  = t(data.frame(aic_age, sat_age, bic_age, gen_age, uni_age))
  
  
  # * CREATE NUMERIC PLOTTABLES
  aic_hpd = matrix(as.numeric(unlist(strsplit(age_df[1,], ","))), ncol=2, byrow=TRUE)
  aic_hpd = cbind(aic_hpd, (aic_hpd[,1] + aic_hpd[,2])/2)
  
  sat_hpd = matrix(as.numeric(unlist(strsplit(age_df[2,], ","))), ncol=2, byrow=TRUE)
  sat_hpd = cbind(sat_hpd, (sat_hpd[,1] + sat_hpd[,2])/2)
  
  bic_hpd = matrix(as.numeric(unlist(strsplit(age_df[3,], ","))), ncol=2, byrow=TRUE)
  bic_hpd = cbind(bic_hpd, (bic_hpd[,1] + bic_hpd[,2])/2)
  
  gen_hpd = matrix(as.numeric(unlist(strsplit(age_df[4,], ","))), ncol=2, byrow=TRUE)
  gen_hpd = cbind(gen_hpd, (gen_hpd[,1] + gen_hpd[,2])/2)
  
  uni_hpd = matrix(as.numeric(unlist(strsplit(age_df[5,], ","))), ncol=2, byrow=TRUE)
  uni_hpd = cbind(uni_hpd, (uni_hpd[,1] + uni_hpd[,2])/2)
  
  
  # * PLOT ALL PAIRWISE COMBINATIONS
  
  par(mfrow=c(5,5), mar=c(0, 0, 0, 0), oma = c(4, 4, 4, 1))
  
  # * SATURATED
  plot(sat_hpd[,3], sat_hpd[,3], pch=19, col="orange2", xlab="", ylab="", axes=FALSE)
  abline(a=0, b=1, col="orange2")
  box(lty=1, col="black")
  mtext("saturated", side=2, line=1, adj=0.5, cex=1)
  
  plot(sat_hpd[,3], aic_hpd[,3], pch=19, col="gray22", xlab="", ylab="", axes=FALSE)
  abline(a=0, b=1, col="snow3", lty=3)
  box(lty=1, col="black")
  
  plot(sat_hpd[,3], bic_hpd[,3], pch=19, col="gray22", xlab="", ylab="", axes=FALSE)
  abline(a=0, b=1, col="snow3", lty=3)
  box(lty=1, col="black")
  
  plot(sat_hpd[,3], gen_hpd[,3], pch=19, col="gray22", xlab="", ylab="", axes=FALSE)
  abline(a=0, b=1, col="snow3", lty=3)
  box(lty=1, col="black")
  
  plot(sat_hpd[,3], uni_hpd[,3], pch=19, col="gray22", xlab="", ylab="", axes=FALSE)
  abline(a=0, b=1, col="snow3", lty=3)
  box(lty=1, col="black")
  
  
  # * AIC
  
  plot(aic_hpd[,3], sat_hpd[,3], pch=19, col="gray22", xlab="", ylab="", axes=FALSE)
  abline(a=0, b=1, col="snow3", lty=3)
  box(lty=1, col="black")
  mtext("aic", side=2, line=1, adj=0.5, cex=1)
  
  
  plot(aic_hpd[,3], aic_hpd[,3], pch=19, col="orange2", xlab="", ylab="", axes=FALSE)
  abline(a=0, b=1, col="orange2")
  box(lty=1, col="black")
  
  plot(aic_hpd[,3], bic_hpd[,3], pch=19, col="gray22", xlab="", ylab="", axes=FALSE)
  abline(a=0, b=1, col="snow3", lty=3)
  box(lty=1, col="black")
  
  plot(aic_hpd[,3], gen_hpd[,3], pch=19, col="gray22", xlab="", ylab="", axes=FALSE)
  abline(a=0, b=1, col="snow3", lty=3)
  box(lty=1, col="black")
  
  plot(aic_hpd[,3], uni_hpd[,3], pch=19, col="gray22", xlab="", ylab="", axes=FALSE)
  abline(a=0, b=1, col="snow3", lty=3)
  box(lty=1, col="black")
  
  
  # * BIC
  
  plot(bic_hpd[,3], sat_hpd[,3], pch=19, col="gray22", xlab="", ylab="", axes=FALSE)
  abline(a=0, b=1, col="snow3", lty=3)
  box(lty=1, col="black")
  mtext("bic", side=2, line=1, adj=0.5, cex=1)
  
  
  plot(bic_hpd[,3], aic_hpd[,3], pch=19, col="gray22", xlab="", ylab="", axes=FALSE)
  abline(a=0, b=1, col="snow3", lty=3)
  box(lty=1, col="black")
  
  
  plot(bic_hpd[,3], bic_hpd[,3], pch=19, col="orange2", xlab="", ylab="", axes=FALSE)
  abline(a=0, b=1, col="orange2")
  box(lty=1, col="black")
  
  plot(bic_hpd[,3], gen_hpd[,3], pch=19, col="gray22", xlab="", ylab="", axes=FALSE)
  abline(a=0, b=1, col="snow3", lty=3)
  box(lty=1, col="black")
  
  plot(bic_hpd[,3], uni_hpd[,3], pch=19, col="gray22", xlab="", ylab="", axes=FALSE)
  abline(a=0, b=1, col="snow3", lty=3)
  box(lty=1, col="black")
  
  
  # * GENE
  
  plot(gen_hpd[,3], sat_hpd[,3], pch=19, col="gray22", xlab="", ylab="", axes=FALSE)
  abline(a=0, b=1, col="snow3", lty=3)
  box(lty=1, col="black")
  mtext("gene", side=2, line=1, adj=0.5, cex=1)
  
  
  plot(gen_hpd[,3], aic_hpd[,3], pch=19, col="gray22", xlab="", ylab="", axes=FALSE)
  abline(a=0, b=1, col="snow3", lty=3)
  box(lty=1, col="black")
  
  plot(gen_hpd[,3], bic_hpd[,3], pch=19, col="gray22", xlab="", ylab="", axes=FALSE)
  abline(a=0, b=1, col="snow3", lty=3)
  box(lty=1, col="black")
  
  plot(gen_hpd[,3], gen_hpd[,3], pch=19, col="orange2", xlab="", ylab="", axes=FALSE)
  abline(a=0, b=1, col="orange2")
  box(lty=1, col="black")
  
  plot(gen_hpd[,3], uni_hpd[,3], pch=19, col="gray22", xlab="", ylab="", axes=FALSE)
  abline(a=0, b=1, col="snow3", lty=3)
  box(lty=1, col="black")
  
  
  # * UNIFORM
  
  plot(uni_hpd[,3], sat_hpd[,3], pch=19, col="gray22", xlab="", ylab="", axes=FALSE)
  abline(a=0, b=1, col="snow3", lty=3)
  box(lty=1, col="black")
  mtext("uniform", side=2, line=1, adj=0.5, cex=1)
  mtext("saturated", side=1, line=1, adj=0.5, cex=1)
  
  
  plot(uni_hpd[,3], aic_hpd[,3], pch=19, col="gray22", xlab="", ylab="", axes=FALSE)
  abline(a=0, b=1, col="snow3", lty=3)
  box(lty=1, col="black")
  mtext("aic", side=1, line=1, adj=0.5, cex=1)
  
  
  plot(uni_hpd[,3], bic_hpd[,3], pch=19, col="gray22", xlab="", ylab="", axes=FALSE)
  abline(a=0, b=1, col="snow3", lty=3)
  box(lty=1, col="black")
  mtext("bic", side=1, line=1, adj=0.5, cex=1)
  
  
  plot(uni_hpd[,3], gen_hpd[,3], pch=19, col="gray22", xlab="", ylab="", axes=FALSE)
  abline(a=0, b=1, col="snow3", lty=3)
  box(lty=1, col="black")
  mtext("gene", side=1, line=1, adj=0.5, cex=1)
  
  
  plot(uni_hpd[,3], uni_hpd[,3], pch=19, col="orange2", xlab="", ylab="", axes=FALSE)
  abline(a=0, b=1, col="orange2")
  box(lty=1, col="black")
  mtext("uniform", side=1, line=1, adj=0.5, cex=1)
  
  
  title(main="Divergence Time Comparisons for Shared Nodes\nUnder 5 Partitioning Schemes", sub="117", cex.main=2, outer=TRUE)
  
}


