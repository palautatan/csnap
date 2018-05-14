#' csnapPairs()
#'
#' Create a Compare Shared Node Ages Plot (CSNAP) between each pairwise sharing (AIC and BIC share these nodes / Gene and AIC share different nodes)
#' @param tre_files_path A path to where your .tre files (assumes you have AIC, Saturated, BIC, Gene, and Uniform). Future edits will make this generalized.
#' @keywords ancestry
#' @export

csnapPairs = function(tre_files_path, study="") {
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
  aic_bic = intersect(aic_species, bic_species)
  aic_gen = intersect(aic_species, gen_species)
  aic_uni = intersect(aic_species, uni_species)

  sat_bic = intersect(sat_species, bic_species)
  sat_gen = intersect(sat_species, gen_species)
  sat_uni = intersect(sat_species, uni_species)

  bic_gen = intersect(bic_species, gen_species)
  bic_uni = intersect(bic_species, uni_species)

  gen_uni = intersect(gen_species, uni_species)

  intersections = list(aic_sat, aic_bic, aic_gen, aic_uni, sat_bic, sat_gen, sat_uni, bic_gen, bic_uni, gen_uni)
  intersections = lapply(intersections, function(x) x[-sapply(x, is.null)])

  num_tips = length(aic_species[sapply(aic_species, is.null)])
  total_nodes = length(aic_species[(num_tips+1):length(aic_species)])


  # * FRACTION SHARED
  shared_fractions = c()

  for (intersection in intersections) {
    length   = length(intersection)
    fraction = round(length / total_nodes, 4)
    shared_fractions = c(shared_fractions, fraction)
  }

  ranks = round(round((((shared_fractions/max(shared_fractions)))*10),2))
  rank_mat = cbind(shared_fractions, ranks)


  # * MATCH UP INTERSECTS
  aic_sat_ix = matchNodes2(aic_species, sat_species, aic_sat, num_tips)
  aic_bic_ix = matchNodes2(aic_species, bic_species, aic_bic, num_tips)
  aic_gen_ix = matchNodes2(aic_species, gen_species, aic_gen, num_tips)
  aic_uni_ix = matchNodes2(aic_species, uni_species, aic_uni, num_tips)

  sat_bic_ix = matchNodes2(sat_species, bic_species, sat_bic, num_tips)
  sat_gen_ix = matchNodes2(sat_species, gen_species, sat_gen, num_tips)
  sat_uni_ix = matchNodes2(sat_species, uni_species, sat_uni, num_tips)

  bic_gen_ix = matchNodes2(bic_species, gen_species, bic_gen, num_tips)
  bic_uni_ix = matchNodes2(bic_species, uni_species, bic_uni, num_tips)

  gen_uni_ix = matchNodes2(gen_species, uni_species, gen_uni, num_tips)



  # * APPEND COMMENTS
  aic_anc   = addComments(tre_files[1], aic_anc)
  sat_anc   = addComments(tre_files[2], sat_anc)
  bic_anc   = addComments(tre_files[3], bic_anc)
  gen_anc   = addComments(tre_files[4], gen_anc)
  uni_anc   = addComments(tre_files[5], uni_anc)


  # * GET AGE HPD'S FOR ALL NODES
  aic_age = unlist(sapply(1:length(aic_species), function(x) aic_anc[[x]]$age_hpd))
  sat_age = unlist(sapply(1:length(aic_species), function(x) sat_anc[[x]]$age_hpd))
  bic_age = unlist(sapply(1:length(aic_species), function(x) bic_anc[[x]]$age_hpd))
  gen_age = unlist(sapply(1:length(aic_species), function(x) gen_anc[[x]]$age_hpd))
  uni_age = unlist(sapply(1:length(aic_species), function(x) uni_anc[[x]]$age_hpd))
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



  # SAMPLE PLOTS

  bg_col = c("gray100", "gray92", "gray84", "gray76", "gray68", "gray60", "gray52", "gray44", "gray36", "gray28", "gray20")
  bg_col = rev(bg_col)
  bg_col = bg_col[ranks]
# 
#   par(mfrow=c(1,1))
# 
#   plot.new()
#   box(lty=1, col="white")
#   abline(v=seq(0.1,0.9,0.1), col=bg_col[1], lwd=300)
#   abline(a=0, b=1, col="black", lty=3)
#   points(sat_hpd[(aic_sat_ix[,2]-num_tips),3], aic_hpd[(aic_sat_ix[,1]-num_tips),3], pch=19, col="gray22", xlab="", ylab="")
#   text(0.08, 1, paste0(as.character(rank_mat[1,1]*100),"%"))
# 
#   plot.new()
#   box(lty=1, col="white")
#   abline(v=seq(0.1,0.9,0.1), col=bg_col[7], lwd=300)
#   abline(a=0, b=1, col="black", lty=3)
#   points(sat_hpd[(sat_uni_ix[,1]-num_tips),3], uni_hpd[(sat_uni_ix[,2]-num_tips),3], pch=19, col="gray22", xlab="", ylab="")
# 
#   plot.new()
#   box(lty=1, col="white")
#   abline(v=seq(0.1,0.9,0.1), col=bg_col[6], lwd=300)
#   abline(a=0, b=1, col="black", lty=3)
#   points(sat_hpd[(sat_gen_ix[,1]-num_tips),3], gen_hpd[(sat_gen_ix[,2]-num_tips),3], pch=19, col="gray22", xlab="", ylab="")
# 


  # * PLOT ALL PAIRWISE COMBINATIONS

  par(mfrow=c(5,5), mar=c(0, 0, 0, 0), oma = c(4, 4, 4, 1))

  # * SATURATED
  plot.new()
  box(lty=1, col="white")
  abline(v=seq(0.1,0.9,0.1), col="white", lwd=300)
  abline(a=0, b=1, col="black", lty=3)
  points(sat_hpd[,3], sat_hpd[,3], pch=19, col="gray22", xlab="", ylab="")
  mtext("saturated", side=2, line=1, adj=0.5, cex=1)

  plot.new()
  box(lty=1, col="white")
  abline(v=seq(0.1,0.9,0.1), col=bg_col[1], lwd=300)
  abline(a=0, b=1, col="black", lty=3)
  text(0.08, 1, paste0(as.character(rank_mat[1,1]*100),"%"))
  points(sat_hpd[(aic_sat_ix[,2]-num_tips),3], aic_hpd[(aic_sat_ix[,1]-num_tips),3], pch=19, col="gray22", xlab="", ylab="")

  plot.new()
  box(lty=1, col="white")
  abline(v=seq(0.1,0.9,0.1), col=bg_col[5], lwd=300)
  abline(a=0, b=1, col="black", lty=3)
  text(0.08, 1, paste0(as.character(rank_mat[5,1]*100),"%"))
  points(sat_hpd[(sat_bic_ix[,1]-num_tips),3], bic_hpd[(sat_bic_ix[,2]-num_tips),3], pch=19, col="gray22", xlab="", ylab="")

  plot.new()
  box(lty=1, col="white")
  abline(v=seq(0.1,0.9,0.1), col=bg_col[6], lwd=300)
  abline(a=0, b=1, col="black", lty=3)
  text(0.08, 1, paste0(as.character(rank_mat[6,1]*100),"%"))
  points(sat_hpd[(sat_gen_ix[,1]-num_tips),3], gen_hpd[(sat_gen_ix[,2]-num_tips),3], pch=19, col="gray22", xlab="", ylab="")

  plot.new()
  box(lty=1, col="white")
  abline(v=seq(0.1,0.9,0.1), col=bg_col[7], lwd=300)
  abline(a=0, b=1, col="black", lty=3)
  text(0.08, 1, paste0(as.character(rank_mat[7,1]*100),"%"))
  points(sat_hpd[(sat_uni_ix[,1]-num_tips),3], uni_hpd[(sat_uni_ix[,2]-num_tips),3], pch=19, col="gray22", xlab="", ylab="")



  # * AIC

  plot.new()
  box(lty=1, col="white")
  abline(v=seq(0.1,0.9,0.1), col=bg_col[1], lwd=300)
  abline(a=0, b=1, col="black", lty=3)
  points(aic_hpd[(aic_sat_ix[,1]-num_tips),3], sat_hpd[(aic_sat_ix[,2]-num_tips),3], pch=19, col="gray22", xlab="", ylab="")
  text(0.08, 1, paste0(as.character(rank_mat[1,1]*100),"%"))
  mtext("aic", side=2, line=1, adj=0.5, cex=1)

  plot.new()
  box(lty=1, col="white")
  abline(v=seq(0.1,0.9,0.1), col="white", lwd=300)
  abline(a=0, b=1, col="black", lty=3)
  points(aic_hpd[,3], aic_hpd[,3], pch=19, col="gray22", xlab="", ylab="")

  plot.new()
  box(lty=1, col="white")
  abline(v=seq(0.1,0.9,0.1), col=bg_col[2], lwd=300)
  abline(a=0, b=1, col="black", lty=3)
  text(0.08, 1, paste0(as.character(rank_mat[2,1]*100),"%"))
  points(aic_hpd[(aic_bic_ix[,1]-num_tips),3], bic_hpd[(aic_bic_ix[,2]-num_tips),3], pch=19, col="gray22", xlab="", ylab="")

  plot.new()
  box(lty=1, col="white")
  abline(v=seq(0.1,0.9,0.1), col=bg_col[3], lwd=300)
  abline(a=0, b=1, col="black", lty=3)
  text(0.08, 1, paste0(as.character(rank_mat[3,1]*100),"%"))
  points(aic_hpd[(aic_gen_ix[,1]-num_tips),3], gen_hpd[(aic_gen_ix[,2]-num_tips),3], pch=19, col="gray22", xlab="", ylab="")

  plot.new()
  box(lty=1, col="white")
  abline(v=seq(0.1,0.9,0.1), col=bg_col[4], lwd=300)
  abline(a=0, b=1, col="black", lty=3)
  text(0.08, 1, paste0(as.character(rank_mat[4,1]*100),"%"))
  points(aic_hpd[(aic_uni_ix[,1]-num_tips),3], uni_hpd[(aic_uni_ix[,2]-num_tips),3], pch=19, col="gray22", xlab="", ylab="")



  # * BIC
  plot.new()
  box(lty=1, col="white")
  abline(v=seq(0.1,0.9,0.1), col=bg_col[5], lwd=300)
  abline(a=0, b=1, col="black", lty=3)
  points(bic_hpd[(sat_bic_ix[,2]-num_tips),3], sat_hpd[(sat_bic_ix[,1]-num_tips),3], pch=19, col="gray22", xlab="", ylab="")
  text(0.08, 1, paste0(as.character(rank_mat[5,1]*100),"%"))
  mtext("bic", side=2, line=1, adj=0.5, cex=1)

  plot.new()
  box(lty=1, col="white")
  abline(v=seq(0.1,0.9,0.1), col=bg_col[2], lwd=300)
  abline(a=0, b=1, col="black", lty=3)
  text(0.08, 1, paste0(as.character(rank_mat[2,1]*100),"%"))
  points(bic_hpd[(aic_bic_ix[,2]-num_tips),3], aic_hpd[(aic_bic_ix[,1]-num_tips),3], pch=19, col="gray22", xlab="", ylab="")

  plot.new()
  box(lty=1, col="white")
  abline(v=seq(0.1,0.9,0.1), col="white", lwd=300)
  abline(a=0, b=1, col="black", lty=3)
  points(bic_hpd[,3], bic_hpd[,3], pch=19, col="gray22", xlab="", ylab="")

  plot.new()
  box(lty=1, col="white")
  abline(v=seq(0.1,0.9,0.1), col=bg_col[8], lwd=300)
  abline(a=0, b=1, col="black", lty=3)
  text(0.08, 1, paste0(as.character(rank_mat[8,1]*100),"%"))
  points(bic_hpd[(bic_gen_ix[,1]-num_tips),3], gen_hpd[(bic_gen_ix[,2]-num_tips),3], pch=19, col="gray22", xlab="", ylab="")


  plot.new()
  box(lty=1, col="white")
  abline(v=seq(0.1,0.9,0.1), col=bg_col[9], lwd=300)
  abline(a=0, b=1, col="black", lty=3)
  text(0.08, 1, paste0(as.character(rank_mat[9,1]*100),"%"))
  points(bic_hpd[(bic_uni_ix[,1]-num_tips),3], uni_hpd[(bic_uni_ix[,2]-num_tips),3], pch=19, col="gray22", xlab="", ylab="")


  # * GENE

  plot.new()
  box(lty=1, col="white")
  abline(v=seq(0.1,0.9,0.1), col=bg_col[6], lwd=300)
  abline(a=0, b=1, col="black", lty=3)
  points(gen_hpd[(sat_gen_ix[,2]-num_tips),3], sat_hpd[(sat_gen_ix[,1]-num_tips),3], pch=19, col="gray22", xlab="", ylab="")
  text(0.08, 1, paste0(as.character(rank_mat[6,1]*100),"%"))
  mtext("gene", side=2, line=1, adj=0.5, cex=1)

  plot.new()
  box(lty=1, col="white")
  abline(v=seq(0.1,0.9,0.1), col=bg_col[3], lwd=300)
  abline(a=0, b=1, col="black", lty=3)
  text(0.08, 1, paste0(as.character(rank_mat[3,1]*100),"%"))
  points(gen_hpd[(aic_gen_ix[,2]-num_tips),3], aic_hpd[(aic_gen_ix[,1]-num_tips),3], pch=19, col="gray22", xlab="", ylab="")

  plot.new()
  box(lty=1, col="white")
  abline(v=seq(0.1,0.9,0.1), col=bg_col[8], lwd=300)
  abline(a=0, b=1, col="black", lty=3)
  text(0.08, 1, paste0(as.character(rank_mat[8,1]*100),"%"))
  points(gen_hpd[(bic_gen_ix[,2]-num_tips),3], bic_hpd[(bic_gen_ix[,1]-num_tips),3], pch=19, col="gray22", xlab="", ylab="")


  plot.new()
  box(lty=1, col="white")
  abline(v=seq(0.1,0.9,0.1), col="white", lwd=300)
  abline(a=0, b=1, col="black", lty=3)
  points(gen_hpd[,3], gen_hpd[,3], pch=19, col="gray22", xlab="", ylab="")


  plot.new()
  box(lty=1, col="white")
  abline(v=seq(0.1,0.9,0.1), col=bg_col[10], lwd=300)
  abline(a=0, b=1, col="black", lty=3)
  text(0.08, 1, paste0(as.character(rank_mat[10,1]*100),"%"))
  points(gen_hpd[(gen_uni_ix[,1]-num_tips),3], uni_hpd[(gen_uni_ix[,2]-num_tips),3], pch=19, col="gray22", xlab="", ylab="")


  # * UNIFORM

  plot.new()
  box(lty=1, col="white")
  abline(v=seq(0.1,0.9,0.1), col=bg_col[7], lwd=300)
  abline(a=0, b=1, col="black", lty=3)
  points(uni_hpd[(sat_uni_ix[,2]-num_tips),3], sat_hpd[(sat_uni_ix[,1]-num_tips),3], pch=19, col="gray22", xlab="", ylab="")
  text(0.08, 1, paste0(as.character(rank_mat[7,1]*100),"%"))
  mtext("uniform", side=2, line=1, adj=0.5, cex=1)
  mtext("saturated", side=1, line=1, adj=0.5, cex=1)


  plot.new()
  box(lty=1, col="white")
  abline(v=seq(0.1,0.9,0.1), col=bg_col[4], lwd=300)
  abline(a=0, b=1, col="black", lty=3)
  points(uni_hpd[(aic_uni_ix[,2]-num_tips),3], aic_hpd[(aic_uni_ix[,1]-num_tips),3], pch=19, col="gray22", xlab="", ylab="")
  text(0.08, 1, paste0(as.character(rank_mat[4,1]*100),"%"))
  mtext("aic", side=1, line=1, adj=0.5, cex=1)


  plot.new()
  box(lty=1, col="white")
  abline(v=seq(0.1,0.9,0.1), col=bg_col[9], lwd=300)
  abline(a=0, b=1, col="black", lty=3)
  points(uni_hpd[(bic_uni_ix[,2]-num_tips),3], bic_hpd[(bic_uni_ix[,1]-num_tips),3], pch=19, col="gray22", xlab="", ylab="")
  text(0.08, 1, paste0(as.character(rank_mat[9,1]*100),"%"))
  mtext("bic", side=1, line=1, adj=0.5, cex=1)


  plot.new()
  box(lty=1, col="white")
  abline(v=seq(0.1,0.9,0.1), col=bg_col[10], lwd=300)
  abline(a=0, b=1, col="black", lty=3)
  points(uni_hpd[(gen_uni_ix[,2]-num_tips),3], gen_hpd[(gen_uni_ix[,1]-num_tips),3], pch=19, col="gray22", xlab="", ylab="")
  text(0.08, 1, paste0(as.character(rank_mat[10,1]*100),"%"))
  mtext("gene", side=1, line=1, adj=0.5, cex=1)

  plot.new()
  box(lty=1, col="white")
  abline(v=seq(0.1,0.9,0.1), col="white", lwd=300)
  abline(a=0, b=1, col="black", lty=3)
  points(uni_hpd[,3], uni_hpd[,3], pch=19, col="gray22", xlab="", ylab="")
  mtext("uniform", side=1, line=1, adj=0.5, cex=1)



  title(main="Divergence Time Comparisons for Shared Nodes\nUnder 5 Partitioning Schemes", sub=study, cex.main=1.5, outer=TRUE)

}
