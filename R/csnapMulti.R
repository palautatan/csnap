#' csnapMulti()
#'
#' Create a Compare Shared Node Ages Plot (CSNAP) between each pairwise sharing (AIC and BIC share these nodes / Gene and AIC share different nodes)
#' @param tre_files_1 A vector of 5 .tre files in comparison order
#' @param tre_files_2 Another vector of 5 .tre files in comparison order
#' @keywords ancestry
#' @export

csnapMulti = function(tre_files_1, tre_files_2) {
  
  # * 1A. READ TRE_FILES_1
  anc_1a   = createAncestry(tre_files_1[1])
  anc_2a   = createAncestry(tre_files_1[2])
  anc_3a   = createAncestry(tre_files_1[3])
  anc_4a   = createAncestry(tre_files_1[4])
  anc_5a   = createAncestry(tre_files_1[5])
  
  # * 1A. READ TRE_FILES_2
  anc_1b   = createAncestry(tre_files_2[1])
  anc_2b   = createAncestry(tre_files_2[2])
  anc_3b   = createAncestry(tre_files_2[3])
  anc_4b   = createAncestry(tre_files_2[4])
  anc_5b   = createAncestry(tre_files_2[5])
  
  
  
  # * 2A. CREATE LIST OF SPECIES FOR TRE_FILES_1 & SORT
  anc_1a_species = lapply(sapply(1:length(anc_1a), function(x) anc_1a[[x]]$species_under), sort)
  anc_2a_species = lapply(sapply(1:length(anc_2a), function(x) anc_2a[[x]]$species_under), sort)
  anc_3a_species = lapply(sapply(1:length(anc_3a), function(x) anc_3a[[x]]$species_under), sort)
  anc_4a_species = lapply(sapply(1:length(anc_4a), function(x) anc_4a[[x]]$species_under), sort)
  anc_5a_species = lapply(sapply(1:length(anc_5a), function(x) anc_5a[[x]]$species_under), sort)
  
  
  # * 2B. CREATE LIST OF SPECIES FOR TRE_FILES_2 & SORT
  anc_1b_species = lapply(sapply(1:length(anc_1b), function(x) anc_1b[[x]]$species_under), sort)
  anc_2b_species = lapply(sapply(1:length(anc_2b), function(x) anc_2b[[x]]$species_under), sort)
  anc_3b_species = lapply(sapply(1:length(anc_3b), function(x) anc_3b[[x]]$species_under), sort)
  anc_4b_species = lapply(sapply(1:length(anc_4b), function(x) anc_4b[[x]]$species_under), sort)
  anc_5b_species = lapply(sapply(1:length(anc_5b), function(x) anc_5b[[x]]$species_under), sort)
  
  num_tips = length(anc_1a_species[sapply(anc_1a_species, is.null)])
  total_nodes = length(anc_1a_species[(num_tips+1):length(anc_1a_species)])
  
  
  
  
  # * 3. INTERSECT TO GET SHARED NODES
  
  # * ROW 1
  anc1a_anc1b = intersect(anc_1a_species, anc_1b_species)
  anc1a_anc2b = intersect(anc_1a_species, anc_2b_species)
  anc1a_anc3b = intersect(anc_1a_species, anc_3b_species)
  anc1a_anc4b = intersect(anc_1a_species, anc_4b_species)
  anc1a_anc5b = intersect(anc_1a_species, anc_5b_species)
  
  
  # * ROW 2
  anc2a_anc1b = intersect(anc_2a_species, anc_1b_species)
  anc2a_anc2b = intersect(anc_2a_species, anc_2b_species)
  anc2a_anc3b = intersect(anc_2a_species, anc_3b_species)
  anc2a_anc4b = intersect(anc_2a_species, anc_4b_species)
  anc2a_anc5b = intersect(anc_2a_species, anc_5b_species)
  
  
  # * ROW 3
  anc3a_anc1b = intersect(anc_3a_species, anc_1b_species)
  anc3a_anc2b = intersect(anc_3a_species, anc_2b_species)
  anc3a_anc3b = intersect(anc_3a_species, anc_3b_species)
  anc3a_anc4b = intersect(anc_3a_species, anc_4b_species)
  anc3a_anc5b = intersect(anc_3a_species, anc_5b_species)  
  
  
  # * ROW 4
  anc4a_anc1b = intersect(anc_4a_species, anc_1b_species)
  anc4a_anc2b = intersect(anc_4a_species, anc_2b_species)
  anc4a_anc3b = intersect(anc_4a_species, anc_3b_species)
  anc4a_anc4b = intersect(anc_4a_species, anc_4b_species)
  anc4a_anc5b = intersect(anc_4a_species, anc_5b_species)  
  
  # * ROW 5
  anc5a_anc1b = intersect(anc_5a_species, anc_1b_species)
  anc5a_anc2b = intersect(anc_5a_species, anc_2b_species)
  anc5a_anc3b = intersect(anc_5a_species, anc_3b_species)
  anc5a_anc4b = intersect(anc_5a_species, anc_4b_species)
  anc5a_anc5b = intersect(anc_5a_species, anc_5b_species)  
  
  
  intersections = list(anc1a_anc1b, anc1a_anc2b, anc1a_anc3b, anc1a_anc4b, anc1a_anc5b, 
                       anc2a_anc1b, anc2a_anc2b, anc2a_anc3b, anc2a_anc4b, anc2a_anc5b, 
                       anc3a_anc1b, anc3a_anc2b, anc3a_anc3b, anc3a_anc4b, anc3a_anc5b, 
                       anc4a_anc1b, anc4a_anc2b, anc4a_anc3b, anc4a_anc4b, anc4a_anc5b, 
                       anc5a_anc1b, anc5a_anc2b, anc5a_anc3b, anc5a_anc4b, anc5a_anc5b)
  
  intersections = lapply(intersections, function(x) x[-sapply(x, is.null)])
  
  
  
  # * 4. FRACTION SHARED
  shared_fractions = c()
  
  for (intersection in intersections) {
    length   = length(intersection)
    fraction = round(length / total_nodes, 4)
    shared_fractions = c(shared_fractions, fraction)
  }
  
  # * RANK BY FRACTION
  ranks = round(round((((shared_fractions/max(shared_fractions)))*10),2))
  rank_mat = cbind(shared_fractions, ranks)
  
  
  
  # * 5. MATCH UP INTERSECTS
  
  # * ROW 1
  anc1a_anc1b_ix = matchNodes2(anc_1a_species, anc_1b_species, anc1a_anc1b, num_tips)
  anc1a_anc2b_ix = matchNodes2(anc_1a_species, anc_2b_species, anc1a_anc2b, num_tips)
  anc1a_anc3b_ix = matchNodes2(anc_1a_species, anc_3b_species, anc1a_anc3b, num_tips)
  anc1a_anc4b_ix = matchNodes2(anc_1a_species, anc_4b_species, anc1a_anc4b, num_tips)
  anc1a_anc5b_ix = matchNodes2(anc_1a_species, anc_5b_species, anc1a_anc5b, num_tips)
  
  
  # * ROW 2
  anc2a_anc1b_ix = matchNodes2(anc_2a_species, anc_1b_species, anc2a_anc1b, num_tips)
  anc2a_anc2b_ix = matchNodes2(anc_2a_species, anc_2b_species, anc2a_anc2b, num_tips)
  anc2a_anc3b_ix = matchNodes2(anc_2a_species, anc_3b_species, anc2a_anc3b, num_tips)
  anc2a_anc4b_ix = matchNodes2(anc_2a_species, anc_4b_species, anc2a_anc4b, num_tips)
  anc2a_anc5b_ix = matchNodes2(anc_2a_species, anc_5b_species, anc2a_anc5b, num_tips)
  
  
  # * ROW 3
  anc3a_anc1b_ix = matchNodes2(anc_3a_species, anc_1b_species, anc1a_anc1b, num_tips)
  anc3a_anc2b_ix = matchNodes2(anc_3a_species, anc_2b_species, anc1a_anc2b, num_tips)
  anc3a_anc3b_ix = matchNodes2(anc_3a_species, anc_3b_species, anc1a_anc3b, num_tips)
  anc3a_anc4b_ix = matchNodes2(anc_3a_species, anc_4b_species, anc1a_anc4b, num_tips)
  anc3a_anc5b_ix = matchNodes2(anc_3a_species, anc_5b_species, anc1a_anc5b, num_tips)
  
  
  # * ROW 4
  anc4a_anc1b_ix = matchNodes2(anc_4a_species, anc_1b_species, anc4a_anc1b, num_tips)
  anc4a_anc2b_ix = matchNodes2(anc_4a_species, anc_2b_species, anc4a_anc2b, num_tips)
  anc4a_anc3b_ix = matchNodes2(anc_4a_species, anc_3b_species, anc4a_anc3b, num_tips)
  anc4a_anc4b_ix = matchNodes2(anc_4a_species, anc_4b_species, anc4a_anc4b, num_tips)
  anc4a_anc5b_ix = matchNodes2(anc_4a_species, anc_5b_species, anc4a_anc5b, num_tips)
  
  # * ROW 5
  anc5a_anc1b_ix = matchNodes2(anc_5a_species, anc_1b_species, anc5a_anc1b, num_tips)
  anc5a_anc2b_ix = matchNodes2(anc_5a_species, anc_2b_species, anc5a_anc2b, num_tips)
  anc5a_anc3b_ix = matchNodes2(anc_5a_species, anc_3b_species, anc5a_anc3b, num_tips)
  anc5a_anc4b_ix = matchNodes2(anc_5a_species, anc_4b_species, anc5a_anc4b, num_tips)
  anc5a_anc5b_ix = matchNodes2(anc_5a_species, anc_5b_species, anc5a_anc5b, num_tips)
  
  
  
  # * 6. APPEND COMMENTS
  anc_1a   = addComments(tre_files_1[1], anc_1a)
  anc_2a   = addComments(tre_files_1[2], anc_2a)
  anc_3a   = addComments(tre_files_1[3], anc_3a)
  anc_4a   = addComments(tre_files_1[4], anc_4a)
  anc_5a   = addComments(tre_files_1[5], anc_5a)
  
  anc_1b   = addComments(tre_files_1[1], anc_1b)
  anc_2b   = addComments(tre_files_1[2], anc_2b)
  anc_3b   = addComments(tre_files_1[3], anc_3b)
  anc_4b   = addComments(tre_files_1[4], anc_4b)
  anc_5b   = addComments(tre_files_1[5], anc_5b)
  
  
  # * 7. AGE HPD'S FOR ALL NODES
  anc1a_age = unlist(sapply(1:length(anc_1a_species), function(x) anc_1a[[x]]$age_hpd))
  anc2a_age = unlist(sapply(1:length(anc_2a_species), function(x) anc_2a[[x]]$age_hpd))
  anc3a_age = unlist(sapply(1:length(anc_3a_species), function(x) anc_3a[[x]]$age_hpd))
  anc4a_age = unlist(sapply(1:length(anc_4a_species), function(x) anc_4a[[x]]$age_hpd))
  anc5a_age = unlist(sapply(1:length(anc_5a_species), function(x) anc_5a[[x]]$age_hpd))
  
  anc1b_age = unlist(sapply(1:length(anc_1b_species), function(x) anc_1b[[x]]$age_hpd))
  anc2b_age = unlist(sapply(1:length(anc_2b_species), function(x) anc_2b[[x]]$age_hpd))
  anc3b_age = unlist(sapply(1:length(anc_3b_species), function(x) anc_3b[[x]]$age_hpd))
  anc4b_age = unlist(sapply(1:length(anc_4b_species), function(x) anc_4b[[x]]$age_hpd))
  anc5b_age = unlist(sapply(1:length(anc_5b_species), function(x) anc_5b[[x]]$age_hpd))  

  age_df  = t(data.frame(anc1a_age, anc2a_age, anc3a_age, anc4a_age, anc5a_age,
                         anc1b_age, anc2b_age, anc3b_age, anc4b_age, anc5b_age))
  
  
  
  # * 8. CREATE PLOTTABLE MEAN VECTORS
  
  # * ROW 1
  anc1a_hpd = matrix(as.numeric(unlist(strsplit(age_df[1,], ","))), ncol=2, byrow=TRUE)
  anc1a_hpd = cbind(anc1a_hpd, (anc1a_hpd[,1] + anc1a_hpd[,2])/2)
  
  anc2a_hpd = matrix(as.numeric(unlist(strsplit(age_df[2,], ","))), ncol=2, byrow=TRUE)
  anc2a_hpd = cbind(anc2a_hpd, (anc2a_hpd[,1] + anc2a_hpd[,2])/2)
  
  anc3a_hpd = matrix(as.numeric(unlist(strsplit(age_df[3,], ","))), ncol=2, byrow=TRUE)
  anc3a_hpd = cbind(anc3a_hpd, (anc3a_hpd[,1] + anc3a_hpd[,2])/2)
  
  anc4a_hpd = matrix(as.numeric(unlist(strsplit(age_df[4,], ","))), ncol=2, byrow=TRUE)
  anc4a_hpd = cbind(anc4a_hpd, (anc4a_hpd[,1] + anc4a_hpd[,2])/2)
  
  anc5a_hpd = matrix(as.numeric(unlist(strsplit(age_df[5,], ","))), ncol=2, byrow=TRUE)
  anc5a_hpd = cbind(anc5a_hpd, (anc5a_hpd[,1] + anc5a_hpd[,2])/2)
  
  # * ROW 2
  anc1b_hpd = matrix(as.numeric(unlist(strsplit(age_df[6,], ","))), ncol=2, byrow=TRUE)
  anc1b_hpd = cbind(anc1b_hpd, (anc1b_hpd[,1] + anc1b_hpd[,2])/2)
  
  anc2b_hpd = matrix(as.numeric(unlist(strsplit(age_df[7,], ","))), ncol=2, byrow=TRUE)
  anc2b_hpd = cbind(anc2b_hpd, (anc2b_hpd[,1] + anc2b_hpd[,2])/2)
  
  anc3b_hpd = matrix(as.numeric(unlist(strsplit(age_df[8,], ","))), ncol=2, byrow=TRUE)
  anc3b_hpd = cbind(anc3b_hpd, (anc3b_hpd[,1] + anc3b_hpd[,2])/2)
  
  anc4b_hpd = matrix(as.numeric(unlist(strsplit(age_df[9,], ","))), ncol=2, byrow=TRUE)
  anc4b_hpd = cbind(anc4b_hpd, (anc4b_hpd[,1] + anc4b_hpd[,2])/2)
  
  anc5b_hpd = matrix(as.numeric(unlist(strsplit(age_df[10,], ","))), ncol=2, byrow=TRUE)
  anc5b_hpd = cbind(anc5b_hpd, (anc5b_hpd[,1] + anc5b_hpd[,2])/2)
