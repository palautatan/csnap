# tre_files_path = "/Volumes/treehouse3/Dropbox/Effect_of_partitioning/step_3_forward/output/RevBayes/Kaffenberger_2011"
# csnapPairs("/Volumes/treehouse3/Dropbox/Effect_of_partitioning/step_3_forward/output/RevBayes/Kaffenberger_2011")

folders = list.files("/Volumes/treehouse3/Dropbox/backup/2018-04-02/csnap_2/RevBayes", full.names=TRUE)

setwd("/Volumes/treehouse3/Dropbox/backup/2018-04-02/csnap_2")

for (each_folder in folders[4:length(folders)]) {
  study = strsplit(each_folder, "\\/")[[1]][length(strsplit(each_folder, "\\/")[[1]])]
  file_name = paste0("examples/",study,"_pairs.png")
  png(file_name, width=900, height=900)
  csnapPairs(each_folder)
  dev.off()
}
