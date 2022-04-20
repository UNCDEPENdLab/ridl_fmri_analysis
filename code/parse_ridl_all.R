# function to read RiDL behav from multiple files
parse_ridl_all <- function(subj_folder, force=FALSE, matlab_dir = "/nas/longleaf/apps/matlab/2021b/bin") {
  require(glue)
  checkmate::assert_directory_exists(subj_folder)
  folder_list <- list.dirs(path = subj_folder, recursive = FALSE)
  
  trial_df <- parallel::mclapply(folder_list, function(ff) {
    parse_ridl(ridl_dir = ff, force=force, matlab_dir = matlab_dir)
  }, mc.cores = 8)
  
  return(dplyr::bind_rows(trial_df))
}