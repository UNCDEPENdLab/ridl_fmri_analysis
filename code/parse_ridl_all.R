# function to read RiDL behav from multiple files
parse_ridl_all <- function(subj_folder) {
  require(glue)
  checkmate::assert_directory_exists(subj_folder)
  folder_list <- list.dirs(path = subj_folder, recursive = FALSE)
  folder_list <- sapply(folder_list, function(x) {
    if (!file.exists(file.path(x, "log_file.mat"))) {
      subdir <- list.dirs(path = x, recursive = FALSE)[1]
      message(glue("Using first subdirectory in folder: {subdir}"))
    } else {
      return(x)
    }
  })
  matlab_dir <- "/nas/longleaf/apps/matlab/2021b/bin"
  trial_df <- lapply(folder_list, function(ff) {
    parse_ridl(ridl_dir = ff, matlab_dir = matlab_dir)
  })
  browser()
}