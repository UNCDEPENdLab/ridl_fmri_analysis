parse_ridl <- function(sub_id, ridl_dir, matlab_dir = "/nas/longleaf/apps/matlab/2021b/bin") {
  pacman::p_load(checkmate, data.table, dplyr, matlabr, tidyr, glue)
  checkmate::assert_directory_exists(ridl_dir)
  sub_dir <- file.path(ridl_dir, sub_id)
  checkmate::assert_directory_exists(sub_dir)

  # four files in each dir -- only need log at present
  log_file <- file.path(sub_dir, "log_file.mat")
  if (!checkmate::test_file_exists(log_file)) {
    stop(glue("Cannot find log file: {log_file}"))
  }

  trial_data_fname <- tempfile(fileext = ".csv")
  trial_info_fname <- tempfile(fileext = ".csv")
  accrej_info_fname <- tempfile(fileext = ".csv")

  if (!is.null(matlab_dir)) {
    cat("Setting matlab.path\n")
    checkmate::assert_directory_exists(matlab_dir)
    options(matlab.path = matlab_dir)
  }
  run_matlab_code(
    c(
      glue("load '{log_file}';"),
      "log.trial_information{:,'id'} = {log.subjectID};", # add subject id as a field
      glue("writetable(log.data.trialdata, '{trial_data_fname}');"),
      glue("writetable(log.trial_information, '{trial_info_fname}');"),
      glue("writetable(struct2table(structfun(@transpose,log.accept_reject,'UniformOutput',false)), '{accrej_info_fname}');") # not using currently
    ),
    endlines = FALSE
  )

  trial_data <- data.table::fread(trial_data_fname) %>% select(-ends_with("_time")) # these are timestamps in ms code (not human readable)
  trial_info <- data.table::fread(trial_info_fname)
  trial_df <- trial_data %>%
    merge(trial_info, by = c("block", "trial")) %>%
    dplyr::select(-trial_duration) %>% # , -trial_started) %>% # not currently used
    mutate(
      trial = trial + 1, # switch to one-based indexing
      # type = if_else(stim2 == -1000, "accrej", "newlearn")
      trial_type = if_else(feedback == 0, "accrej", "new_learn") # human-readable trial code
    ) %>%
    dplyr::rename(iti = isi_length) %>% # ITI in seconds (starts trial in current scheme)
    group_by(block) %>%
    mutate(
      run_number = cur_group_id(),
      next_start = dplyr::lead(trial_started)
    ) %>%
    ungroup() %>%
    dplyr::select(id, run_number, trial, block, everything())

  # where this is going -- in progress for conversion to long events format
  #    pivot_longer(cols = c(stim_time, choice_time, feedback_time), names_to = "event_type", values_to = "onset") %>%

  return(trial_df)
}

ridl_dir <- "/proj/mnhallqlab/users/michael/ridl_fmri_analysis/data/momentum"
trial_df <- parse_ridl(sub_id = "220256", ridl_dir = ridl_dir)

#txt <- data.table::fread("/proj/mnhallqlab/studies/momentum/clpipe/data_onsets/220256/Subject_220256.txt")

# blows up
#x <- read.mat("/proj/mnhallqlab/studies/momentum/clpipe/data_onsets/220256/log_file.mat")

# checks
calc_rt <- trial_df$choice_time - trial_df$stim_time
rt_diff <- trial_df$reaction_time - calc_rt # looks reasonable
summary(rt_diff)
calc_feedback <- trial_df$choice_time + trial_df$feedback_lag
feedback_diff <- trial_df$feedback_time - calc_feedback

summary(feedback_diff)

# I would have guessed that if feedback_time is the onset of trial feedback, then the next start should be around 1.0s, with maybe a bit of an extension for 
onset_diff <- trial_df$next_start - trial_df$feedback_time

summary(onset_diff)

# about 40% of trials do have a ~1.0s difference
prop.table(table(abs(onset_diff - 1) < .1))

# 40% of trials have ~ 0s difference
prop.table(table(onset_diff < .1))

# 17% have differences > 1.2s
prop.table(table(onset_diff > 1.2))

# there is some dependency in the offsets -- when onset_diff is 1.0, feedback_diff is often 0.0 and vice versa
cbind(feedback_diff, onset_diff) %>% head(n=20)

fwrite(trial_df, file="../output/example_trial_data.csv")