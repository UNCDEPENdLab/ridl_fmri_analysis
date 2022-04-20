parse_ridl <- function(sub_id = NULL, ridl_dir, force=FALSE, matlab_dir = "/nas/longleaf/apps/matlab/2021b/bin") {
  pacman::p_load(checkmate, data.table, dplyr, matlabr, tidyr, glue)
  checkmate::assert_directory_exists(ridl_dir)
  if (is.null(sub_id)) {
    # assume that ridl_dir is the full path to the subject folder
    sub_dir <- ridl_dir
  } else {
    sub_dir <- file.path(ridl_dir, sub_id)
  }
  
  message(glue("Working on subject directory: {sub_dir}"))
  
  checkmate::assert_directory_exists(sub_dir)
  
  if (!file.exists(file.path(sub_dir, "log_file.mat"))) {
    odirs <- list.dirs(path = sub_dir, recursive = FALSE)[1]
    message(glue("Using first subdirectory in folder: {odirs}"))
    sub_dir <- odirs
  }
  
  expected_file <- Sys.glob(file.path(sub_dir, "*_timing.csv.gz"))[1L]
  if (file.exists(expected_file) && isFALSE(force)) {
    message(glue("Parsed file already exists: {expected_file}"))
    return(data.table::fread(expected_file, data.table = FALSE))
  }

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
  
  # note that in recent versions after Oct 2021, feedback_lag (isi duration) was removed from log.trial_information. For
  # consistency across datasets and checks on isi durations, add this back in if it is missing.
  run_matlab_code(
    glue(
      "
        load '<log_file>';
        log.trial_information{:,'id'} = {log.subjectID}; % add subject id as a field
        log.trial_information{:,'experiment_date'} = {log.Experiment_Start_Time}; % add experiment date as field
        log.trial_information{:,'experiment_build'} = {log.parameters_used.experiment_build}; % add experiment date as field
        if ~ismember('feedback_lag',log.trial_information.Properties.VariableNames)
          log.trial_information.row_order = (1:size(log.trial_information,1))';
          spin_data = reshape(log.parameters_used.subject.spins',[],1);
          spin_table = table(spin_data, repmat(0:83, 1, 4)', reshape(repmat(1000:1003, 84, 1), [], 1), 'VariableNames', {'feedback_lag', 'trial', 'block'});
          aa = outerjoin(log.trial_information, spin_table, 'Type','left','MergeKeys',true);
          aa = sortrows(aa,{'row_order'}); % restore original row ordering (join resorts data)
          aa.row_order=[]; % remove sort variable
          log.trial_information = aa;
        end

        writetable(log.data.trialdata, '<trial_data_fname>');
        writetable(log.trial_information, '<trial_info_fname>');
        writetable(struct2table(structfun(@transpose,log.accept_reject,'UniformOutput',false)), '<accrej_info_fname>');
        ", .open="<", .close=">"     
    ),
    endlines = FALSE
  )
  
#"if ~ismember('feedback_lag',log.trial_information.Properties.VariableNames), log.trial_information.feedback_lag = reshape(log.parameters_used.subject.spins',[],1); end",

  trial_data <- data.table::fread(trial_data_fname) %>% select(-ends_with("_time")) # these are timestamps in ms code (not human readable)
  trial_info <- data.table::fread(trial_info_fname)
  if (is.null(trial_info$times_too_slow)) {
    # this will also occur for new 
    message("No times_too_slow in trial_info data. Assuming that there were no slow trials and adding zeros for this column.")
    trial_info$times_too_slow <- 0
  }
  
  # convert experiment time to date/time using anytime package. Note that we need "20" to be prepended to year to make "2021"
  trial_info$experiment_date <- anytime::anytime(paste0("20", trial_info$experiment_date))
  
  # for a few subjects (early task version?), RT was not recorded
  if (!"reaction_time" %in% names(trial_info)) { 
    message(glue("Filling in reaction_time based on choice_time - stim_time for subject directory: {sub_dir}"))
    trial_info$reaction_time <- trial_info$choice_time - trial_info$stim_time
  }
  
  trial_df <- trial_data %>%
    merge(trial_info, by = c("block", "trial")) %>%
    dplyr::select(-trial_duration) %>% # , -trial_started) %>% # not currently used
    mutate(
      site = if_else(id >= 540000 & id < 600000, 'unc', 'pitt'), #simple, but effective?
      trial = trial + 1, # switch to one-based indexing
      # type = if_else(stim2 == -1000, "accrej", "newlearn")
      trial_type = if_else(feedback == 0, "accrej", "new_learn"), # human-readable trial code
      
      # if choice_time is NA, reaction_time should be, too (currently it shows 0.0)
      reaction_time = if_else(is.na(choice_time), NA_real_, reaction_time) 
    ) %>%
    dplyr::rename(
      iti_ideal_lag = isi_length,  # ITI in seconds (starts trial in current scheme)
      feedback_isi = feedback_lag # ISI between choice and feedback (wheel spinning)
    ) %>%
    group_by(block) %>%
    mutate(
      run_number = cur_group_id(),
      trial_started_next = dplyr::lead(trial_started, 1, order_by = trial),
      stim_time_next = dplyr::lead(stim_time, 1, order_by = trial),
      iti_ideal = dplyr::lead(iti_ideal_lag, 1, order_by = trial), # this puts iti_ideal as an event after the trial
      too_slow_lead = dplyr::lead(times_too_slow, 1, order_by = trial) # are there repeats on the next trial (bumps out timing)
    ) %>%
    ungroup() %>%
    # empty trials have all 0s on the row that are prepopulated
    mutate(
      bad_trial=rowSums(dplyr::select(., trial_started, stim_time, choice_time, feedback_time)) < 1e-3,
      reaction_time=if_else(bad_trial==TRUE, NA_real_, reaction_time),
      outcome_fac=factor(outcome, levels=c(-1,0,1), labels=c("loss", "nothing", "win"))
    ) %>% # set reaction_time to NA on bad trials
    dplyr::select(id, run_number, trial, block, everything())

  apply_feedback_correction <- FALSE
  if (trial_df$site[1] == "unc" && trial_df$experiment_date[1] < as.POSIXct('2022-03-16')) {
    message("Subtracting 1.0s from feedback_time for new_learn trials for early UNC subject")
    apply_feedback_correction <- TRUE
  } else if (trial_df$site[1] == "pitt" && trial_df$experiment_date[1] < as.POSIXct('2022-03-11') && trial_df$id[1] != 221604) {
    # note that 221604 seems to have the correct feedback times for new learn trials. The subtraction correction leads to 
    # negative times in calculated versus recorded feedback
    message("Subtracting 1.0s from feedback_time for new_learn trials for early Pitt subject")
    apply_feedback_correction <- TRUE
  }
  
  if (isTRUE(apply_feedback_correction)) {
    trial_df <- trial_df %>%
      mutate(
        # for new learn trials, feedback_time is accidentally offset time
        feedback_time = if_else(trial_type == "new_learn", feedback_time - 1.0, feedback_time)
      )  
  }

  # calculate observed ITI after timing adjustment is made to feedback_time
  trial_df <- trial_df %>%
    mutate(iti_actual = stim_time_next - (feedback_time + 1.0))
  
  expected_file <- file.path(sub_dir, paste0(trial_df$id[1L], "_timing.csv.gz"))
  data.table::fwrite(trial_df, file = expected_file)
  return(trial_df)
}
