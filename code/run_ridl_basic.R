library(dplyr)
library(emmeans)
library(fmri.pipeline)
library(readr)

# get trial data
trial_df <- parse_ridl_all("/proj/mnhallqlab/studies/momentum/clpipe/data_onsets")


subject_df <- readRDS("/proj/mnhallqlab/users/michael/fmri.pipeline/inst/example_files/mmclock_subject_data.rds") %>%
  mutate(mr_dir=paste0(mr_dir, "/mni_5mm_aroma")) #make sure we're looking in the right folder

# run_df <- readRDS("/proj/mnhallqlab/users/michael/fmri.pipeline/inst/example_files/mmclock_run_data.rds")

meg_ranefs <- readRDS("/proj/mnhallqlab/projects/clock_analysis/meg/data/MEG_betas_wide_echange_vmax_reward_Nov30_2021.RDS") %>%
  dplyr::select(-contains("corr")) %>%
  mutate(id = as.integer(id)) %>%
  mutate(across(.cols = starts_with("avg_"), .fns = ~ as.vector(scale(.x)))) # z-score betas

gpa <- setup_glm_pipeline(analysis_name="mmclock_nov2021", scheduler="slurm",
  output_directory = "/proj/mnhallqlab/users/michael/mmclock_entropy",
  subject_data=subject_df, trial_data=trial_df, # run_data=run_df,
  tr=1.0, drop_volumes = 2,
  n_expected_runs=8,
  l1_models=NULL, l2_models=NULL, l3_models=NULL,
  fmri_file_regex="nfaswuktm_clock[1-8]_5\\.nii\\.gz",
  fmri_path_regex="clock[0-9]",
  run_number_regex=".*clock([0-9]+)_5.*",
  confound_settings=list(
    motion_params_file = "motion.par",
    confound_input_file="nuisance_regressors.txt",
    confound_input_colnames = c("csf", "dcsf", "wm", "dwm"),
    l1_confound_regressors = c("csf", "dcsf", "wm", "dwm"),
    exclude_run = "max(FD) > 5 | sum(FD > .9)/length(FD) > .10", #this must evaluate to a scalar per run
    exclude_subject = "n_good_runs < 4",
    truncate_run = "(FD > 0.9 & time > last_offset) | (time > last_offset + last_isi)",
    spike_volumes = NULL
  ),
  parallel=list(
    fsl=list(l1_feat_alljobs_time="72:00:00") # give all l1 runs a while to clear the scheduler
  )
)

# gpa <- build_l1_models(gpa, from_spec_file = "entropy_all.yaml")

gpa <- build_l2_models(gpa)
gpa <- build_l3_models(gpa)

gpa <- readRDS("mmclock_nov2021_addmegbetas.rds")

gpa$subject_data <- gpa$subject_data %>% left_join(meg_ranefs, by = "id")

gpa <- build_l3_models(gpa)
saveRDS(gpa, file = "mmclock_nov2021_addmegbetas.rds")

run_glm_pipeline(gpa)


