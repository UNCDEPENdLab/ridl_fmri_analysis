library(dplyr)
library(emmeans)
library(fmri.pipeline)
library(readr)
library(glue)

# get trial data
setwd("/proj/mnhallqlab/projects/ridl_fmri_analysis/code")
trial_df <- read.csv("../output/ridl_combined_fmri.csv.gz") %>%
  mutate(
    id = as.character(id),
    # recode outcome_fac to note the curtain outcome in accept/reject so that feedback trials can be separated by type
    outcome_fac = if_else(trial_type == "accrej", "curtain", outcome_fac)
  )

run_df <- generate_run_data_from_bids("/proj/mnhallqlab/studies/momentum/clpipe/data_postproc2", task_name = "ridl", suffix = "_postproccessed")

# xtabs(~id, run_df)

subject_df <- run_df %>%
  group_by(id) %>%
  filter(row_number() == 1) %>%
  mutate(mr_dir = dirname(run_nifti)) %>%
  select(id, mr_dir) %>%
  ungroup()

# run_df <- readRDS("/proj/mnhallqlab/users/michael/fmri.pipeline/inst/example_files/mmclock_run_data.rds")

gpa <- setup_glm_pipeline(analysis_name="basic_apr2022", scheduler="slurm",
  output_directory = "/proj/mnhallqlab/users/michael/ridl_fmri",
  trial_data=trial_df, subject_data = subject_df, run_data=run_df,
  tr=.635, drop_volumes = 2,
  l1_models=NULL, l2_models=NULL, l3_models=NULL,
  n_expected_runs=4,
  confound_settings=list(
    confound_input_colnames = c("csf", "csf_derivative1", "white_matter", "white_matter_derivative1"), # assumption
    l1_confound_regressors = c("csf", "csf_derivative1", "white_matter", "white_matter_derivative1"),
    na_strings=c("NA", "n/a"), #weird fmriprep-ism for confound files
    #exclude_run = "max(framewise_displacement) > 5 | sum(framewise_displacement > .5)/length(framewise_displacement) > .15", #this must evaluate to a scalar per run
    #exclude_subject = "n_good_runs < 2",
    exclude_run = NULL,
    exclude_subject = NULL,
    #truncate_run = "framewise_displacement > 2 & (time > last_offset + )" # 2 seconds after last offset
    truncate_run = "(framewise_displacement > 0.9 & time > last_offset) | (time > last_offset + last_isi)",
    spike_volumes = NULL
  ),
  parallel=list(
    fsl=list(l1_feat_alljobs_time="72:00:00") # give all l1 runs a while to clear the scheduler
  )
)

#l1_models <- gpa$l1_models

gpa <- build_l1_models(gpa, from_spec_file = "ridl_basic_l1.yaml")
gpa <- build_l2_models(gpa)
gpa <- build_l3_models(gpa)
#gpa$parallel$compute_environment <- c("module use /proj/mnhallqlab/sw/modules", "module load r/4.1.2_depend")
saveRDS(gpa, file = "/proj/mnhallqlab/users/michael/ridl_fmri/gpa_snapshot_20apr2022.rds")

run_glm_pipeline(gpa)

gpa <- readRDS(file = "/proj/mnhallqlab/users/michael/ridl_fmri/basic_apr2022/basic_apr2022.rds")
# for combining results after run completes
res <- combine_feat_l3_to_afni(
  gpa, 
  feat_l3_combined_filename = "{gpa$output_directory}/feat_l3_combined/L1m-{l1_model}/L2m-{l2_model}_L3m-{l3_model}/l2c-{l2_cope_name}_l3c-{l3_cope_name}_stats",
  feat_l3_combined_briknames = "l1c-{l1_cope_name}",
  template_brain = "/proj/mnhallqlab/lab_resources/standard/mni_icbm152_nlin_asym_09c/mni_icbm152_t1_tal_nlin_asym_09c_brain.nii"
)
