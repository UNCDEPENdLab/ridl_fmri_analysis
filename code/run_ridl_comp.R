library(dplyr)
library(emmeans)
library(fmri.pipeline)
library(readr)
library(glue)

# get trial data
setwd("/proj/mnhallqlab/users/michael/ridl_fmri_analysis/code")
trial_df <- read.csv("../output/ridl_combined_fmri.csv.gz") %>%
  mutate(
    id = as.character(id),
    # recode outcome_fac to note the curtain outcome in accept/reject so that feedback trials can be separated by type
    outcome_fac = if_else(trial_type == "accrej", "curtain", outcome_fac)
  ) %>%
  filter(id != 440443) # time = 0 problems

# make some symlinks to make BIDS-compatible
# sdirs <- grep("sub-.*", list.dirs("/proj/mnhallqlab/no_backup/momentum/MR_Proc", recursive = FALSE), value = TRUE)
# for (ss in sdirs) {
#   rdirs <- grep("ridl_run-\\d+", list.dirs(file.path(ss, "func"), recursive = FALSE), value = TRUE)
#   for (rr in rdirs) {
#     expect_nif <- Sys.glob(paste0(rr, "/nfas*_sub-*_6.nii.gz"))
#     if (length(expect_nif) == 1L) {
#       file.symlink(expect_nif, file.path(ss, "func", sub("nfaswuktm_", "", basename(expect_nif))))
#     }

#     expect_confound <- file.path(rr, "nuisance_regressors.txt")
#     if (file.exists(expect_confound)) {
#       file.symlink(expect_confound, file.path(ss, "func", sub("nfaswuktm_(sub-.*)_bold_6.nii.gz", "\\1-confounds.tsv", basename(expect_nif))))
#     }
#   }
# }

run_df <- generate_run_data_from_bids("/proj/mnhallqlab/no_backup/momentum/MR_Proc", task_name = "ridl", suffix = "bold_6")

run_df <- run_df %>%
  filter(!id %in% c(
    32765, # pilot
    110392, # pilot
    440300, # QC clipping failure
    220256, # raw_QA: axial slice head is significantly tilted
    540042, # Raw_QA: Severely blurred,
    440453 # for now, a bunch of choice onsets == 0... need to sort out!
  )) %>%
  mutate(
    motion_params_file = glue::glue_data(., "{mr_dir}/func/ridl_run-{sprintf('%02d', run_number)}/motion.par")
  )

table(file.exists(run_df$motion_params_file))

# xtabs(~id, run_df)

subject_df <- run_df %>%
  group_by(id) %>%
  filter(row_number() == 1) %>%
  mutate(mr_dir = dirname(run_nifti)) %>%
  select(id, mr_dir) %>%
  ungroup()

selfreport_df <- read.csv("/proj/mnhallqlab/projects/ridl_fmri_analysis/output/Momentum-SelfReports_03-27-2023.csv") %>%
  select(record_id, aff_instability_score, neg_rel_score, pai_total_score, imp_score, aff_score, self_score, total_score) %>%
  rename(id = record_id) %>%
  mutate(id = as.character(id))

subject_df <- left_join(subject_df, selfreport_df, by = "id")

gpa <- setup_glm_pipeline(
  analysis_name = "comp_jun2023", scheduler = "slurm",
  output_directory = "/proj/mnhallqlab/no_backup/ridl_fmri",
  trial_data = trial_df, subject_data = subject_df, run_data = run_df,
  tr = .635, drop_volumes = 2,
  l1_models = NULL, l2_models = NULL, l3_models = NULL,
  n_expected_runs = 4,
  confound_settings = list(
    confound_input_colnames = c("csf", "csf_derivative1", "white_matter", "white_matter_derivative1"), # assumption
    l1_confound_regressors = c("csf", "csf_derivative1", "white_matter", "white_matter_derivative1"),
    na_strings = c("NA", "n/a"), # weird fmriprep-ism for confound files
    exclude_run = "max(FD) > 8 | sum(FD > .9)/length(FD) > .10", #this must evaluate to a scalar per run
    # exclude_subject = "n_good_runs < 2",
    exclude_subject = NULL,
    # truncate_run = "framewise_displacement > 2 & (time > last_offset + )" # 2 seconds after last offset
    truncate_run = "(FD > 0.9 & time > last_offset) | (time > last_offset + last_isi)",
    spike_volumes = NULL
  ),
  parallel = list(
    fsl = list(
      l1_feat_alljobs_time = "72:00:00", # give all l1 runs a while to clear the scheduler
      finalize_time <- "10:00:00" # finalize can be slow with truncation behaviors
    )
  )
)

#l1_models <- gpa$l1_models

gpa <- build_l1_models(gpa, from_spec_file = "ridl_comp_l1.yaml")
gpa <- build_l2_models(gpa)
gpa <- build_l3_models(gpa)
#gpa$parallel$compute_environment <- c("module use /proj/mnhallqlab/sw/modules", "module load r/4.1.2_depend")
saveRDS(gpa, file = "/proj/mnhallqlab/users/michael/ridl_fmri_analysis/code/gpa_snapshot_23Jun2023.rds")

run_glm_pipeline(gpa)

gpa <- readRDS(file = "/proj/mnhallqlab/users/michael/ridl_fmri/basic_apr2022/basic_apr2022.rds")
# for combining results after run completes
res <- combine_feat_l3_to_afni(
  gpa, 
  feat_l3_combined_filename = "{gpa$output_directory}/feat_l3_combined/L1m-{l1_model}/L2m-{l2_model}_L3m-{l3_model}/l2c-{l2_cope_name}_l3c-{l3_cope_name}_stats",
  feat_l3_combined_briknames = "l1c-{l1_cope_name}",
  template_brain = "/proj/mnhallqlab/lab_resources/standard/mni_icbm152_nlin_asym_09c/mni_icbm152_t1_tal_nlin_asym_09c_brain.nii"
)
