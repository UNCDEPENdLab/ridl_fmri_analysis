library(dplyr)
library(emmeans)
library(fmri.pipeline)
library(readr)
library(glue)

# get trial data
setwd("/proj/mnhallqlab/projects/ridl_fmri_analysis/code")
trial_df <- read.csv("../output/ridl_combined_fmri.csv.gz") %>%
  mutate(id=as.character(id))

run_data_from_bids <- function(bids_dir, modality="func", type="task", task_name="ridl") {
  checkmate::assert_directory_exists(bids_dir)
  sub_dirs <- list.dirs(bids_dir, recursive = FALSE)
  slist <- lapply(sub_dirs, function(ss) {
    expect_dir <- file.path(ss, modality)
    if (!checkmate::test_directory_exists(expect_dir)) {
      warning(glue("Cannot find expected modality directory: {expect_dir}"))
      return(NULL)
    }

    # I wonder if these could diverge in order
    nii_files <- Sys.glob(glue("{expect_dir}/sub*_{type}-{task_name}*_postproccessed.nii.gz"))
    if (length(nii_files) == 0L) {
      warning(glue("No NIfTI file matches in: {expect_dir}"))
      return(NULL)
    }
    confound_files <- Sys.glob(glue("{expect_dir}/sub*_{type}-{task_name}*-confounds*.tsv"))

    if (length(nii_files) != length(confound_files)) {
      warning(glue("Cannot align nifti and confound files for {expect_dir}"))
      return(NULL)
    }
    id <- sub("^sub-", "", basename(ss))
    run_number <- as.integer(sub(glue(".*sub-.*_{type}-{task_name}(\\d+).*"), "\\1", nii_files, perl = TRUE)) # NB. this is not a BIDS-compliant approach. Needs to be _run-01
    data.frame(id = id, run_number = run_number, run_nifti = nii_files, confound_input_file = confound_files)
  })

  bind_rows(slist)
}

run_df <- run_data_from_bids("/proj/mnhallqlab/studies/momentum/clpipe/data_postproc2")

xtabs(~id, run_df)

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


