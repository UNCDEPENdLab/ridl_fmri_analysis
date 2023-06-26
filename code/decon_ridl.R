library(fmri.pipeline)
library(dplyr)
library(data.table)
library(glue)
library(oro.nifti)
#basedir <- "/proj/mnhallqlab/no_backup/ridl_fmri/comp_jun2023/scheduler_scripts/batch_8c10eee9-72d6-421c-adf6-2ae8b59eccd0"
#load(file.path(basedir, "run_pipeline_cache.RData"))

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


l1_df <- run_df
setDF(l1_df)

atlas <- "/proj/mnhallqlab/lab_resources/parcellation/schaefer_wb_parcellation/Schaefer_444_final_2009c_2.3mm.nii.gz"
metadata <- l1_df %>% select(id, run_number)
l1_niftis <- l1_df$run_nifti # drop 2 TRs

ncpus <- 10 # per job submission
TR <- 0.635

decon_settings = list(
  nev_lr = .01, # neural events learning rate (default in algorithm)
  epsilon = .005, # convergence criterion (default)
  beta = 1, # best for individual differences in MEDuSA sims
  kernel = fmri.pipeline:::spm_hrf(TR)$hrf
) # canonical SPM difference of gammas

runsperproc <- 1 # only one l1 nifti per job processor

# adapted from run_feat_sepjobs
njobs <- ceiling(length(l1_niftis) / (ncpus * runsperproc))
# use length.out on rep to ensure that the vectors align even if chunks are uneven wrt files to run
df <- data.frame(
  l1_niftis = l1_niftis,
  job = rep(1:njobs, each = ncpus * runsperproc, length.out = length(l1_niftis)), stringsAsFactors = FALSE
)

df <- df %>% mutate(metadata = metadata) # add list column
df <- df[order(df$job), ]

# decons for entire 444 in one file are huge and slow to compute. Switch to by network
# submission_id <- basename(tempfile(pattern = "job"))
# joblist <- rep(NA_character_, njobs)
# for (j in seq_len(njobs)) {
#   jind <- which(df$job == j)
#   metadata <- df$metadata[jind,]
#   l1_niftis <- df$l1_niftis[jind]

#   d_batch <- R_batch_job$new(
#     job_name = glue("decon_momentum_{submission_id}"), n_cpus = ncpus, mem_per_cpu = "10g",
#     wall_time = "112:00:00", scheduler = "slurm",
#     # pass relevant vars to the batch
#     input_objects = fmri.pipeline:::named_list(l1_niftis, metadata, atlas, ncpus, TR, decon_settings),
#     r_packages = "fmri.pipeline",
#     r_code = expression({
#       # 0 TR time shift -- deal with dropped volumes and truncation later
#       voxelwise_deconvolution(l1_niftis, metadata,
#         out_dir = "/proj/mnhallqlab/users/michael/momentum_decon",
#         TR = TR, time_offset = 0, atlas_files = atlas, mask = NULL, nprocs = ncpus, save_original_ts = TRUE,
#         algorithm = "bush2011", decon_settings = decon_settings,
#         out_file_expression = expression(paste0("sub", this_subj$id, "_run", this_subj$run_number, "_", atlas_img_name))
#       )
#     })
#   )

#   d_batch$submit()
# }

# by network

atlas_dir <- "/proj/mnhallqlab/lab_resources/parcellation/schaefer_wb_parcellation"
atlas <- "/proj/mnhallqlab/lab_resources/parcellation/schaefer_wb_parcellation/Schaefer_444_final_2009c_2.3mm.nii.gz"
lookup <- "/proj/mnhallqlab/lab_resources/parcellation/schaefer_wb_parcellation/labels/Schaefer_444_region_labels.csv"

orig <- readNIfTI(atlas, reorient=FALSE)
labels <- read.csv(lookup, comment.char="#") %>%
  mutate(network = if_else(roi_num > 400, "Subcortical", network)) # combine all subcort into 'network'
  
#for parallel speed, divide parcellation into one mask per network
for (nn in unique(labels$network)) {
  thisnet <- labels %>% filter(network == nn)

  mod <- orig
  bad_mi <- which(!orig %in% unique(thisnet$roi_num))
  mod[bad_mi] <- 0 #dump

  writeNIfTI(mod, filename = file.path(atlas_dir, "by_network", paste0("Schaefer_444_", nn, "_2.3mm")))
 
}

for (nn in unique(labels$network)) {
  
  submission_id <- basename(tempfile(pattern = "job"))
  joblist <- rep(NA_character_, njobs)

  for (j in seq_len(njobs)) {
    jind <- which(df$job == j)
    metadata <- df$metadata[jind, ]
    l1_niftis <- df$l1_niftis[jind]
    atlas <- file.path(atlas_dir, "by_network", paste0("Schaefer_444_", nn, "_2.3mm.nii.gz"))

    d_batch <- R_batch_job$new(
      job_name = glue("decon_momentum_{nn}_{submission_id}"), n_cpus = ncpus, mem_per_cpu = "10g",
      wall_time = "48:00:00", scheduler = "slurm",
      # pass relevant vars to the batch
      input_objects = fmri.pipeline:::named_list(l1_niftis, metadata, atlas, ncpus, TR, decon_settings),
      r_packages = "fmri.pipeline",
      r_code = expression({
        # 0 TR time shift -- deal with dropped volumes and truncation later
        voxelwise_deconvolution(l1_niftis, metadata,
          out_dir = "/proj/mnhallqlab/users/michael/momentum_decon",
          TR = TR, time_offset = 0, atlas_files = atlas, mask = NULL, nprocs = ncpus, save_original_ts = FALSE,
          algorithm = "bush2011", decon_settings = decon_settings,
          out_file_expression = expression(paste0("sub", this_subj$id, "_run", this_subj$run_number, "_", atlas_img_name))
        )
      })
    )

    d_batch$submit()
  }
}



# d_batch <- R_batch_job$new(
#   job_name = "decon_bsocial", n_cpus = ncpus, mem_per_cpu = '10g',
#   wall_time = "24:00:00", scheduler = "slurm",
#   # pass relevant vars to the batch
#   input_objects = fmri.pipeline:::named_list(l1_niftis, metadata, atlas, ncpus, TR, decon_settings),
#   r_packages = "fmri.pipeline",
#   r_code = expression({
#     # 2 TR time shift
#     voxelwise_deconvolution(l1_niftis, metadata,
#       out_dir = "/proj/mnhallqlab/users/michael/bsocial_decon",
#       TR = TR, time_offset = 2*TR, atlas_files = atlas, mask = NULL, nprocs = ncpus, save_original_ts = TRUE,
#       algorithm = "bush2011", decon_settings=decon_settings,
#       out_file_expression = expression(paste0("sub", this_subj$id, "_run", this_subj$run_number, "_", atlas_img_name))
#     )
#   })
# )

# d_batch$submit()

# aa <- RNifti::readNifti(l1_niftis[1])
# inp <- t(aa[30, 30, 30:31, ])
# system.time(res <- deconvolve_nlreg(inp, decon_settings$kernel, beta = 1))

