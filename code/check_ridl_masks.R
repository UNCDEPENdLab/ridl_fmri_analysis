mask_dir <- "~/Downloads/masks" # MNH copied these off of Longleaf

library(fmri.pipeline)
library(RNifti)
library(ggplot2)
library(dplyr)

repo_dir <- "~/Data_Analysis/ridl_fmri_analysis"

masks <- list.files(mask_dir, pattern=".*brain_mask\\.nii\\.gz", full.names=TRUE, recursive = TRUE)
mask_stats <- generate_mask_diagnostics(masks, 
                                        output_directory = file.path(repo_dir, "output", "mask_diagnostics"),
                                        reference_t1 = file.path(repo_dir, "data", "mni_t1w_2.7mm.nii.gz"),
                                        reference_mask = file.path(repo_dir, "data", "mni_mask_2.7mm.nii.gz"),
                                        generate_run_plots = FALSE) # takes a long time if TRUE!

hist(mask_stats$reference_overlap)
hist(mask_stats$run_not_in_intersect)

# worst overlap with reference
mask_stats %>% slice_min(reference_overlap, n=10)

# smallest difference from intersection mask -- contributing most to deletion
mask_stats %>% slice_min(run_not_in_intersect, n=10)

# smallest masks
mask_stats %>% slice_min(mask_voxels, n=10)

# highest number of voxels not in reference (could be processing problems or poor skull strip)
mask_stats %>% slice_max(run_not_in_reference, n=10)
