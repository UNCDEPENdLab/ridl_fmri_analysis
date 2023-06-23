library(data.table)
library(tidyverse)
repo_dir <- "~/Data_Analysis/ridl_fmri_analysis"
setwd(repo_dir)

source(file.path(repo_dir, "code/parse_ridl_all.R"))
source(file.path(repo_dir, "code/parse_ridl.R"))

mout <- data.table::fread(file.path(repo_dir, "data/modeling_results_table_US_updated.csv.gz"))
fmri_blocks <- mout %>% filter(block %in% c(1000, 1001, 1002, 1003)) %>%
  mutate(
    trial = trial + 1, # avoid 0-based indexing
    outcome = factor(o, levels=c(-1,0,1), labels=c("loss", "nothing", "win")),
    across(c(RPE, Q_chosen, Q_unchosen, Q_chosen_new, delta_Q), 
           .fns=\(x) DescTools::Winsorize(x, probs=c(.005, .995), na.rm=TRUE), .names = "{.col}_wins"),
    run_number = block - 999,
    Q_diff = Q_chosen_wins - Q_unchosen_wins
  )

# View(fmri_blocks)
accrej <- fmri_blocks %>%
  group_by(subject_name, block) %>%
  summarise(accrej_tot=sum(stim2==-1000))

table(accrej$accrej_tot)
accrej %>% filter(accrej_tot == 0)

# oddity of no acc/rej trials?
fmri_blocks %>% filter(subject_name==440612)

# ggplot(fmri_blocks, aes(x=Q_chosen)) + geom_histogram()
# ggplot(fmri_blocks, aes(x=Q_chosen_wins)) + geom_histogram()
# ggplot(fmri_blocks, aes(x=Q_chosen_new)) + geom_histogram()
# ggplot(fmri_blocks, aes(x=Q_unchosen)) + geom_histogram()
# ggplot(fmri_blocks, aes(x=Q_unchosen_wins)) + geom_histogram()
# ggplot(fmri_blocks, aes(x=Q_diff)) + geom_histogram()
# ggplot(fmri_blocks, aes(x=delta_Q)) + geom_histogram()
# ggplot(fmri_blocks, aes(x=RPE)) + geom_histogram()
# ggplot(fmri_blocks, aes(x=RPE, y=delta_Q)) + geom_point()
# ggplot(fmri_blocks, aes(x=RPE_wins)) + geom_histogram()
# table(fmri_blocks$acc)

#oddity of stim2 != -1000 despite having feedback == 0
xtabs(~stim2 + feedback, fmri_blocks)

curious_stim2 <- fmri_blocks %>%
  filter(stim2 != -1000 & feedback == 0)

write.csv(curious_stim2, file = file.path(repo_dir, "output/curious_stim2_feedback0.csv"), row.names=FALSE)
table(curious_stim2$subject_name)

# output for fmri.pipeline
fmri_blocks <- fmri_blocks %>%
  rename(id=subject_name) %>%
  dplyr::select(id, run_number, trial, outcome, RPE, RPE_wins, Q_chosen, Q_chosen_wins, Q_unchosen, Q_unchosen_wins, Q_diff)
  
trial_df <- parse_ridl_all(file.path(repo_dir, "data/momentum"), 
                           matlab_dir="/Applications/MATLAB_R2021b.app/bin")

# divergences between fMRI files and EMA-derived computational model outputs
setdiff(unique(trial_df$id), unique(fmri_blocks$id))
setdiff(unique(fmri_blocks$id), unique(trial_df$id))

length(unique(trial_df$id))
length(unique(fmri_blocks$id))

head(fmri_blocks)

#trial_df <- trial_df %>% left_join(fmri_blocks, by=c("id", "run_number", "trial"))

# For now, use inner join since we have a few subjects in trial_df not present in the computational model outputs
trial_df <- trial_df %>% inner_join(fmri_blocks, by=c("id", "run_number", "trial"))

trial_df_fmri <- trial_df %>% select(
  id, run_number, trial, trial_type, outcome_fac, 
  stim_time, choice_time, reaction_time, feedback_time, feedback_isi, iti_actual,
  RPE_wins, Q_chosen_wins, Q_unchosen_wins, Q_diff
) %>%
  rename(choice_onset = stim_time, feedback_onset = feedback_time) %>%
  mutate(feedback_duration = 1.0)

xtabs(~is.na(RPE_wins) + trial_type, trial_df_fmri)

# anomalous cases where RPE is missing on a new learn trial
# trial_df_fmri %>% filter(is.na(RPE_wins) & trial_type == "new_learn") %>% View()
trial_df_fmri %>% filter(is.na(RPE_wins) & trial_type == "new_learn") %>% pull(id) %>% table()

# anomalous cases where Q_chosen is missing on a trial
trial_df_fmri %>% filter(is.na(Q_chosen_wins)) %>% group_by(run_number, trial) %>% tally() %>% arrange(desc(n))

# anomalous cases where outcome coding diverges between computational model outputs and MATLAB trial data
trial_df %>% filter(outcome.y != outcome_fac) %>% pull(id) %>% table()
outcome_coding_diverges <- trial_df %>% filter(outcome.y != outcome_fac) %>%
  select(id, run_number, trial, block, stim1, stim2, outcome.x, outcome.y, outcome_fac)

write.csv(outcome_coding_diverges, file = file.path(repo_dir, "output/different_outcome_coding.csv"), row.names=FALSE)

xtabs(~outcome_fac + trial_type, trial_df_fmri)

fwrite(trial_df_fmri, file=file.path(repo_dir, "output/ridl_combined_fmri.csv.gz"))