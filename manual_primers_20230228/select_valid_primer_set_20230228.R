library(dplyr)
library(magrittr)
library(readr)
library(tidyr)
library(ggplot2)
library(stringr)
options(dplyr.summarise.inform = FALSE)
library(openPrimeR)
library(seqinr)
source("functions_primer_selection.R")

# generating format for primer-dimer input
df_primer <- read_csv("manual_primers_20230228/manual_primers_20230228.csv")
df_primers_formatted <- df_primer %>% 
  mutate(fasta_header = paste(Target, oligo, Rank, sep="|")) %>% 
  select(fasta_header, seq)
write.fasta(as.list(df_primers_formatted$seq), as.list(df_primers_formatted$fasta_header), "manual_primers_20230228/primer-dimer_input.fasta")

# primer-dimer.com -> copy and paste fasta file as input -> select "multiplex analysis" and "excel summary report"

# read in primer-dimer predictions
df_dimer <- read_csv("manual_primers_20230228/dimer_prediction_report_20230228.csv")

# creating symmetric matrix of interactions in long df form
df_dimer_lower <- df_dimer %>% 
  separate_wider_delim("Forward Primer Name", "|", names=c("id_target_a", "direction_a", "id_pair_a")) %>% 
  separate_wider_delim("Reverse Primer Name", "|", names=c("id_target_b", "direction_b", "id_pair_b")) %>% 
  select(id_target_a, direction_a, id_pair_a, id_target_b, direction_b, id_pair_b, dG)

df_dimer_upper <- df_dimer_lower %>% 
  rename(id_target_a = id_target_b,
         id_target_b = id_target_a,
         id_pair_a = id_pair_b,
         id_pair_b = id_pair_a,
         direction_a = direction_b,
         direction_b = direction_a) %>% 
  select(id_target_a, direction_a, id_pair_a, id_target_b, direction_b, id_pair_b, dG) # ensure column order matches lower tri

df_dimer_symmetric <- rbind(df_dimer_lower, df_dimer_upper)

# validating everything is there
df_dimer_symmetric %>% 
  filter((id_target_a == 12 & direction_a == "RIGHT" & id_pair_a == 3 & id_target_b == 7 & direction_b == "LEFT" & id_pair_b == 0) |
         (id_target_a == 7 & direction_a == "LEFT" & id_pair_a == 0 & id_target_b == 12 & direction_b == "RIGHT" & id_pair_b == 3)) %>% 
  View()

# performing hard filtering to identify bad primers
df_valid_primers <- df_dimer_symmetric %>% 
  group_by(id_target_a, direction_a, id_pair_a) %>% 
  summarize(n = n(),
            avg = mean(dG, na.rm = TRUE),
            n_below_neg10 = sum(dG <= -10),
            n_below_neg7 = sum(dG <= -7),
            n_below_neg6 = sum(dG <= -6)) %>%
  filter(avg >= -7) %>% 
  select(id_target_a, direction_a, id_pair_a) %>% 
  rename(target = id_target_a,
         direction = direction_a,
         pair = id_pair_a) %>% 
  mutate(retain = TRUE)

# removing bad primers from both a and b sets
df_valid_pairs <- df_dimer_symmetric %>% 
  left_join(df_valid_primers,
            by=c("id_target_a" = "target",
                 "direction_a" = "direction",
                 "id_pair_a" = "pair"),
            multiple="all") %>% 
  filter(retain) %>% 
  left_join(df_valid_primers,
            by=c("id_target_b" = "target",
                 "direction_b" = "direction",
                 "id_pair_b" = "pair"),
            multiple="all",
            suffix=c(".sym", ".val")) %>% 
  filter(retain.val) %>% 
  select(-retain.sym, -retain.val)

# running first test
# throwing a warning now...
func_recurs_filter_primers(df_valid_pairs, -7) %>% 
  func_count_targets()

# generated a "best guess" first set of primers, with redundant options for some targets
df_selected_set <- func_get_best_run(df_valid_pairs, 300, -9, random=TRUE)

write_csv(df_selected_set, "manual_primers_20230228/init_set_20230301.csv")
df_selected_set <- read_csv("manual_primers_20230228/init_set_20230301.csv")

# filling in missing primers manually
df_missing_oneway <- df_selected_set %>% 
  group_by(id_target_a) %>% 
  summarize(n=n()) %>% 
  filter(n == 1) %>% 
  ungroup

# manually identify all missing target-directions:
# target | direction
# 9 | LEFT
# 13 | RIGHT
# 15 | BOTH
# 14 | BOTH
# 16 | LEFT
# 20 | BOTH

# pull all valid primers for those missing above
df_missing_options <- df_valid_primers %>% 
  filter((target == 9 & direction == "RIGHT") |
         (target == 13 & direction == "RIGHT") |
         (target == 16 & direction == "LEFT") |
         (target == 14) |
         (target == 15) |
         (target == 20))

# options (reformatted results of above):
# 9 | RIGHT | 0,1
# 13 | RIGHT | 0,1,2,3,4
# 15 | LEFT | 0,1
# 15 | RIGHT | 1
# 14 | LEFT | 0,1
# 14 | RIGHT | 0
# 16 | LEFT | none
# 20 | LEFT | 0,1
# 20 | RIGHT | 0,1

# get full set of these options for filling in joined with the original selected set
df_full_set_to_eval <- tibble(id_target_a = c(9,13,13,13,13,13,15,15,15,14,14,14,20,20,16),
                              direction_a = c("LEFT", "RIGHT", "RIGHT", "RIGHT", "RIGHT", "RIGHT", "LEFT", "LEFT", "RIGHT", "LEFT", "LEFT", "RIGHT", "LEFT", "RIGHT","LEFT"),
                              id_pair_a = c(0,0,1,2,3,4,0,1,1,0,1,0,1,0,0),
                              n = 10) %>% 
  rbind(df_selected_set) %>% 
  group_by(id_target_a, direction_a) %>% 
  rename(id_target = id_target_a,
         direction = direction_a,
         id_pair = id_pair_a) %>% 
  mutate(retain=TRUE,
         id_target = as.character(id_target),
         id_pair = as.character(id_pair)) %>% 
  group_by(id_target, direction) %>% # choose best ranked primer for each target-direction pair
  slice_min(order_by=id_pair)
  

# identify interactions across reduced set
df_dimer_symmetric %>% 
  left_join(df_full_set_to_eval, by=c("id_target_a" = "id_target", "direction_a" = "direction", "id_pair_a" = "id_pair")) %>% 
  filter(retain) %>% 
  select(-retain) %>% 
  left_join(df_full_set_to_eval, by=c("id_target_b" = "id_target", "direction_b" = "direction", "id_pair_b" = "id_pair")) %>% 
  filter(retain) %>%
  select(-retain) %>% 
  arrange(dG) %>% 
  View()

# overview of decisions:
# ultimately chose 20L1

# remaining problematic interactions:
# 1R2 <-> 3L0
# 12R1 <-> 14R0
# 8R2 <-> 14R0
# 10R1 <-> 19R0
# 6R0 <-> 13R0
# 9L0 <-> 5R0

# alternate choices for the above, based on df_selected_set:
# 1R4, 3L1, 3L3
# 12R2
# no alt option for 14R0
# 8R3, 8R4
# 10R2, 10R3, 19R2 -> choosing 10R2 resolves
# 6R2, 13R1,2,3,4 -> all of these options led to worse interactions, best to keep as-is
# 9L1, 5R1 -> use 9L1 but keep 5R0, still have interaction but it's slightly less bad

# implement desired replacements to improve selected set
df_replacement <- tibble(id_target = c("3", "10", "8", "1", "4", "9"),
                         direction = c("LEFT", "RIGHT", "RIGHT", "LEFT", "RIGHT", "LEFT"),
                         id_pair = c("1", "2", "3", "0", "3", "1"),
                         n = 10,
                         retain = TRUE) %>% 
  rbind(df_full_set_to_eval) %>% 
  filter(!(id_target == "3" & direction == "LEFT" & id_pair == "0"),
         # !(id_target == "14" & direction == "RIGHT" & id_pair == "0"),
         !(id_target == "10" & direction == "RIGHT" & id_pair == "1"),
         !(id_target == "8" & direction == "RIGHT" & id_pair == "2"),
         !(id_target == "4" & direction == "RIGHT" & id_pair == "1"),
         !(id_target == "1" & direction == "LEFT" & id_pair == "1"),
         !(id_target == "9" & direction == "LEFT" & id_pair == "0"),
         # !(id_target == "5" & direction == "RIGHT" & id_pair == "0"),
         ) 
  # filter(!(id_target == "1" & direction == "RIGHT" & id_pair == "2"))

df_dimer_symmetric %>% 
  left_join(df_replacement, by=c("id_target_a" = "id_target", "direction_a" = "direction", "id_pair_a" = "id_pair")) %>% 
  filter(retain) %>% 
  select(-retain) %>% 
  left_join(df_replacement, by=c("id_target_b" = "id_target", "direction_b" = "direction", "id_pair_b" = "id_pair")) %>% 
  filter(retain) %>%
  select(-retain) %>% 
  arrange(dG) %>% 
  View()

write_csv(df_replacement, "manual_primers_20230228/chosen_primer_set_4.csv")
df_chosen_primers <- read_csv("manual_primers_20230228/chosen_primer_set_1.csv")
df_primer <- read_csv("manual_primers_20230228/manual_primers_20230228.csv")

df_chosen_primers %>% 
  mutate(retain = TRUE) %>% 
  right_join(df_primer, by = c("id_target" = "Target", "direction" = "oligo", "id_pair" = "Rank")) %>% 
  filter(retain) %>% 
  View()

# problems
# 14R has to be 0 or 1 -> set to 0, and change 8R2 -> 8R3 (no solution to help interaction with 12R1; using 12R2 kept interaction and 12R3 had bad interactions with others)
# 16L has to be 1,2,3 -> 
## 16L3 has -7.04 for everything except with 18R1 at -7.14
## 16L1 has -7.01 with one interaction with 3L1 at -7.17
## 16L2 has -7.24
## 16L0 has -6.99, but has one interaction with 13R0 at -11.7
# 20R has to be 0,1

# if set 20R0, need to change 5R0 and 14R0 (but can't change 14R0). 5R1 is only other option and it has numerous bad interactions
# if set 20R1, need to change 4R1, 15R1, 1L1
# changing 4R1 to ...
## 4R0 <-> 3L1 -7.7; 
## 4R2 <-> 5L1 -7.53; 
## 4R3 <-> 4L4 -7.05;
## 4R4 <-> 5L1 -7 and 4R4 <-> 18R1 -8.13
# 15R: 15R0 is only other option and it has -10.8 for almost everything
# 1L: 1L1 resolves 20R1 interaction

# let's just try some combos in openPrimeR i guess?
# unfortunately can't set up openPrimeR environments anymore -- dependency for cross-dimerization calculations no longer supported

# visualize boxplots of different sets' dGs
source("functions_check_opr_constraints.R")
func_eval_primers_in_dir("manual_primers_20230228/chosen_primer_sets", "manual_primers_20230228/manual_primers_20230228.csv", "manual_primers_20230228/chosen_primer_fastas")


# all sets passed in openPrimeR! picking my favorite:
chosen_primers_1 <- read_csv("manual_primers_20230228/chosen_primer_sets/chosen_primer_set_1.csv")
chosen_primers_2 <- read_csv("manual_primers_20230228/chosen_primer_sets/chosen_primer_set_2.csv")
chosen_primers_3 <- read_csv("manual_primers_20230228/chosen_primer_sets/chosen_primer_set_3.csv")
chosen_primers_4 <- read_csv("manual_primers_20230228/chosen_primer_sets/chosen_primer_set_4.csv")

df_dimer_symmetric %>% 
  mutate(id_target_a = as.numeric(id_target_a),
         id_target_b = as.numeric(id_target_b),
         id_pair_a = as.numeric(id_pair_a),
         id_pair_b = as.numeric(id_pair_b)) %>% 
  left_join(chosen_primers_3, by=c("id_target_a" = "id_target", "direction_a" = "direction", "id_pair_a" = "id_pair")) %>% 
  filter(retain) %>% 
  select(-retain) %>% 
  left_join(chosen_primers_3, by=c("id_target_b" = "id_target", "direction_b" = "direction", "id_pair_b" = "id_pair")) %>% 
  filter(retain) %>%
  select(-retain) %>% 
  arrange(dG) %>% 
  View()
