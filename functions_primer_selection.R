# functions for choosing preferred primer set from pairwise dimer dG values


func_return_higher_id <- function(id_pair_a, id_pair_b, random=FALSE){
  #' accepts two integers and boolean for random state
  #' returns string corresponding to higher int position, or random choice

  if (random) return(sample(c("a", "b"), 1))
  
  if (id_pair_a > id_pair_b) {
    return("a")
  } else if (id_pair_a < id_pair_b) {
    return("b")
  } else {
    return(sample(c("a", "b"), 1))
  }
}

Vfunc_return_higher_id <- Vectorize(func_return_higher_id)

func_find_primers_to_drop <- function(df_valid_pairs, random=FALSE){
  #' identify worst primer for each target-direction

  df_valid_pairs %>% 
    group_by(id_target_a, direction_a, id_pair_a) %>% # choose worst interaction for each primer
    slice_min(order_by=dG) %>%
    group_by(id_target_a, direction_a) %>% # choose worst interaction for each target/dir combo
    slice_min(order_by=dG) %>% 
    sample_n(1) %>% # if still have multiple id_pairs per target/dir, select one randomly
    mutate(primer_to_drop = Vfunc_return_higher_id(id_pair_a, id_pair_b, random)) %>% # for each row, determine which is higher: id_pair_a or id_pair_b
    rename_with(.cols=starts_with('id_'), .fn=str_remove, pattern="id_") %>% 
    pivot_longer(-c(dG, primer_to_drop), 
                 names_to = c(".value", "Var"), 
                 names_sep="_" ) %>% # pivoting so i can easily grab the row we want to drop
    filter(Var == primer_to_drop) %>% 
    select(target, direction, pair) %>% 
    rename(target = target,
           direction = direction,
           pair = pair) %>% 
    mutate(to_drop = TRUE)
}

func_remove_primers <- function(df_valid_pairs, df_primers_to_drop){
  #' remove all occurrences of any primer to drop from the symmetric df of valid primer pairs

  df_valid_pairs %>% 
    left_join(df_primers_to_drop,
              by=c("id_target_a" = "target",
                   "direction_a" = "direction",
                   "id_pair_a" = "pair"),
              multiple="all") %>% 
    filter(is.na(to_drop)) %>%
    select(-to_drop) %>% 
    left_join(df_primers_to_drop,
              by=c("id_target_b" = "target",
                   "direction_b" = "direction",
                   "id_pair_b" = "pair"),
              multiple="all"
    ) %>% 
    filter(is.na(to_drop)) %>% 
    select(-to_drop)
}

func_count_targets <- function(df_filtered_pairs){
  #' count how many primer options remain for each target-direction
  df_filtered_pairs %>% 
    group_by(id_target_a, direction_a) %>% 
    summarize(n = n()) %>% 
    nrow()
}

func_recurs_filter_primers <- function(df_valid_pairs, threshold, random=FALSE){
  #' Iteratively remove primers from the worst interaction pairs until all remaining interactions exceed threshold
  if (min(df_valid_pairs$dG) >= threshold){
    df_valid_pairs
  } 
  
  else {
    func_find_primers_to_drop(df_valid_pairs, random) %>% 
      func_remove_primers(df_valid_pairs, .) %>%
      func_recurs_filter_primers(., threshold, random)
  }
}

helper_func_count_targets <- function(x, array_all_out){
  array_all_out[,x] %>% 
    as.data.frame() %>% 
    func_count_targets()
}

func_get_best_run <- function(df_valid_pairs, i, threshold, random=FALSE){
  #' identify primer sets through recursion i times, then return the set with the most targets included
  
  array_all_out <- replicate(i, func_recurs_filter_primers(df_valid_pairs, threshold, random))
  
  n_targets <- sapply(seq(1,i), helper_func_count_targets, array_all_out)
  
  array_all_out[,which.max(n_targets)] %>% 
    as.data.frame() %>% 
    group_by(id_target_a, direction_a, id_pair_a) %>% 
    summarize(n = n())
}