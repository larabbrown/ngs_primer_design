# constraint functions

func_join_primer_to_seq <- function(chosen_set_path, primer_full_path){
  #' DANGER ZONE: designed to overwrite files
  #' Only run once, to add sequences to the filtered primer set (only because i was a fool and did not do this while creating them)

  read_csv(chosen_set_path) %>% 
    mutate(retain = TRUE) %>% 
    right_join(read_csv(primer_full_path), by = c("id_target" = "Target", "direction" = "oligo", "id_pair" = "Rank")) %>% 
    filter(retain) %>% 
    write_csv(chosen_set_path)
}

func_get_primer_fasta <- function(chosen_set_path, primer_full_path, destination_dir){
  #' convert primer3 CSV to fasta
  fasta_out_name <- paste0(destination_dir, "/", tools::file_path_sans_ext(basename(chosen_set_path)), ".fasta")
  
  df <- read_csv(chosen_set_path) %>% 
    mutate(header_str = paste(paste(id_target, direction, sep="_"), id_pair, sep="|"),
           lower_seq = tolower(seq))
  
  seqinr::write.fasta(as.list(df$lower_seq), as.list(df$header_str), fasta_out_name)
  
  if(!file.exists(fasta_out_name)){

    df <- read_csv(chosen_set_path) %>%
      mutate(header_str = paste(paste(id_target, direction, sep="_"), id_pair, sep="|"),
             lower_seq = tolower(seq))

    seqinr::write.fasta(as.list(df$lower_seq), as.list(df$header_str), fasta_out_name)
  }
  
}

func_analysis_settings <- function(){
  #' set settings for OpenPrimeR analyses, notable cross_dimerization threshold
  
  list.files(system.file("extdata", "settings", package = "openPrimeR"), pattern = "*\\.xml")
  settings.xml <- system.file("extdata", "settings", 
                              "C_Taq_PCR_high_stringency.xml", package = "openPrimeR")
  settings <- read_settings(settings.xml)
  design.settings <- settings
  constraints(design.settings)[["self_dimerization"]] = c("min"=-6, "max"=0)
  constraints(design.settings)[["cross_dimerization"]] = c("min"=-7, "max"=0)
  out.file <- tempfile("interaction_settings", fileext = ".xml")
  write_settings(design.settings, out.file)
  analysis.settings <- read_settings(out.file)
  setting_list = list("settings" = settings, "analysis_settings" = analysis.settings)
  return(setting_list)
}

func_eval_primers <- function(fasta_primers, analysis_settings){
  #' Compare a set of primers from a fasta file against the constraints set for OpenPrimeR

  primers <- read_primers(fasta_primers, fw.id = "_LEFT", rev.id = "_RIGHT")
  templates <- read_templates("templates.fasta")
  
  check_constraints(primers, templates, analysis_settings, active.constraints = c("self_dimerization", "cross_dimerization"))
  
}

func_eval_primers_in_dir <- function(primer_directory, primer_full_path, primer_fasta_destination_dir){
  #' Accept directory of primer csvs, a path to file containing initial primer3 output, and path to store primer fastas
  #' Generate fasta files from primer CSV
  #' Evaluate how well each primer fasta meets desired constraints
  #' Visualize dimerization constraint box plots

  purrr::map(list.files(path=primer_directory, full.names = TRUE), func_get_primer_fasta, primer_full_path, primer_fasta_destination_dir)
  
  setting_list <- func_analysis_settings()
  
  constraint_list <- purrr::map(list.files(path=primer_fasta_destination_dir, full.names = TRUE), func_eval_primers, setting_list$analysis_settings)
  
  plot_constraint(constraint_list, setting_list$settings, active.constraints = c("self_dimerization", "cross_dimerization"))
  
}