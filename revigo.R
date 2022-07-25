# R script for programtic access to Revigo. Run it with (last output file name is optional):
# Rscript revigo.R example.csv result.csv

library(httr)
library(rvest)

revigo <- function(df, cutoff, valueType, species){
  
  ## Reformat dataframe into a csv-like string
  csv <- df %>% 
    tidyr::unite("combined", sep = "\t") %>% 
    dplyr::pull(combined) %>% 
    stringr::str_c(collapse = "\n")
  
  
  ## Maximum number of columns = 2 (GO ID and associated score)
  if (!(ncol(df) %in% c(1,2))) {
    stop("Your input dataframe must contain 1 or 2 columns: GO IDs and/or associated scores")
  }
  
  
  ## Cutoff info
  cutoff_opts <- c("large", "medium", "small", "tiny")
  cutoff_num_opts <- c(0.9, 0.7, 0.5, 0.4)
  if (missing(cutoff)) {
    cutoff_ <- "0.7"
    warning("Your cutoff is set to 'Medium (0.7)'. \n Other options could be \n - 'Large' - 0.9 \n - 'Small' - 0.5 \n - 'Tiny' - 0.4 \n To set the cutoff use key word or number")
  }
  else if (tolower(cutoff) %in% cutoff_opts) {
    cutoff_ <- as.character(cutoff_num_opts[match(tolower(cutoff), cutoff_opts)])
    # message("cutoff is OK")
  }
  else if (cutoff %in% cutoff_num_opts) {
    cutoff_ <- as.character(cutoff)
    # message("cutoff is OK")
  }
  else {
    cutoff_ <- "0.7"
    warning("Your cutoff was set wrong. Now it is set to 'Medium (0.7)'. \n Other options could be \n - 'Large' - 0.9 \n - 'Small' - 0.5 \n - 'Tiny' - 0.4 \n To set the cutoff use key word or number")
  }
  
  
  ## Value type info
  valueType_opts <- c("pvalue", "higher", "lower", "higherabsolute", "higherabslog2")
  if (ncol(df) > 1) {
    ## Notify if no valueType was specified
    if (missing(valueType)) {
      valueType_ <- "pvalue"
      warning("Your valueType is set to 'PValue'. \n Other options could be \n - 'Higher' - 'Higher value is better' \n - 'Lower' - Lower value is better \n - 'HigherAbsolute' - Higher absolute value is better \n - 'HigherAbsLog2' - Higher absolute log2 value is better")
    }
    ## Warn if wrong key was used
    else if (!(tolower(valueType) %in% valueType_opts)) {
      valueType_ <- "pvalue"
      warning("Your valueType was set wrong. Now it is set to 'PValue'. \n Other options could be \n - 'Higher' - 'Higher value is better' \n - 'Lower' - Lower value is better \n - 'HigherAbsolute' - Higher absolute value is better \n - 'HigherAbsLog2' - Higher absolute log2 value is better")
    }
    else {
      valueType_ <- tolower(valueType)
      # message("valueType is OK")
    }
  }
  
  
  ## Organism taxon info
  if(missing(species)){
    speciesTaxon_ <- "0"
  }
  else{
    if(species ==  "human"){
      speciesTaxon_ <- "9606"
    }
    if(species == "mouse"){
      speciesTaxon_ <- "10090" 
    }
    if(!(species %in% c("human", "mouse"))){
      speciesTaxon_ <- "0"
      warning("Your species should be 'human' or 'mouse'. Now it changed to Whole UniProt database (default)")
    }
  }
  

  # Submit job to Revigo
  httr::POST(
    url = "http://revigo.irb.hr/Revigo.aspx",
    body = list(
      cutoff = cutoff_,
      valueType = valueType_,
      speciesTaxon = speciesTaxon_,
      measure = "SIMREL",
      goList = csv
    ),
    # application/x-www-form-urlencoded
    encode = "form"
  ) -> res
  
  dat <- httr::content(res, encoding = "UTF-8")
  dat <- rvest::html_elements(dat, "#BiologicalProcess") 
  dat <- rvest::html_table(dat)
  
  dat[[1]]

}