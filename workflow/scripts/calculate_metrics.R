##########################################
## This script is to summarize and generate new metric from the hap.py results_extended.csv.  
## This script assumes hap.py was run with and without target-regions = dip.bed and stratification 
## with GIAB v2 stratification BEDs.
############################################

library(tidyverse)
library(here)
library(fs)

############################################
## Input parameter parsing #################
############################################

## TODO - add code for taking targeted and non-targeted directory paths from
## command line, and snakemake

# ## directory paths for results_extended.csv files
# data_dir_targeted <- here("results", "happy", "targeted")
# data_dir_nontargeted <- here("results", "happy", "nontargeted")

## Stratifications used in metric calculations
strats <- c("*" ,
             "GRCh38_HG002_GIABv4.1_notin_complexandSVs_alldifficultregions.bed.gz" ,
             "GRCh38_notinalldifficultregions.bed.gz",
             "GRCh38_notinalllowmapandsegdupregions.bed.gz",
             "GRCh38_notinsegdups.bed.gz",
             "GRCh38_notinAllTandemRepeatsandHomopolymers_slop5.bed.gz",
             "GRCh38_notinAllHomopolymers_gt6bp_imperfectgt10bp_slop5.bed.gz",
             "GRCh38_segdups.bed.gz",
             "GRCh38_AllHomopolymers_gt6bp_imperfectgt10bp_slop5.bed.gz",
             "GRCh38_MHC.bed.gz")

## Path and file to save metric results to 
# outfile <- here("results", "asm_benchmarking_metrics.tsv")

######################################################################## 
### Functions for loading targeted and non-targeted benchmarking results
########################################################################

## TODO - fix read_csv load error with titv column formatting

#' Load non-targeted extended.csv benchmarking results
#' Reads in only the metrics needed to calculate non-targeted metrics for one or
#' more assemblies. Assumes the assembly id is the part of the extended.csv file
#' name, before the first `_`.
#' 
#' @param dir character string with relative path to non-targeted benchmarking results
#'
#' @return tibble 
#' @export
#'
#' @examples read_nt_ext_results("results/happy/nontargeted") 
read_nontargeted_ext_results <- function(dir) {
  extended_files <- fs::dir_ls(dir, regexp = "extended.csv$")
  ext_results_df <- extended_files %>%  
    ## Extracting assembly ID
    set_names(str_remove(basename(.),"_.*")) %>% 
    map_dfr(read_csv, .id="asm_id") 
  
  ## Only including metric columns and row relevant to assembly benchmarking
  ## metric calculations
  ext_results_tidy <- ext_results_df %>%
    filter(Filter %in% c("PASS", "ALL"), Subtype == "*")  %>%
    select(
      asm_id,
      Filter,
      Type,
      Subtype,
      Subset,
      TRUTH.TOTAL.het,
      TRUTH.TOTAL.homalt,
      TRUTH.TOTAL
    ) %>%
    # Providing unique names for non-targeted metrics
    rename(
      ntTRUTH.TOTAL.homalt = TRUTH.TOTAL.homalt,
      ntTRUTH.TOTAL.het    = TRUTH.TOTAL.het,
      ntTRUTH.TOTAL        = TRUTH.TOTAL
    )

  return(ext_results_tidy)
}


#' Load targeted extended.csv benchmarking results
#' Reads in only the metrics needed to calculate targeted metrics for one or
#' more assemblies. Assumes the assembly id is the part of the extended.csv file
#' name, before the first `_`.
#' 
#' @param dir character string with relative path to targeted benchmarking results
#'
#' @return tibble 
#' @export
#'
#' @examples read_targeted_ext_results("results/happy/targeted") 
read_targeted_ext_results <- function(dir) {
  extended_files <- fs::dir_ls(dir, regexp = "extended.csv$")
  ext_results_df <- extended_files %>%  
    ## Extracting assembly ID
    set_names(str_remove(basename(.),"_.*")) %>% 
    map_dfr(read_csv, .id="asm_id")
  
  ## Only including metric columns and row relevant to assembly benchmarking
  ## metric calculations
  ext_results_tidy <- ext_results_df %>%
    select(
      asm_id,
      Filter,
      Type,
      Subtype,
      Subset,
      QUERY.FP,
      Subset.IS_CONF.Size,
      FP.gt,
      FP.al,
      METRIC.Recall,
      TRUTH.TOTAL,
      TRUTH.TP,
      TRUTH.TP.het,
      TRUTH.TP.homalt,
      TRUTH.TOTAL.het,
      TRUTH.TOTAL.homalt
    )

  return(ext_results_tidy)
}


#' Reformatting benchmarking results for assembly metric calculations
#' This function performs three rols
#' 1) subset to desired stratifications 
#' 2) rename colnames to specific metric for Type 
#' 3) combine metric types into single table where all type.metric cols are present. 
#' 
#' @param combined_results data frame with combined targeted and non-targeted assembly metrics
#' @param strats vector of stratifications to subset the extended results
#'
#' @return
#' @export
#'
#' @examples
reformat_results <- function(combined_results, strats) {
  ## TODO - add check that strats are in combined_results$Subset
  
  SNP_ext_results <- combined_results %>%
    filter(Type == "SNP", 
           Subset %in% strats) %>%
    rename(
      SNP.QUERY.FP = QUERY.FP,
      SNP.Subset.IS_CONF.Size = Subset.IS_CONF.Size,
      SNP.FP.gt = FP.gt,
      SNP.FP.al = FP.al,
      SNP.Recall = METRIC.Recall,
      SNP.TRUTH.TOTAL = TRUTH.TOTAL,
      SNP.TRUTH.TP = TRUTH.TP,
      SNP.TRUTH.TP.het = TRUTH.TP.het,
      SNP.TRUTH.TP.homalt = TRUTH.TP.homalt,
      SNP.TRUTH.TOTAL.het = TRUTH.TOTAL.het,
      SNP.TRUTH.TOTAL.homalt = TRUTH.TOTAL.homalt,
      SNP.ntTRUTH.TOTAL.het = ntTRUTH.TOTAL.het,
      SNP.ntTRUTH.TOTAL.homalt = ntTRUTH.TOTAL.homalt,
      SNP.ntTRUTH.TOTAL = ntTRUTH.TOTAL
    ) %>%
    select(-Type, -Subtype)
  
  INDEL_ext_results <- combined_results %>%
    filter(Type == "INDEL",
           Subtype == "*",
           Subset %in% strats) %>%
    rename(
      INDEL.QUERY.FP = QUERY.FP,
      INDEL.Subset.IS_CONF.Size = Subset.IS_CONF.Size,
      INDEL.FP.gt = FP.gt,
      INDEL.FP.al = FP.al,
      INDEL.Recall = METRIC.Recall,
      INDEL.TRUTH.TOTAL = TRUTH.TOTAL,
      INDEL.TRUTH.TP = TRUTH.TP,
      INDEL.TRUTH.TP.het = TRUTH.TP.het,
      INDEL.TRUTH.TP.homalt = TRUTH.TP.homalt,
      INDEL.TRUTH.TOTAL.het = TRUTH.TOTAL.het,
      INDEL.TRUTH.TOTAL.homalt = TRUTH.TOTAL.homalt,
      INDEL.ntTRUTH.TOTAL.het = ntTRUTH.TOTAL.het,
      INDEL.ntTRUTH.TOTAL.homalt = ntTRUTH.TOTAL.homalt,
      INDEL.ntTRUTH.TOTAL = ntTRUTH.TOTAL
    ) %>%
    select(-Type, -Subtype)
  
  ext_results_combined_renamed <- left_join(
    INDEL_ext_results,
    SNP_ext_results
  )
  
  return(ext_results_combined_renamed)
}


#' Calculate Assembly Benchmarking Metrics
#' Function calculates new summary metrics using values form extendded results
#' 
#' @param results 
#'
#' @return
#' @export
#'
#' @examples summary_metrics(combined_results)
calculate_asm_metrics <- function(results_df) {                                                                    

    ## Calculating QV (correctness) metrics all QV metric denominators
    ## multiplied by 2 to be more representative of diploid value
    qv_metrics <- results_df %>% 
      mutate(QV_dip_snp = 
               -10*log10(
                  SNP.QUERY.FP/(2 * SNP.Subset.IS_CONF.Size))
             ) %>%
      mutate(QV_dip_indel = 
               -10*log10(
                 INDEL.QUERY.FP/(2 * INDEL.Subset.IS_CONF.Size))
             ) %>%
      mutate(QV_dip_snp_indel = 
               -10*log10(
                  (SNP.QUERY.FP + INDEL.QUERY.FP)/(2 * SNP.Subset.IS_CONF.Size))
             ) %>%
      mutate(QV_hap_snp = 
              -10*log10(
                    (SNP.QUERY.FP - SNP.FP.gt -SNP.FP.al)/
                    (2 * SNP.Subset.IS_CONF.Size)
                )
             ) %>%
      mutate(QV_ignoreGT_snp = 
               -10*log10(
                    (SNP.QUERY.FP - SNP.FP.gt)/
                    (2 * SNP.Subset.IS_CONF.Size)
                 )
             ) %>%
      mutate(QV_hap_indel = 
               -10*log10(
                    (INDEL.QUERY.FP - INDEL.FP.gt - INDEL.FP.al)/
                    (2 * INDEL.Subset.IS_CONF.Size)
                 )
             ) %>%
      mutate(QV_ignoreGT_indel = 
               -10*log10(
                 (INDEL.QUERY.FP - INDEL.FP.gt)/
                 (2 * INDEL.Subset.IS_CONF.Size)
                )
             ) %>%
      mutate(QV_hap_snp_indel = 
               -10*log10(
                  ((SNP.QUERY.FP - SNP.FP.gt - SNP.FP.al) + 
                     (INDEL.QUERY.FP - INDEL.FP.gt - INDEL.FP.al))/
                  (2 * SNP.Subset.IS_CONF.Size)
               )
             ) %>%
      mutate(QV_ignoreGT_snp_indel = 
               -10*log10(
                  ((SNP.QUERY.FP - SNP.FP.gt) + (INDEL.QUERY.FP - INDEL.FP.gt))/
                                 (2 * SNP.Subset.IS_CONF.Size)
                )
             )
    
    ## Calculating recall or completeness metrics
    recall_metrics <- results_df %>% 
      mutate(SNP.Recall_ignoreGT = (SNP.TRUTH.TP + SNP.FP.gt) / 
                                    SNP.TRUTH.TOTAL
             ) %>%
      mutate(INDEL.Recall_ignoreGT = (INDEL.TRUTH.TP + INDEL.FP.gt) / 
                                      INDEL.TRUTH.TOTAL
             ) %>%
      mutate(SNP.Recall.het = SNP.TRUTH.TP.het / 
                              SNP.TRUTH.TOTAL.het
             ) %>%
      mutate(SNP.Recall.hom = SNP.TRUTH.TP.homalt /
                              SNP.TRUTH.TOTAL.homalt
             ) %>%
      ## TODO - JM and JZ Question - changed from SNP assume want INDEL counts
      mutate(INDEL.Recall.het = INDEL.TRUTH.TP.het / 
                                INDEL.TRUTH.TOTAL.het
             ) %>%
      mutate(INDEL.Recall.hom = INDEL.TRUTH.TP.homalt / 
                                INDEL.TRUTH.TOTAL.homalt
             ) %>%
      mutate(SNP.Recall.het.TOTALnontarget = SNP.TRUTH.TP.het / 
                                             SNP.ntTRUTH.TOTAL.het
             ) %>%
      mutate(SNP.Recall.hom.TOTALnontarget = SNP.TRUTH.TP.homalt / 
                                             SNP.ntTRUTH.TOTAL.homalt
             ) %>%
      mutate(INDEL.Recall.het.TOTALnontarget = SNP.TRUTH.TP.het / 
                                               SNP.ntTRUTH.TOTAL.het
             ) %>%
      mutate(INDEL.Recall.hom.TOTALnontarget = SNP.TRUTH.TP.homalt / 
                                               SNP.ntTRUTH.TOTAL.homalt
             ) %>%
      mutate(SNP.Recall_ignoreGT.TOTALnontarget = (SNP.TRUTH.TP + SNP.FP.gt) / 
                                                    SNP.ntTRUTH.TOTAL
             ) %>%
      mutate(SNP.Recall.TOTALnontarget = (SNP.TRUTH.TP) / 
                                          SNP.ntTRUTH.TOTAL)

  left_join(qv_metrics, recall_metrics)
}

#' Make assembly benchmarking metrics table Calculates assembly benchmarking
#' results metrics table from targeted and non-targeted hap.py benchmarking
#' results table
#'
#' @param targeted_ext_dir path to directory with targeted benchmarking extended.csv 
#' @param nontargeted_ext_dir path to directory with nontargeted benchmarking extended.csv 
#' @param strats list of stratifications to subset for metrics table
#'
#' @return tibble with assembly benchmarking metrics
#' @export
#'
#' @examples make_asm_metrics_table("results/happy/nontargeted", "results/happy/targeted")
make_asm_metrics_table <- function(targeted_ext_dir, nontargeted_ext_dir, strats = strats){
  ## Loading targeted and non-targeted files
  targeted_results <- read_targeted_ext_results(targeted_ext_dir)
  nontargeted_results <- read_nontargeted_ext_results(nontargeted_ext_dir)
  
  ## Defining stratifications used in subsetting
  ### TODO - modify to allow for benchmarking GRCh37 or GRCh38
  strats <- c("*" ,
              "GRCh38_HG002_GIABv4.1_notin_complexandSVs_alldifficultregions.bed.gz" ,
              "GRCh38_notinalldifficultregions.bed.gz",
              "GRCh38_notinalllowmapandsegdupregions.bed.gz",
              "GRCh38_notinsegdups.bed.gz",
              "GRCh38_notinAllTandemRepeatsandHomopolymers_slop5.bed.gz",
              "GRCh38_notinAllHomopolymers_gt6bp_imperfectgt10bp_slop5.bed.gz",
              "GRCh38_segdups.bed.gz",
              "GRCh38_AllHomopolymers_gt6bp_imperfectgt10bp_slop5.bed.gz",
              "GRCh38_MHC.bed.gz")
  
  ## Combining and reformatting results
  combined_results <- left_join(targeted_results, 
                                nontargeted_results) %>%
    reformat_results(strats)
  
  
  ## Calculating Metrics
  benchmarking_summary <- calculate_asm_metrics(combined_results)
  
  benchmarking_summary
}



################################################################################
## Saving Metrics Table ########################################################
################################################################################
# write.csv(benchmarking_summary, outfile)
