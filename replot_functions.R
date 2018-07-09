library(magrittr)
library(data.table)
library(gtable)
library(grid)
library(gridExtra)
library(ggplot2)
library(cowplot)
library(seqsetvis)

#' loads data from GSEA output for replotting
#'
#' @param rnk_file path to rank file
#' @param sets_file path to sets file
#' @param results_file path to results file
#'
#' @return list of rank, sets, and results ready for plotting
#' @export 
#'
#' @examples
load_gsea = function(rnk_file, sets_file, results_file){
  #load rank
  rnk_df = read.table(rnk_file, stringsAsFactors = F)
  rnk = rnk_df[order(rnk_df$V2), ]$V1
  
  #load sets
  sets_raw = read.table(sets_file, sep= "\n", stringsAsFactors = F)$V1
  sets = lapply(strsplit(sets_raw, "\t"), function(x){
    x[-(1:2)]
  })
  names(sets) = sapply(strsplit(sets_raw, "\t"), function(x){
    x[1]
  })
  
  #load scores
  scores_raw = read.table(results_file, sep = "\n", stringsAsFactors = F)$V1
  scores_raw = scores_raw[grepl("^  <DTG", scores_raw)]
  
  GENESET = scores_raw %>% sub(".+GENESET=", "", .) %>% sub(" ES=.+", "", .) %>% sub(".+#", "", .)
  ES = scores_raw %>% sub(".+ ES=", "", .) %>% sub(" NES=.+", "", .)
  NES = scores_raw %>% sub(".+ NES=", "", .) %>% sub(" NP=.+", "", .)
  NP = scores_raw %>% sub(".+ NP=", "", .) %>% sub(" FDR=.+", "", .)
  FDR = scores_raw %>% sub(".+ FDR=", "", .) %>% sub(" FWER=.+", "", .)
  FWER = scores_raw %>% sub(".+ FWER=", "", .) %>% sub(" RND_ES=.+", "", .)
  RND_ES = scores_raw %>% sub(".+ RND_ES=", "", .) %>% sub(" HIT_INDICES=.+", "", .)
  HIT_INDICES = scores_raw %>% sub(".+ HIT_INDICES=", "", .) %>% sub(" ES_PROFILE=.+", "", .)
  ES_PROFILE = scores_raw %>% sub(".+ ES_PROFILE=", "", .) %>% sub(" RANK_AT_ES=.+", "", .)
  RANK_AT_ES = scores_raw %>% sub(".+ RANK_AT_ES=", "", .) %>% sub(" RANK_SCORE_AT_ES=.+", "", .)
  RANK_SCORE_AT_ES = scores_raw %>% sub(".+ RANK_SCORE_AT_ES=", "", .) %>% sub("/>", "", .)
  
  scores_dt = data.table(GENESET, ES, NES, NP, FDR, FWER, RND_ES, HIT_INDICES, ES_PROFILE, RANK_AT_ES, RANK_SCORE_AT_ES)
  
  list(rnks = rnk, sets = sets, scores = scores_dt)
}

#' Combines loaded GSEA data lists so they can be plotted simultaneously
#'
#' @param gsea1 loaded set 1
#' @param gsea2 loaded set 2
#'
#' @return combination of set 1 and set 2
#' @export
#'
#' @examples
combine_gsea = function(gsea1, gsea2){
  if(!all(gsea1$rnks == gsea2$rnks)){
    stop("ranks don't match and should!")
  }
  out = list()
  out$rnks = gsea1$rnks
  out$sets = append(gsea1$sets, gsea2$sets)
  out$scores = rbind(gsea1$scores, gsea2$scores)
  return(out)
}

#' Extracts data necessary for plotting for set_name from gsea_data
#'
#' @param set_name perfect match GSEA set name
#' @param gsea_data loaded GSEA data
#'
#' @return data.table with profile and hit_index data for plotting
#' @export
#'
#' @examples
extract_set = function(set_name, gsea_data){
  soi_dt = gsea_data$scores[GENESET == set_name]
  prof = as.numeric(strsplit(soi_dt$ES_PROFILE, " ")[[1]])
  # plot(prof)
  hi = as.numeric(strsplit(soi_dt$HIT_INDICES, " ")[[1]])
  # plot(hi, rep(1, length(hi)))
  rbind(
    data.table(x = hi, y = prof, set = set_name, type = "profile"),
    data.table(x = hi, y = 0, set = set_name, type = "hit_index")
  )
}