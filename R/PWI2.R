#
# PWI2.R
# Author: Dr. Robin Haunschild
# Version: 0.0.3
# Date: 10/20/2023
#

#' @title Function to calculate prize winner indices based on bibliometric data
#'
#' @description
#' This function calculates prize winner indices based on bibliometric data.
#' The default prize winners are the recipients of the Derek de Solla Price Memorial Medal.
#' Users can provide recipients of other prizes.
#'
#' @details
#' PWI2(bib_df=bibliographic_dataframe, pw_pattern = '(BOYACK K)|(KLAVANS R)|(BORNMANN L)|(BAR-ILAN J)|(BARILAN J)|
#' (WALTMAN L)|(THELWALL M)|(CRONIN B)|(PERSSON O)|(VINKLER P)|(MCCAIN K)|(INGWERSEN P)|
#' (LEYDESDORFF L)|(ROUSSEAU R)|(EGGHE L)|(GLANZEL W)|(GLAENZEL W)|(MOED H)|(IRVINE J)|
#' (MARTIN B)|(GRIFFITH B)|(VAN RAAN A)|(VANRAAN A)|(MERTON R)|(SCHUBERT A)|(BROOKES B)|
#' (NARIN F)|(NALIMOV V)|(BRAUN T)|(MORAVCSIK M)|(GARFIELD E)',
#' method=1, verbosity=1)
#' Only the argument bib_df is necessary. All other arguments are optional.
#'
#' Literature:
#'
#' - Bornmann, L. & Haunschild, R. (in preparation). "Prize Winner Index".
#'
#' @examples
#' \donttest{
#' bib_df <- bibliometrix::convert2df('http://andreas-thor.github.io/cre/data/savedrecs_JOI2.txt')
#' JoI <- PWI2(bib_df)
#' head(JoI)
#' }
#' @param bib_df bibliographic dataframe variable from \link{convert2df}
#' @param pw_pattern character variable (optional parameter) that is passed as
#' search pattern to the \link{grep} function to identify the prize winners in
#' the data set
#' @param method integer variable (optional parameter) that determines if only
#' the authors in the dataset with number of papers and co-authors is returned
#' or if the prize winner index is calculated
#' 0: return only a list with authors, number of papers, and number of co-authorships
#' 1: calculate the prize winner index and return it alongside with number of papers
#' and number of co-authorships
#' @param verbosity level of verbosity (0=quiet and 1=informative)
#'
#' @returns data frame of researcher names, PWI value, number of papers, and number of co-authors
#'
#' @export
#'
#' @importFrom bibliometrix biblioNetwork
#' @importFrom bibliometrix biblioAnalysis
#' @importFrom igraph graph.adjacency
#' @importFrom igraph as_edgelist
#' @importFrom igraph vertex_attr
#' @importFrom igraph get.shortest.paths
#' @importFrom igraph V
#' @importFrom stats aggregate
#' @importFrom utils setTxtProgressBar
#' @importFrom utils txtProgressBar

PWI2 <- function(bib_df, pw_pattern = '(BOYACK K)|(KLAVANS R)|(BORNMANN L)|(...',
                method = 1, verbosity=1) {

 if(pw_pattern == '(BOYACK K)|(KLAVANS R)|(BORNMANN L)|(...') {
   pw_pattern <- '(BOYACK K)|(KLAVANS R)|(BORNMANN L)|(BAR-ILAN J)|(BARILAN J)|(WALTMAN L)|(THELWALL M)|(CRONIN B)|(PERSSON O)|(VINKLER P)|(MCCAIN K)|(INGWERSEN P)|(LEYDESDORFF L)|(ROUSSEAU R)|(EGGHE L)|(GLANZEL W)|(GLAENZEL W)|(MOED H)|(IRVINE J)|(MARTIN B)|(GRIFFITH B)|(VAN RAAN A)|(VANRAAN A)|(MERTON R)|(SCHUBERT A)|(BROOKES B)|(NARIN F)|(NALIMOV V)|(BRAUN T)|(MORAVCSIK M)|(GARFIELD E)'
 }
 df <- bib_df
 if(verbosity==1) { print(paste('Read', length(unique(df$UT)), 'unique documents.')) }
 NetMat <- bibliometrix::biblioNetwork(df, analysis = 'collaboration', network = 'authors')

# requireNamespace('igraph')
 bsk.network <- igraph::graph.adjacency(NetMat,mode = "undirected")

 v <- igraph::V(bsk.network)
 bib_anal <- bibliometrix::biblioAnalysis(df)
 au <- bib_anal$Authors
 df_au <- as.data.frame(au)
 bsk.network_d <- igraph::graph.adjacency(NetMat,mode = "directed")
 el <- unique(as.data.frame(igraph::as_edgelist(bsk.network_d)))
 el$cnt <- 1
 el_agg <- aggregate(cnt ~ V1, data = el, FUN = sum)
 df_au <- merge(df_au, el_agg, by.x='AU', by.y='V1')
 colnames(df_au) <- c('NAME', 'PAPERS', 'COAUTHORS')
 df_au$COAUTHORS <- df_au$COAUTHORS-1
 if(method == 0) {
   df_au <- df_au[order(-df_au$PAPERS),]
   df_au$PW <- 'default'
   df_au[grepl(pw_pattern, df_au$NAME),]$PW <- 'YES'
   df_au[!grepl(pw_pattern, df_au$NAME),]$PW <- 'NO'
   return(df_au)
 } else {
  pw_ids <- grep(pw_pattern, v$name)
  df_dist <- data.frame(node = 'name', dist = 10, freq = 0)

  if(verbosity==1) { pb <- txtProgressBar(0, length(vertex_attr(bsk.network)$name), style=3) }
  for(i in seq(1, length(vertex_attr(bsk.network)$name))) {
    if(verbosity==1) { setTxtProgressBar(pb, i) }
    for(j in pw_ids) {
      p <- igraph::get.shortest.paths(bsk.network, j, i)
      p_occ <- data.frame(link = 'dummy', strength = 0)
      if(length(p$vpath[[1]])>1) {
        for(k in seq(1, length(p$vpath[[1]])-1)) {
          tmp <- data.frame(link = paste(p$vpath[[1]][k]$name, p$vpath[[1]][k+1]$name, sep=' -- '), strength = NetMat[rownames(NetMat) %in% p$vpath[[1]][k]$name,colnames(NetMat) %in% p$vpath[[1]][k+1]$name])
          p_occ <- rbind(p_occ, tmp)
        }
        p_occ <- p_occ[p_occ$link != 'dummy',]
      } else {
        if(length(p$vpath[[1]])>0) {
          p_occ <- data.frame(link = paste(p$vpath[[1]][1]$name, p$vpath[[1]][1]$name, sep=' -- '), strength = NetMat[rownames(NetMat) %in% p$vpath[[1]][1]$name,colnames(NetMat) %in% p$vpath[[1]][1]$name])
        } else {
          p_occ <- data.frame(link = paste(V(bsk.network)[i], V(bsk.network)[j], sep=' -- '), strength = 0)
        }
      }
      tmp <- data.frame(node = igraph::vertex_attr(bsk.network)$name[i], dist = length(V(bsk.network)[p[1]$vpath[[1]]]), freq = min(p_occ$strength))
      df_dist <- rbind(df_dist, tmp)
    }
  }

  if(verbosity==1) { close(pb) }

  df_dist <- df_dist[df_dist$node != 'name',]
  df_dist$weight <- 1/(2^(df_dist$dist-1))
  df_dist_freq_weight_sum_agg <- aggregate(weight*freq ~ node, data = df_dist, FUN = sum)
  colnames(df_dist_freq_weight_sum_agg) <- c('NAME', 'PWI')
  df_dist_freq_weight_sum_agg <- merge(df_dist_freq_weight_sum_agg, df_au, by='NAME')

# Add minimum and maximum publication years
  df_py <- data.frame(NAME = 'dummy', MIN_PY=-1, MAX_PY=-2)
  for(i in df_dist_freq_weight_sum_agg$NAME) {
     min_py <- min(df[grep(i, df$AU),]$PY, rm.NA=TRUE)
     max_py <- max(df[grep(i, df$AU),]$PY, rm.NA=TRUE)
     df_tmp <- data.frame(NAME=i, MIN_PY = min_py, MAX_PY = max_py)
     df_py <- rbind(df_py, df_tmp)
  }
    
  df_dist_freq_weight_sum_agg <- merge(df_dist_freq_weight_sum_agg, df_py, by='NAME')
  df_dist_freq_weight_sum_agg <- df_dist_freq_weight_sum_agg[order(-df_dist_freq_weight_sum_agg$PWI),]
  df_dist_freq_weight_sum_agg$PW <- 'default'
  df_dist_freq_weight_sum_agg[grepl(pw_pattern, df_dist_freq_weight_sum_agg$NAME),]$PW <- 'YES'
  df_dist_freq_weight_sum_agg[!grepl(pw_pattern, df_dist_freq_weight_sum_agg$NAME),]$PW <- 'NO'
  return(df_dist_freq_weight_sum_agg)
 }
}

