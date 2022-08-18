#
# PWI.R
# Author: Dr. Robin Haunschild
# Version: 0.0.1
# Date: 07/29/2022
#

#' @title Function to calculate prize winner indices based on bibliometric data
#'
#' @description
#' This function calculates prize winner indices based on bibliometric data.
#' The default prize winners are the recipients of the Derek de Solla Price Memorial Medal.
#' Users can provide recipients of other prizes.
#'
#' @details
#' PWI(files=bibliographic_files, pw_pattern = '(BORNMANN L)|(BAR-ILAN J)|(WALTMAN L)|
#' (THELWALL M)|(CRONIN B)|(PERSSON O)|(VINKLER P)|(MCCAIN K)|(INGWERSEN P)|
#' (LEYDESDORFF L)|(ROUSSEAU R)|(GLANZEL W)|(GLAENZEL W)|(MOED H)|(IRVINE J)|(MARTIN B)
#' |(GRIFFITH B)|(VAN RAAN A)|(VANRAAN A)|(MERTON R)|(SCHUBERT A)|(BROOKES B)|
#' (NARIN F)|(NALIMOV V)|(BRAUN T)|(MORAVCSIK M)|(GARFIELD E)',
#' method=1, verbosity=1)
#' Only the argument files is necessary. All other arguments are optional.
#'
#' Literature:
#'
#' - Bornmann, L. & Haunschild, R. (in preparation). "Prize Winner Index".
#'
#' @examples
#' \donttest{
#' JoI <- PWI('http://andreas-thor.github.io/cre/data/savedrecs_JOI2.txt')
#' head(JoI)
#' }
#' @param files character variable or list of character variables that contain(s)
#' file name(s) of bibliographic data file(s) that are supported by the package
#' \link{bibliometrix}
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
#' @importFrom bibliometrix convert2df
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

PWI <- function(files, pw_pattern = '(BORNMANN L)|(...',
                method = 1, verbosity=1) {
# requireNamespace('bibliometrix')

# attach(countries)
# on.exit('detach(countries)')

 if(pw_pattern == '(BORNMANN L)|(...') {
   pw_pattern <- '(BORNMANN L)|(BAR-ILAN J)|(WALTMAN L)|(THELWALL M)|(CRONIN B)|(PERSSON O)|(VINKLER P)|(MCCAIN K)|(INGWERSEN P)|(LEYDESDORFF L)|(ROUSSEAU R)|(GLANZEL W)|(GLAENZEL W)|(MOED H)|(IRVINE J)|(MARTIN B)|(GRIFFITH B)|(VAN RAAN A)|(VANRAAN A)|(MERTON R)|(SCHUBERT A)|(BROOKES B)|(NARIN F)|(NALIMOV V)|(BRAUN T)|(MORAVCSIK M)|(GARFIELD E)'
 }
 df <- bibliometrix::convert2df(files)
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
  df_dist_freq_weight_sum_agg <- df_dist_freq_weight_sum_agg[order(-df_dist_freq_weight_sum_agg$PWI),]
  return(df_dist_freq_weight_sum_agg)
 }
}

