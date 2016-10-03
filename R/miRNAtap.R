
#' @docType package
#' @name miRNAtap
#' @title miRNAtap: microRNA Targets - Aggregated Predictions.
#' @description It is a package with tools to facilitate implementation of
#' workflows
#' requiring miRNA prediction through access to multiple prediction results 
#' (DIANA, Targetscan, PicTar, Miranda, and miRDB) and their aggregation. 
#' Three aggregation methods are available: minimum, maximum and geometric mean,
#' additional parameters provide further tuning of the results. 
#' Predictions are available for Homo sapiens, Mus musculus 
#' and Rattus norvegicus (the last one through homology translation).
#' @import AnnotationDbi RSQLite DBI stringr sqldf plyr methods
#' @author Maciej Pajak \email{m.pajak@@sms.ed.ac.uk}, Ian Simpson 
#' @examples
#' #direct targets in mouse aggregated from all sources:
#' targets_mouse <- getPredictedTargets('let-7a',species='mmu', method='geom') 
#' #homology-translated targets in rat aggregated from all sources
#' targets_rat <- getPredictedTargets('let-7a',species='rno', method='geom') 
NULL


####################
####################
# ANALYTICS
####################
####################



#####################
# MAIN FUNCTION
#####################

#' @title Get aggregated ordered list of predicted targets for miRNA
#'
#' @description This method performs aggregation of target lists from multiple 
#' sources. Aggregated list is more accurate than any list from a single 
#' source. Multiple aggregation methods are available.Direct target data from 
#' five sources for Human and Mouse is supplied through \code{miRNAtap.db} 
#' package, for Rat targets are derived through homology translations whenever 
#' direct ones are not available.
#'
#' @usage getPredictedTargets(mirna, sources = c("pictar", "diana", 
#' "targetscan", "miranda","mirdb"), species = "mmu", min_src = 2, 
#' method = "geom", promote = TRUE, synonyms = TRUE, both_strands = FALSE, ...)
#' @param mirna miRNA in a standard format
#' @param sources a list of sources to use for aggregation, 
#' default is all five sources, i.e. 
#' \code{c('pictar','diana','targetscan','miranda','mirdb')}
#' @param species species in a standard three-letter acronym, \code{'mmu'} 
#' and \code{'hsa'} 
#' available as direct targets, \code{'rno'} as homology translations, 
#' default \code{'mmu'}
#' @param min_src minimum number of sources required for a target to be 
#' considered, default 2
#' @param method method of aggregation - choose from \code{'min'}, 
#' \code{'max'}, and \code{'geom'}; 
#' \code{'min'} is a minimum of ranks, \code{'max'} is a maximum of ranks,
#' and default \code{'geom'} 
#' is based on geometric mean of the ranks which proves to be the most 
#' accurate method.
#' @param promote add weights to improve accuracy of the method, default TRUE
#' @param synonyms when searching for -3p miRNA automatically also searches for
#' miRNA with the same name but ending with * (some databases list -3p miRNA
#' this way) and other way around, similarly for -5p miRNA, default TRUE
#' @param both_strands overrides \code{synonyms} and searches for targets of
#' both -5p and -3p strands together
#' @param ... any optional arguments
#' @return \code{data.frame} object where row names are entrez IDs of target
#' genes, ranks from individual sources and aggregated rank are shown 
#' in columns.
#' If no targets are found in any of the sources \code{NULL} and a warning
#' are returned.
#' @details Tuning \code{min_src} parameter is an easy way of prioritising 
#' precision at the top of the list (high values) or total recall (low values).
#' For the five default input sources, recommended values are 2, 3, or 4.
#' @export
#' @author Maciej Pajak \email{m.pajak@@sms.ed.ac.uk}
#' @references
#' Agarwal V, Bell GW, Nam J, Bartel DP. Predicting effective microRNA 
#' target sites in mammalian mRNAs. eLife, 4:e05005, (2015).
#' @references
#' Griffiths-Jones, S., Saini, H. K., van Dongen, S., and Enright, A. J. 
#' (2008). miRBase: tools for microRNA genomics. Nucleic acids research, 
#' 36(Database issue):D154-8.
#' @references
#' Lall, S., Grun, D., Krek, A., Chen, K., Wang, Y.-L., Dewey, C. N., ... 
#' Rajewsky, N. (2006). A genome-wide map of conserved microRNA targets in 
#' C. elegans. Current biology : CB, 16(5):460-71.
#' @references
#' Paraskevopoulou MD, Georgakilas G, Kostoulas N, Vlachos IS, Vergoulis 
#' T, Reczko M, Filippidis C, Dalamagas T, Hatzigeorgiou AG., 
#' "DIANA-microT web server v5.0: service integration into miRNA 
#' functional analysis workflows.", Nucleic Acids Res. 2013 
#' Jul;41(Web Server issue):W169-73.
#' @references 
#' Wong N and Wang X (2015) miRDB: an online resource for microRNA 
#' target prediction and functional annotations. Nucleic Acids Research. 
#' 43(D1):D146-152.
#' @examples
#' targets <- getPredictedTargets('let-7a',species='hsa', method = 'min') 
#' head(targets) #top of the list with minimum aggregation
#' targets2 <- getPredictedTargets('let-7a',species='hsa', method='geom') 
#' head(targets2) #top of the list with geometric mean aggregation
getPredictedTargets <- function(mirna,
            sources=c('pictar','diana','targetscan','miranda','mirdb'), 
            species = 'mmu', min_src = 2, method = 'geom',
            promote = TRUE, synonyms = TRUE, both_strands = FALSE, ...) {
    
    if (!(species %in% c('mmu','hsa','rno','dme'))) {
        warning(paste('species ',species,
                ' not supported in the current version',sep=''))
        return(NULL)
    }
    n_sources <- length(sources)
    
    if (n_sources<min_src) {
        warning('min_src > n_sources!')
        return(NULL)
    }
    
    if (min_src<=0) {
        warning('min_src changed to 1')
        min_src <- 1
    }
    
    if (nchar(mirna)<3) {
        warning('unrecognised miRNA')
        return(NULL)
    }
    
    l_outputs <- list()
    
    cols <- c('GeneID','score')
    
    if (substr(mirna,1,3)==species) {
        mirna=substr(mirna,5,15)
    }
    
    for (src in sources) {
        l_outputs[[src]] <- getTargetsFromSource(mirna, species, source = src,
                                                 synonyms = synonyms,
                                                 both_strands = both_strands)
    }
    
    #it creates non-unique colnames, hence warning surpression
    merged.scores <- suppressWarnings(Reduce(
                function(...) .mergeRename(..., by='GeneID', all=TRUE), 
                l_outputs))
    
    if (is.null(merged.scores)) { # | dim(merged.scores)[1]<1 
        warning(paste('no targets found for mirna ', mirna, sep = ''))
        return(NULL)
    }
    
    valid_srcs <- (1:(n_sources+1))[colSums(merged.scores,na.rm = TRUE)>0]
    
    if (length(valid_srcs)<2) { # less than 2 cos the first column isnt src 
        warning(paste('no targets found for mirna ', mirna, sep = ''))
        return(NULL)
    }
    
    n_valid_srcs <- length(valid_srcs)-1
    
    if (n_valid_srcs<min_src) {
        warning(paste(
        'sources which returned a target list<min_src, min_src reduced to ',
        n_valid_srcs,sep=''))
        min_src <- n_valid_srcs
    }
    
    merged.scores <- merged.scores[,valid_srcs]
    
    #take over by ranking function.
    result <- aggregateRanks(merged.scores, 
                n_valid_srcs, min_src, method = method, promote=promote)
    
    return(result)
}


####################
# Ranking function
####################
#' @title Aggreagate ranks from multiple sources with various methods
#'
#' @description This function performs aggregation phase of target 
#' prediction for \code{\link{getPredictedTargets}}. 
#' Consensus ranking is derived from multiple individual rankings. 
#' Available methods include minimum, maximum and geometric mean with further 
#' tuning parameters which promote true positives at the top of the final 
#' ranking 
#' @export
#' @param ranks \code{data.frame} with ordered scores
#' @param n_valid_srcs number of valid sources in the dataset
#' @param min_src minimum acceptable number fo sources
#' @param method \code{'min'},\code{'max'}, or \code{'geom'},
#' default \code{'geom'}
#' @param promote add weights to improve accuracy of the method, default 
#' \code{TRUE}
#' @return \code{data.frame} object with ranks per source and aggregate ranks
#' @author Maciej Pajak \email{m.pajak@@sms.ed.ac.uk}
#' @examples
#' data = data.frame(GeneID=c("15364", "56520", "57781", "58180", "18035"),
#'                 source1scores=c(0.9,0.5,0.3,NA,NA),
#'                 source2scores=c(0.7,NA,0.8,0.6,0.5),
#'                 source3scores=c(0.5,NA,0.3,0.1,0.2))
#' data #dataframe with scores
#' aggregateRanks(data, n_valid_srcs=3, min_src=2, method='geom')
#' #note how gene 56520 is eliminated as it appeared in fewer than 2 sources
aggregateRanks <- function(ranks, n_valid_srcs, min_src,
                            method = 'geom', promote=TRUE) {

    if (n_valid_srcs>1) {
        target_found <- rowSums(!is.na(ranks[,2:(n_valid_srcs+1)])) #n_sources
    } else {
        target_found <- 1*(!is.na(ranks[,2:(n_valid_srcs+1)])) #n_sources
    }
    
    ranks <- ranks[target_found>=min_src,]
    ranks <- ranks[!duplicated(ranks$GeneID),]
    ranks <- ranks[!is.na(ranks$GeneID),]
    
    row.names(ranks) <- ranks$GeneID
    
    if (n_valid_srcs==1) {
        data1 <- data.frame(source1scores=ranks[,2], row.names = row.names(ranks))
    } else {
        data1 <- ranks[,2:(n_valid_srcs+1)] #n_sources
        data1 <- as.data.frame(data1)
    }
    num.gene <- dim(data1)[1]
    
    if (num.gene<1) {
        warning('not enough targets overlapping between sources,
reduce min_srcs parameter or add sources')
        return(NULL)
    } else if (num.gene==1) {
        warning('only one target found')
        result <- data.frame(matrix(1,nrow=1,ncol=n_valid_srcs+2))
        colnames(result) <- c(paste('source_',1:n_valid_srcs,sep=''),
                            'rank_product','rank_final')
        row.names(result) <- row.names(data1)
        return(result)
    }
    
    #we need to reverse rank them
    rank_ind <- apply(data1, 2, 
    function(x) (num.gene+1) - rank(x, ties.method='average', na.last = FALSE))
    
    
    if (method=='geom') {
        result <- .aggregateGeom(data1, rank_ind, promote)
    }  else if (method=='max') {
        result <- .aggregateMinMax(data1, rank_ind, minmax=method)
    } else { #even if it's something invalid, still use min
        result <- .aggregateMinMax(data1, rank_ind, minmax=method)
    } #more methods possible to add, so far geom suffices and scores well
    
    colnames(result) <- c(paste('source_',1:n_valid_srcs,sep=''),
                            'rank_product','rank_final') #coating
    return(result)
}



#####################
#  Auxiliary
#####################

.aggregateGeom <- function(data1, rank_ind, promote=TRUE) {
    rank_ind[is.na(data1)] <- 1
    num.rank <- apply(is.na(data1) == FALSE, 1, sum)
    if (promote==TRUE) {
        rank_mean <- (1/num.rank)*(apply(rank_ind, 1, prod))^(1/num.rank)
    } else {
        rank_mean <- (apply(rank_ind, 1, prod))^(1/num.rank)
    }
    rank_mean[num.rank == 0] <- NA
    rank_final <- rank(rank_mean)
    rank_ind[is.na(data1)] <- NA
    result <- cbind(rank_ind,rank_mean,rank_final)  
    result <- result[order(result[,'rank_final']),]
    return(result)
}

.aggregateMinMax <- function(data1, rank_ind, minmax='min') {
    rank_ind[is.na(data1)] <- NA
    if (minmax=='max') {
        rank_mean <- apply(rank_ind, 1, max, na.rm=TRUE)
    } else {
        rank_mean <- apply(rank_ind, 1, min, na.rm=TRUE)
    }
    rank_final <- rank(rank_mean)
    result <- cbind(rank_ind,rank_mean,rank_final)
    result <- result[order(result[,'rank_final']),]
}


.mergeRename <- function(...) {
    merged <- merge(...)
    colnames(merged) <- c(colnames(merged)[1],
                        paste('V',1:(length(names(merged))-1),sep=''))
    return(merged)
}

.mirnaSynonyms <- function(mirna) {
    if (str_sub(mirna,start=-1)=='*') {
        return(c(mirna, paste0(str_sub(mirna, start=1, end=-2),'-3p')))
    } else if (str_sub(mirna,start=-3)=='-3p') {
        return(c(mirna, paste0(str_sub(mirna, start=1, end=-4),'*')))
    } else if (str_sub(mirna,start=-3)=='-5p') {
        return(c(mirna, str_sub(mirna, start=1, end=-4)))
    } else  {
        return(c(mirna, paste0(mirna,'-5p')))
    }
}

.allBothStrands <- function(mirna) {
    base <- str_match(mirna, '[a-zA-Z]+-[0-9a-z]+')[1,] 
    if (!is.na(base)) {
       return(paste0(base,c('','-5p','*','-3p')))
    } else {
      return(mirna)
    }
  
}

.justBaseName <- function(mirna) {
  base <- str_match(mirna, '[a-zA-Z\\-]+-[0-9a-z]+')[1,] 
  if (!is.na(base)) {
    return(base)
  } else {
    return(mirna)
  }
  
}



