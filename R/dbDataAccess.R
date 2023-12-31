##############
#    Database access
##############

#' @title Get target list from a single source
#'
#' @description This function queries precompiled annotation SQLite database
#' which contains miRNA - target gene associations with their respective scores.
#' @export
#' @param mirna miRNA in a standard format
#' @param source a source target prediction algorithm table to query, default 
#' \code{'diana'}, other possible values are \code{'miranda'}, 
#' \code{'targetscan'}, and \code{'pictar'}.
#' @param species species in a standard three-letter acronym, default 
#' \code{'mmu'}
#' @param synonyms when searching for -3p miRNA automatically also searches for
#' miRNA with the same name but ending with * (some databases list -3p miRNA
#' this way) and other way around, similarly for -5p miRNA, default TRUE
#' @param both_strands overrides \code{synonyms} and searches for targets of
#' both -5p and -3p strands together
#' @return \code{data.frame} object with entrez IDs of target genes and their
#' scores, if there are no targets found for a given miRNA in a given 
#' table then an empty 
#' @author Maciej Pajak \email{m.pajak@@sms.ed.ac.uk}
#' @references
#' Friedman, R. C., Farh, K. K.-H., Burge, C. B., and Bartel, D. P. (2009). 
#' Most mammalian mRNAs are conserved targets of microRNAs. Genome research, 
#' 19(1):92-105.
#' @references
#' Griffiths-Jones, S., Saini, H. K., van Dongen, S., and Enright, A. J. 
#' (2008). miRBase: tools for microRNA genomics. Nucleic acids research, 
#' 36(Database issue):D154-8.
#' @references
#' Lall, S., Grun, D., Krek, A., Chen, K., Wang, Y.-L., Dewey, C. N., ... 
#' Rajewsky, N. (2006). A genome-wide map of conserved microRNA targets in 
#' C. elegans. Current biology : CB, 16(5):460-71.
#' @references
#' Maragkakis, M., Vergoulis, T., Alexiou, P., Reczko, M., Plomaritou, K., 
#' Gousis, M., ... Hatzigeorgiou, A. G. (2011). DIANA-microT Web server 
#' upgrade supports Fly and Worm miRNA target prediction and bibliographic 
#' miRNA to disease association. Nucleic Acids Research, 39(Web Server issue), 
#' W145-8.
#' @examples
#' targets <- getTargetsFromSource('let-7a', species='hsa', source='targetscan')
#' head(targets) 
#' #top of the listof human targets of let-7a from TargetScan only
getTargetsFromSource <- function(mirna, species = 'mmu', source = 'diana',
                                 synonyms = TRUE, both_strands = FALSE) {
    tryCatch( {
        require(miRNAtap.db)
        
        cols <- c('GeneID','score')
        tabname <- paste(source,'_mapped_',species,sep='')
        tabnameUC <- toupper(tabname)
        
        #res - unique names of mirna of rbind of results of all synonyms!
        #if both strands then add wildcard to keys=mirna
        
        if (synonyms) {
            syn_mirnas <- .mirnaSynonyms(mirna)
        } else {
            syn_mirnas <- mirna
        }
         
        if (both_strands) { # sql wildcard at the end of base name
          syn_mirnas <- paste0(.justBaseName(mirna),'%')
        }
        
        if (tabnameUC %in% keytypes(miRNAtap.db)) {
        res <- Reduce(rbind, lapply(syn_mirnas, 
                function(x)  
                    select(miRNAtap.db,keys=x,columns=NULL,
                                                keytype=tabnameUC)[,cols]
          ))
           
        } else if (species=='rno') {
        tabname <- paste(source,'_mapped_mmu', sep='')
        tabnameUC <- toupper(tabname)
        if (tabnameUC %in% keytypes(miRNAtap.db)) {
          res <- Reduce(rbind, lapply(syn_mirnas, 
                                      function(x)  
                                        select(miRNAtap.db,keys=x,columns=NULL,
                                               keytype=tabnameUC)
          ))  
            res <- translate(res,from='mmu',to=species)[,cols]
        } else {
            return(NULL)
        }
        } else {
        return(NULL)
        }

        res <- res[order(-res$score),]
        res <- res[!duplicated(res$GeneID),]
        
        return(res)
    } , error = function(cond) {
        message(cond)
        return(NULL)
    }, warning = function(cond) {
        message(cond)
        return(NULL)
    })
    
    
}




#####################
#     homology translation
#####################


#' @title Homology transfer for miRNAtap
#'
#' @description 
#' This function maps gene entrez ID between species using homology information 
#' from Homologene.
#'
#' @param entrezes data.frame with entrez Gene IDs and their scores
#' @param from origin species, default \code{'mmu'}, Mus musculus
#' @param to target species, default
#' @param ... any optional arguments
#' @export
#' @return data.frame object with orthologous genes' entrez IDs and 
#' corresponding scores
#' @author Maciej Pajak \email{m.pajak@@sms.ed.ac.uk}    
#' @examples
#' mouse_genes <- data.frame(GeneID = 
#'         c("15364", "56520", "57781", "58180", "18035", "239857"))
#' translate(mouse_genes, from='mmu', to='rno')
translate <- function(entrezes, from = 'mmu', to = 'rno', ...) {
    tryCatch( {
        require(miRNAtap.db)
        
        conL <- dbconn(miRNAtap.db)
        species_map <- data.frame(code=c('mmu','rno','hsa'),
                    tax_id=c('10090','10116','9606'), stringsAsFactors=FALSE)

        from.id <- species_map[species_map$code==from,'tax_id']
        to.id <- species_map[species_map$code==to,'tax_id']

        one_string <- paste("'", paste(entrezes$GeneID,collapse="','"),
                            "'", sep='')

        query <- paste(
    "SELECT h1.family AS family, h1.entrez AS GeneID, h2.entrez AS to_entrez 
        FROM homologene_raw AS h1, homologene_raw AS h2 
        WHERE h1.tax_id=",from.id," AND h1.entrez IN (",one_string,") 
        AND h2.tax_id=",to.id," AND h1.family=h2.family;",
            sep="")
                
                
        mappings <- dbGetQuery(conL, query)
        entrezes <- join(entrezes,mappings,by='GeneID')

        if ('mirna' %in% colnames(entrezes) & 'score' %in% colnames(entrezes)){
        entrezes <- entrezes[,c('mirna','score','to_entrez')]
        colnames(entrezes) <- c('mirna','score','GeneID')
        } else {
        colnames(entrezes) <- c(paste('GeneID_',from,sep=''),'family','GeneID')
        }
    
    #dbDisconnect(conL)
    return(entrezes) #temp
    } , error = function(cond) {
            message(cond)
            return(NULL)
    }, warning = function(cond) {
            message(cond)
            return(NULL)
    } )
    
    
}

