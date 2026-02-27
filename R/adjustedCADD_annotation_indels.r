adjustedCADD.annotation.indels <- function(x, variant.scores = NULL, cores = 10, verbose = T, path.data, build = c("b37", "b38")){
  if(bedr::check.binary(x = "bedtools", verbose = F)==F) stop("'bedtools' is not available and need to be installed on the system")

  if("adjCADD" %in% colnames(x@snps)){
    warning("'adjCADD' already exists and will be replaced")
    x@snps <- x@snps[,-which(colnames(x@snps)=="adjCADD")]
  }
  if(missing(path.data)) stop("the directory 'path.data' to download and use the necessary files for RAVA-FIRST analysis should be provided")
  if(is.null(variant.scores))  stop("Annotation of indels, 'variant.scores' should be provided with Phred scores 1.4 for Indels")
  if(!(build %in% c("b37", "b38"))) {stop("Wrong build provided")}
  if(build == "b37"){
    if(!file.exists(paste0(path.data, "/AdjustedCADD_v1.4_202204_indels.tsv.gz"))){
        if(verbose){
          cat("Downloading adjusted CADD scores of indels in b37 in ", path.data, "\n")
          curl_download("https://lysine.univ-brest.fr/RAVA-FIRST/AdjustedCADD_v1.4_202204_indels.tsv.gz", destfile = paste0(path.data, "/AdjustedCADD_v1.4_202204_indels.tsv.gz"), quiet = F)
          curl_download("https://lysine.univ-brest.fr/RAVA-FIRST/AdjustedCADD_v1.4_202204_indels.tsv.gz.tbi", destfile = paste0(path.data, "/AdjustedCADD_v1.4_202204_indels.tsv.gz.tbi"))
        }else{
          curl_download("https://lysine.univ-brest.fr/RAVA-FIRST/AdjustedCADD_v1.4_202204_indels.tsv.gz", destfile = paste0(path.data, "/AdjustedCADD_v1.4_202204_indels.tsv.gz"))
          curl_download("https://lysine.univ-brest.fr/RAVA-FIRST/AdjustedCADD_v1.4_202204_indels.tsv.gz.tbi", destfile = paste0(path.data, "/AdjustedCADD_v1.4_202204_indels.tsv.gz.tbi"))
        }
    }
    CADDfile = paste0(path.data, "/AdjustedCADD_v1.4_202204_indels.tsv.gz")
    col.CADD = "PHRED_1.4"
  }
  if(build == "b38"){
    if(!file.exists(paste0(path.data, "/AdjustedCADD_v1.6_202602_indels.tsv.gz"))){
        if(verbose){
          cat("Downloading adjusted CADD scores of indels in b38 in ", path.data, "\n")
          curl_download("https://lysine.univ-brest.fr/RAVA-FIRST/AdjustedCADD_v1.6_202602_indels.tsv.gz", destfile = paste0(path.data, "/AdjustedCADD_v1.6_202602_indels.tsv.gz"), quiet = F)
          curl_download("https://lysine.univ-brest.fr/RAVA-FIRST/AdjustedCADD_v1.6_202602_indels.tsv.gz.tbi", destfile = paste0(path.data, "/AdjustedCADD_v1.6_202602_indels.tsv.gz.tbi"))
        }else{
          curl_download("https://lysine.univ-brest.fr/RAVA-FIRST/AdjustedCADD_v1.6_202602_indels.tsv.gz", destfile = paste0(path.data, "/AdjustedCADD_v1.6_202602_indels.tsv.gz"))
          curl_download("https://lysine.univ-brest.fr/RAVA-FIRST/AdjustedCADD_v1.6_202602_indels.tsv.gz.tbi", destfile = paste0(path.data, "/AdjustedCADD_v1.6_202602_indels.tsv.gz.tbi"))
        }
    }
    CADDfile = paste0(path.data, "/AdjustedCADD_v1.6_202602_indels.tsv.gz")
    col.CADD = "PHRED_1.6"
  }
    
    
    ##Annotation with adjusted CADD
    x.pos <- paste0(x@snps$chr, ":", format(x@snps$pos-1, scientific = F, trim = T), "-", format(x@snps$pos, scientific = F, trim = T))
    ##Sort and merge positions
    x.pos.sort <- bedr.sort.region(x.pos, check.chr=F, verbose = F, method="natural")
    x.pos.merged <- bedr.merge.region(x.pos.sort, check.chr=F, verbose = F)
    ##Split pos for parallelisation
    end.indexes <- round(seq(1, length(x.pos.merged), length.out = cores+1)[-1])
    start.indexes <- c(1, end.indexes[-length(end.indexes)]+1)
    ##Get CADD score for corresponding positions
    tmp.scores <- do.call(rbind, mclapply(1:cores, function(z) tabix(x.pos.merged[start.indexes[z]:end.indexes[z]], CADDfile, check.chr=F, verbose = F), mc.cores = cores))
    colnames(tmp.scores) <- c("chr", "pos", "A1", "A2", col.CADD, "adjCADD", "SubRegion")
    
    x@snps$order <- 1:nrow(x@snps) #To make sure that SNPs are in the right order
    scores.notflip.tmp <- merge(x@snps, tmp.scores[,c("chr", "pos", "A1", "A2", "adjCADD")], by = c("chr", "pos", "A1", "A2"), all.x = T)
    scores.notflip <- scores.notflip.tmp$adjCADD[order(scores.notflip.tmp$order)]
    #Inverse alleles to merge CADD scores when maf != @p
    scores.flip.tmp <- merge(x@snps, tmp.scores[,c("chr", "pos", "A1", "A2", "adjCADD")], by.x = c("chr", "pos", "A2", "A1"), by.y = c("chr", "pos", "A1", "A2"), all.x = T)
    scores.flip <- scores.flip.tmp$adjCADD[order(scores.flip.tmp$order)]
    x@snps$adjCADD <- scores.notflip
    x@snps$adjCADD[which(is.na(x@snps$adjCADD))] <- scores.flip[which(is.na(x@snps$adjCADD))]
    
    ###For indels that are new: attribute the adjusted value of the nearest CADD value  
    if(sum(is.na(x@snps$adjCADD))>0){
      if(!(all(c("chr", "pos", "A1", "A2", col.CADD) %in% colnames(variant.scores)))) stop(paste0("'variant.scores' should contain the columns 'chr', 'pos', 'A1', 'A2' and "), col.CADD)      
      
      ###For indels that are new: attribute the adjusted value of the nearest CADD value  
      #Annotation of PHRED by looking at A1->A2 and A2->A1
      scores.notflip <- merge(x@snps[,c("chr", "pos", "A1", "A2", "SubRegion", "adjCADD", "order")], variant.scores, by = c("chr", "pos", "A1", "A2"), all.x = T)
      #Inverse alleles to merge CADD scores when maf != @p
      scores.flip <- merge(x@snps[,c("chr", "pos", "A1", "A2", "SubRegion", "adjCADD", "order")], variant.scores, by.x = c("chr", "pos", "A2", "A1"), by.y = c("chr", "pos", "A1", "A2"), all.x = T)
      scores.flip <- scores.flip[order(scores.flip$order),]
      tmp.x <- x@snps
      tmp.x$adjCADD <- scores.notflip$adjCADD[order(scores.notflip$order)]
      tmp.x[,col.CADD] <- scores.notflip[order(scores.notflip$order),col.CADD]
      tmp.x[which(is.na(tmp.x[,col.CADD])),col.CADD] <- scores.flip[which(is.na(tmp.x[,col.CADD])),col.CADD]
      tmp.x$SubRegion <- as.character(tmp.x$SubRegion)
      #Select scores of new variants
      variant.scores.new <- tmp.x[which(is.na(x@snps$adjCADD)),]
      if(build=="b37"){
        scores.adj.indels <- fread(paste0(path.data, "/AdjustedCADD_v1.4_202204_indels.tsv.gz"), select=c(5,6,7))
      }
      if(build=="b38"){
        scores.adj.indels <- fread(paste0(path.data, "/AdjustedCADD_v1.6_202602_indels.tsv.gz"), select=c(5,6,7))
      }
      scores.byarea <- lapply(c("Coding", "Regulatory", "Intergenic"), function(z) subset(scores.adj.indels, scores.adj.indels$SubRegion==z)) 
      names(scores.byarea) <- c("Coding", "Regulatory", "Intergenic")
      #Find the nearest neighbour
      adj.scores <- mclapply(1:nrow(variant.scores.new), function(z) if(is.na(variant.scores.new[z,col.CADD])){NA}else{as.numeric(scores.byarea[[variant.scores.new[z,"SubRegion"]]][which.min(abs(as.numeric(variant.scores.new[z,col.CADD]) - scores.byarea[[variant.scores.new[z,"SubRegion"]]][,col.CADD])), "adjCADD"])}, mc.cores = cores)
      x@snps$adjCADD[which(is.na(x@snps$adjCADD))] <- unlist(adj.scores)
      x@snps$adjCADD <- as.numeric(x@snps$adjCADD)
    }
    x@snps <- x@snps[,-which(colnames(x@snps)=="order")]
    x
}
  
    
  
