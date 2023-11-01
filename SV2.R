###########################
#Obtain SV candidate reads 
#with softclip and indels
#Now also report the longest 
#indel locations for better 
#sensitivity in case blast 
#does not detect it
##########################

get_candidate_v3 <- function(reads=NULL, BAM, softclip_length) {
  library(ShortRead)
  if(is.null(reads)){
    print("Importing bam file...")
    read <- readGAlignments(BAM, use.names = T)
    print("Done importing bam.")
  } else {
    read <- reads
  }
  
  #Reads with large softclip
  dum <- read@cigar
  bad1 <- cigarRangesAlongQuerySpace(read@cigar, ops="S")
  names(bad1) <- names(read)
  foo <- unlist(bad1)
  foo1 <- foo[width(foo)>=softclip_length]
  bad1_good <- unique(names(foo1))
  
  #Also include large insertion in the middle
  bad2_r <- cigarRangesAlongReferenceSpace(read@cigar, ops="I", pos = read@start)
  bad2_q <- cigarRangesAlongQuerySpace(read@cigar, ops="I")
  
  names(bad2_r) <- names(read)
  names(bad2_q) <- names(read)
  foo_r <- unlist(bad2_r)
  foo_q <- unlist(bad2_q)
  
  ind <- which(width(foo_q)>=softclip_length)
  loc <- start(foo_r)[ind]
  read_name <- names(foo_q)[ind]
  len <- width(foo_q)[ind]
  chr <- as.character(read@seqnames)[match(read_name, names(read))]
  bad2_result <- cbind(chr, loc, len, read_name, "insertion")
  
  
  #Also include large deletion in the middle
  bad3_r <- cigarRangesAlongReferenceSpace(read@cigar, ops="D", pos = read@start)
  #bad3_q <- cigarRangesAlongQuerySpace(read@cigar, ops="D")
  
  names(bad3_r) <- names(read)
  foo_r <- unlist(bad3_r)
  
  ind <- which(width(foo_r)>=softclip_length)
  loc <- start(foo_r)[ind]
  read_name <- names(foo_r)[ind]
  len <- width(foo_r)[ind]
  chr <- as.character(read@seqnames)[match(read_name, names(read))]
  bad3_result <- cbind(chr, loc, len, read_name, "deletion")
  
  
  bad <- unique(c(bad1_good, bad2_result[,4], bad3_result[,4]))
  #browser()
  indels <- rbind(bad2_result, bad3_result)
  #indels_read <- c(names(read)[as.integer(names(bad3))], names(read)[as.integer(names(bad4))])
  write.csv(indels, file="indels.csv")
  print(paste("There are", length(bad), "reads with long softclips or large indels."))
  return(bad)
}

get_candidate_v4 <- function(BAM, softclip_length, chunk= block_size) {
  library(ShortRead)
  BAM_chunk <- BamFile(BAM, yieldSize=chunk)
  open(BAM_chunk)
  print("Importing bam file...")
  
  bad <- list()
  indels <- list()
  count <- 0
  total_read <- 0
  while(length(read <- readGAlignments(BAM_chunk, use.names = T))) {
    
    count <- count + 1
    total_read <- total_read+length(read)
    print(paste("Reading", total_read, "alignments"))
    
    #Reads with large softclip
    dum <- read@cigar
    bad1 <- cigarRangesAlongQuerySpace(read@cigar, ops="S")
    names(bad1) <- names(read)
    foo <- unlist(bad1)
    foo1 <- foo[width(foo)>=softclip_length]
    bad1_good <- unique(names(foo1))
    
    #Also include large insertion in the middle
    bad2_r <- cigarRangesAlongReferenceSpace(read@cigar, ops="I", pos = read@start)
    bad2_q <- cigarRangesAlongQuerySpace(read@cigar, ops="I")
    
    names(bad2_r) <- names(read)
    names(bad2_q) <- names(read)
    foo_r <- unlist(bad2_r)
    foo_q <- unlist(bad2_q)
    
    ind <- which(width(foo_q)>=softclip_length)
    loc <- start(foo_r)[ind]
    read_name <- names(foo_q)[ind]
    len <- width(foo_q)[ind]
    chr <- as.character(read@seqnames)[match(read_name, names(read))]
    bad2_result <- cbind(chr, loc, len, read_name, "insertion")
    
    
    #Also include large deletion in the middle
    bad3_r <- cigarRangesAlongReferenceSpace(read@cigar, ops="D", pos = read@start)
    #bad3_q <- cigarRangesAlongQuerySpace(read@cigar, ops="D")
    
    names(bad3_r) <- names(read)
    foo_r <- unlist(bad3_r)
    
    ind <- which(width(foo_r)>=softclip_length)
    loc <- start(foo_r)[ind]
    read_name <- names(foo_r)[ind]
    len <- width(foo_r)[ind]
    chr <- as.character(read@seqnames)[match(read_name, names(read))]
    bad3_result <- cbind(chr, loc, len, read_name, "deletion")
    
    
    bad[[count]] <- unique(c(bad1_good, bad2_result[,4], bad3_result[,4]))
    #browser()
    indels[[count]] <- rbind(bad2_result, bad3_result)
  } 
  close(BAM_chunk)
  write.csv(do.call("rbind", indels), file="indels.csv")
  print(paste("There are", length(unique(unlist(bad))), "reads with long softclips or large indels."))
  return(unique(unlist(bad)))
}

####################
#This is using scanBam
####################
get_candidate_v5 <- function(read, softclip_length, chunk= block_size) {
  library(ShortRead)

    #Reads with large softclip
    #read <- read[which(read$qwidth>=2500)]
    dum <- read$cigar
    bad1 <- cigarRangesAlongQuerySpace(read$cigar, ops="S")
    names(bad1) <- read$qname
    foo <- unlist(bad1)
    foo1 <- foo[width(foo)>=softclip_length]
    bad1_good <- unique(names(foo1))
    
    #Also include large insertion in the middle
    bad2_r <- cigarRangesAlongReferenceSpace(read$cigar, ops="I", pos = read$pos)
    bad2_q <- cigarRangesAlongQuerySpace(read$cigar, ops="I")
    
    names(bad2_r) <- read$qname
    names(bad2_q) <- read$qname
    foo_r <- unlist(bad2_r)
    foo_q <- unlist(bad2_q)
    
    ind <- which(width(foo_q)>=softclip_length)
    loc <- start(foo_r)[ind]
    read_name <- names(foo_q)[ind]
    len <- width(foo_q)[ind]
    chr <- as.character(read$rname)[match(read_name, read$qname)]
    bad2_result <- cbind(chr, loc, len, read_name, "insertion")
    
    
    #Also include large deletion in the middle
    bad3_r <- cigarRangesAlongReferenceSpace(read$cigar, ops="D", pos = read$pos)
    #bad3_q <- cigarRangesAlongQuerySpace(read@cigar, ops="D")
    
    names(bad3_r) <- read$qname
    foo_r <- unlist(bad3_r)
    
    ind <- which(width(foo_r)>=softclip_length)
    loc <- start(foo_r)[ind]
    read_name <- names(foo_r)[ind]
    len <- width(foo_r)[ind]
    chr <- as.character(read$rname)[match(read_name, read$qname)]
    bad3_result <- cbind(chr, loc, len, read_name, "deletion")
    
    bad <- unique(c(bad1_good, bad2_result[,4], bad3_result[,4]))
    #browser()
    indels <- rbind(bad2_result, bad3_result)

  write.csv(indels, file="indels.csv")
  print(paste("There are", length(bad), "reads with long softclips or large indels."))
  return(bad)
}

################################
#This function convert alignment
#to dot plot like image for easy 
#visualization of translocation 
#event
################################
generate_image_translocation <- function(aligned_list, output=T, cores=cores) {
  #library(rBLAST)
  if(length(aligned_list)==0) {
    print("Nothing to plot")
    return(NULL)
  }
  #for(j in 1:length(aligned_list)) {
  mclapply(1:length(aligned_list), function(j) {
    if(j %% 1000 == 0) {
      print(j)
    }
    
    dum_filtered <- as.data.frame(aligned_list[[j]])
    #dum_filtered[,c(2,3,4,6,7,9:12)] <- sapply(dum_filtered[,c(2,3,4,6,7,9:12)], as.numeric)
    dum_filtered[,3:13] <- sapply(dum_filtered[,3:13], as.numeric)
    
    if(length(unique(dum_filtered$sseqid))!=2) {
      #translocation
      print("Not translocation.")
      next
    }
    if(output) {
      png(paste0("Candidate_", j, "_plot.png"))
    }
    par(mfrow=c(1,2))
    for(j in unique(dum_filtered$sseqid)) {
      temp <- dum_filtered[dum_filtered$sseqid==j,]	
      for(i in 1:nrow(temp)) {
        min_X <- min(c(temp$qstart, temp$qend))
        max_X <- max(c(temp$qstart, temp$qend))
        min_Y <- min(c(temp$sstart, temp$send))
        max_Y <- max(c(temp$sstart, temp$send))
        length <- abs(temp$send[i] - temp$sstart[i])
        start <- temp$qstart[i]
        end <- start + length
        #print(paste(i, start, end))
        if(i == nrow(temp)) {
          plot(temp$sstart[i]:temp$send[i], start:end, cex=0.25, ylim=c(1, unique(temp$qlen)), xlim=c(min_Y, max_Y), main=names(aligned_list)[j], ylab="Query", xlab=unique(temp$sseqid))
        } else {
          plot(temp$sstart[i]:temp$send[i], start:end, cex=0.25, ylim=c(1, unique(temp$qlen)), xlim=c(min_Y, max_Y), main="", xlab="", ylab="")
        }
        par(new=T)
      }
      par(new=F)
    }
    if(output) {
      dev.off()
    }
  }, mc.cores=)
}

################################
#This function convert alignment
#to dot plot like image for easy 
#visualization of SV event
################################
generate_image <- function(aligned_list, output=T, cores=cores) {
  #library(rBLAST)
  if(length(aligned_list)==0) {
    print("Nothing to plot")
    return(NULL)
  }
  #browser()
  #for(j in 1:length(aligned_list)) {
  mclapply(1:length(aligned_list), function(j) {
    if(j %% 1000 == 0) {
      print(j)
    }
    
    dum_filtered <- as.data.frame(aligned_list[[j]])
    sseqid <- unique(dum_filtered$sseqid)
    dum_filtered[,c(3:13)] <- sapply(dum_filtered[,c(3:13)], as.numeric)
    #dum_filtered[,c(2,3,4,6,7,9:12)] <- sapply(dum_filtered[,c(2,3,4,6,7,9:12)], as.numeric)
    #if(length(unique(dum_filtered$sseqid))>1) {
    # print("More than one chromosome.")
    #  next
    #}
    if(output) {
      png(paste0("Candidate_", j, "_plot.png"))
    }
    
    for(i in 1:nrow(dum_filtered)) {
      
      min_X <- min(c(dum_filtered$qstart, dum_filtered$qend))
      max_X <- max(c(dum_filtered$qstart, dum_filtered$qend))
      min_Y <- min(c(dum_filtered$sstart, dum_filtered$send))
      max_Y <- max(c(dum_filtered$sstart, dum_filtered$send))
      length <- abs(dum_filtered$send[i] - dum_filtered$sstart[i])
      start <- dum_filtered$qstart[i]
      end <- start + length
      #print(paste(i, start, end))
      if(i == nrow(dum_filtered)) {
        plot(dum_filtered$sstart[i]:dum_filtered$send[i], start:end, cex=0.25, ylim=c(1, unique(dum_filtered$qlen)), xlim=c(min_Y, max_Y), main="", ylab=names(aligned_list)[j], xlab=sseqid)
      } else {
        plot(dum_filtered$sstart[i]:dum_filtered$send[i], start:end, cex=0.25, ylim=c(1, unique(dum_filtered$qlen)), xlim=c(min_Y, max_Y), main="", xlab="", ylab="")
      }
      par(new=T)
    }
    if(output) {
      dev.off()
    }
  }, mc.cores=)
}

#############################
#This function converts blast 
#alignment to projection on 
#reference and query
#############################
convert_alignment_new <- function(alignment) {
  
  if(length(unique(alignment$sseqid))>1) {
    print("More than one chromosome.")
    print(alignment)
    return(NULL)
  }
  alignment[,c(3:13)] <- sapply(alignment[,c(3:13)], as.numeric)
  
  #reference
  foo <- alignment[,c("sseqid", "sstart", "send")]
  foo[,c("sstart", "send")] <- t(apply(foo[,c("sstart", "send")], 1, function(X) sort(as.numeric(X), decreasing = F)))
  colnames(foo) <- c("seqname", "start", "end")
  foo1 <- GRanges(foo)
  a1 <- as.integer(coverage(foo1)[[1]])[-c(1:min(foo[,c("start", "end")]))]
  
  #query
  foo <- alignment[,c("qseqid", "qstart", "qend")]
  foo[,c("qstart", "qend")] <- t(apply(foo[,c("qstart", "qend")], 1, function(X) sort(as.numeric(X), decreasing = F)))
  colnames(foo) <- c("seqname", "start", "end")
  foo1 <- GRanges(foo)
  b1 <- as.integer(coverage(foo1)[[1]])[-c(1:min(foo[,c("start", "end")]))]
  
  return(list(a1, b1))
}

######################
#It is the function 
#to summarize detected 
#SV events
######################
Collapse_SV_v2 <- function(result_table) {
  #result_table <- result[result[,5]!="translocation",]
  result_table <- as.data.frame(result_table)
  result_table <- result_table[result_table[,5]!="translocation",] #need a separate solution for TRA
  result_table[,1] <- as.character(result_table[,1])
  result_table[,2:3] <- t(apply(result_table[,2:3], 1, function(X) as.numeric(X)))
  result_table[,3] <- result_table[,2]+result_table[,3]
  result_table[,4] <- as.character(result_table[,4])
  result_table[,5] <- as.character(result_table[,5])
  
  
  colnames(result_table) <- c("seqname", "start", "end", "read_id", "sv_type")
  temp <- split(result_table, f=result_table$sv_type)
  summary_report <- list()
  for(i in 1:length(temp)) {
    result_temp <- temp[[i]]
    result_temp <- result_temp[order(result_temp$seqname, result_temp$start),]
    result_range <- GRanges(result_temp)
    combined <- reduce(result_range, min.gapwidth=25)
    foo <- findOverlaps(result_range, combined)
    combined$support_reads <- table(subjectHits(foo))
    combined$sv_type <- unique(result_temp$sv_type)
    print(table(combined$sv_type))
    summary_report[[i]] <- data.frame(Chr=seqnames(combined), Location=start(combined), Length=width(combined), sv_type=combined$sv_type, supporting_reads=as.numeric(combined$support_reads))
  }
  summary_report <- do.call("rbind", summary_report)
  return(summary_report)
}

##################################
#Using hs-blastn to improve on speed
#Now record the breakpoint at real time
##################################
find_sv_hs_blastn_parallel_v2 <- function(FASTQ, candidates= NULL, min_sv_size = 30, max_dist=10000, blast_ref="Z:/genome/SacCer3/SacCer3_blast", masker_db=NULL, block_size=1000, platform="PacBio CCS", graph=T, node=8) {
  #install.packages('rBLAST', repos = 'https://mhahsler.r-universe.dev')
  library(ShortRead)
  library(parallel)
  #library(tictoc)
  #browser()
  #Sys.setenv(PATH=paste("/opt/hs-blastn/0.0.5/bin/", Sys.getenv("PATH"), sep=":"))
  #tic()
  indels <- read.csv("indels.csv")
  indels <- indels[,-1]
  #colnames(indels) <- ""
  #candidates=indels[,4]
  
  fs <- FastqStreamer(FASTQ, n=block_size)
  aligned_result <- list()
  count <- 0
  streamer_count <-  0
  temp_result <- list()
  print("Start detection process")
  while(length(fq <- yield(fs))) {
    #browser()
    streamer_count <- streamer_count+1
    print(paste("Imported", streamer_count*block_size, "reads."))
    temp_id <- sapply(strsplit((as.character(fq@id)), split=" "), function(X) X[[1]][1])
    ind <- which(is.element(temp_id, candidates)==T)
    if(length(ind)==0) {
      next
    }
    count <- count+length(ind)
    print(paste("Processing", length(ind), "candidates..."))
    a <- DNAStringSet(fq@sread[ind])
    names(a) <- temp_id[ind]
   
    width <- width(a)
    names(width) <- names(a)
    writeXStringSet(a, file="temp.fa")
    #aligned <- predict(ref, a, BLAST_args = paste("-task megablast -window_masker_db /net/nfs-irwrsrchnas01/labs/xwu/genome/hg38/hg38.fa.counts.obinary -num_threads", node), custom_format = "qseqid qlen qstart qend sseqid sstart send sstrand length pident bitscore")
    if(is.null(masker_db)) {
      #this will be too slow
      system(paste("hs-blastn align -query temp.fa -db", blast_ref, "-outfmt 6 -evalue 1e-10 -max_target_seqs 20 -num_threads", node, "> temp.out"))
    } else {
      #Fastest megablast
      system(paste("hs-blastn align -query temp.fa -db", blast_ref, "-window_masker_db", masker_db, "-max_target_seqs 20 -outfmt 6 -evalue 1e-10 -num_threads", node, "> temp.out"), ignore.stderr = T)
    }
    if(file.size("temp.out")) {
      aligned <- read.delim("temp.out", header=F)
    } else {
      next
    }
    aligned$qlen <- width[match(aligned[,1], names(width))]
    colnames(aligned) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen")
    if(platform=="PacBio CCS") {
      #use higher pindet for ccs data, remove chrM by default
      aligned <- aligned[(aligned$qend-aligned$qstart)>=min_sv_size&aligned$pident>=94&aligned$sseqid!="chrM",]
    } else {
      #allow more errors with non-CCS data, e.g. Nanopore or PacBio CLS
      aligned <- aligned[(aligned$qend-aligned$qstart)>=min_sv_size&aligned$pident>=87&aligned$sseqid!="chrM",]
    }
    aligned[,c("qstart", "qend")] <- t(apply(aligned[,c("qstart", "qend")], 1, function(X) sort(as.numeric(X), decreasing = F)))
    aligned <- aligned[order(aligned$qseqid, aligned$qstart),]
    aligned_list <- split(aligned, f=aligned$qseqid)
    print("Finished alignment with Blast.")

    temp_result[[streamer_count]] <- mclapply(1:length(aligned_list), function(i) {
      dum <- aligned_list[[i]]
      dum$sstrand <- apply(dum[,c("sstart", "send")], 1, function(X) ifelse(X[1]< X[2], "pos", "neg"))
      dum$cov <- (dum$qend - dum$qstart +1)/dum$qlen
      
      #Need to refine the subregion to make sure there is no redundant alignments
      if(nrow(dum)>2) {
        dum_temp <- dum
        dum_temp[,c("qstart", "qend")] <- t(apply(dum_temp[,c("qstart", "qend")], 1, function(X) sort(as.numeric(X), decreasing = F)))
        x <- GRanges(seqnames=unique(dum_temp$qseqid), ranges = IRanges(start=dum_temp$qstart, end=dum_temp$qend))
        y <- reduce(x-round(min_sv_size/2))
        foo <- as.data.frame(findOverlaps(x,y, minoverlap = min_sv_size/2))
        segmentation <- list()
        for(j in unique(foo$subjectHits)) {
          segmentation_temp <- dum_temp[foo$queryHits[which(foo$subjectHits==j)],]
          segmentation[[j]] <- segmentation_temp[which(segmentation_temp$bitscore==max(segmentation_temp$bitscore)),]
        }
        dum <- do.call("rbind", segmentation)
        #updated segmentation step to remove redundant, not as good
        #y <- findOverlaps(x,x+min_sv_size/2, type="within") 
        #a1 <- as.matrix(y)
        #b1 <- a1[which(a1[,1]!=a1[,2]),]
        #dum <- x[-unique(b1[,1]),]
      }
      best <- which(dum$bitscore==max(dum$bitscore))
      if((max(dum$cov[best])>=0.99 & max(dum$pident[best])>=97) | length(unique(dum$sseqid))>2) {
        #perfect align or too complicated align
        return(NULL)
      }
      
      if(length(unique(dum$sseqid))==2) {
        #browser()
        #1. Check if one chrom already covers entire read, go to following step for segmentation
        foo <- split(dum, f=dum$sseqid)
        each_list <- list()
        for(z in 1:length(foo)) {
          each <- array(0, dim=c(unique(foo[[z]]$qlen), nrow(foo[[z]])))
          for(j in 1:nrow(foo[[z]])) {
            each[foo[[z]][j, "qstart"]:foo[[z]][j, "qend"],j] <- 1
          }
          each_list[[z]] <- each
        }
        cov_list <- sapply(each_list, function(each) sum(apply(each, 1, sum)>0)/nrow(each))
        if(max(cov_list)>=0.99) {
          #one chrom covers entire read, so just use the chrom for next step and ignore the other chrom
          dum <- foo[[which(cov_list>= 0.99)]]
        } else {
          #Check if fragments on both chrom are contiguous and covers most of the read
          each_combined <- do.call("cbind", each_list)
          cov_combined <- sum(apply(each_combined, 1, sum)>0)/nrow(each_combined)
          if(cov_combined>= 0.99) {
            #translocation
            dum$SV <- "translocation"
            align_good <- dum
            dum <- dum[order(dum$sseqid),]
            count<- count+1
            if(dum$sstrand[1]=="pos") {
              start <- dum$send[1]
            } else {
              start <- dum$sstart[1]
            }  
            
            if(dum$sstrand[2]=="pos") {
              end <- dum$sstart[2]
            } else {
              end <- dum$sstart[2]
            }
            seqname <- paste(unique(dum$sseqid), collapse="-")
            read <- unique(dum$qseqid)
            report_good <- c(seqname, start, paste(start, end, sep=":"), read, "translocation")
            return(list(align_good, report_good))
          } else {
            #return(NULL)
            ind <- which(indels[,4]==unique(dum1$qseqid))
            if(length(ind)>0) {
              dum2 <- list()
              for(i in 1:length(ind)) {
                dum3 <- dum1
                dum3$SV <- indels[ind[i],5]
                dum2[[i]] <- dum3
              }
              return(list(do.call("rbind", dum2), indels[ind,]))
            }
          }
        }
      }
      
      #add in segmentation step to resolve multiple reports involving distant regions
      dum_temp <- dum
      dum_temp[,c("sstart", "send")] <- t(apply(dum_temp[,c("sstart", "send")], 1, function(X) sort(as.numeric(X), decreasing = F)))
      x <- GRanges(seqnames=unique(dum_temp$sseqid), ranges = IRanges(start=dum_temp$sstart, end=dum_temp$send))
      y <- reduce(x+max_dist)
      foo <- as.data.frame(findOverlaps(x, y))
      segmentation <- array(dim=length(unique(foo$subjectHits)))
      for(j in unique(foo$subjectHits)) {
        segmentation_temp <- x[foo$queryHits[which(foo$subjectHits==j)]]
        segmentation[j] <- sum(width(segmentation_temp))
      }
      
      idx <- which(segmentation==max(segmentation))
      dum1 <- dum[foo$queryHits[which(foo$subjectHits==idx)],]
      
      #Now checking for SV 
      if(nrow(dum1)==1) {
        #when one size of query not aligned, it is an insertion, although not completed detected.
        #if((dum1$qlen-dum1$cov*dum1$qlen)>= min_sv_size*10){
        #  dum1$SV <- "BP_other"
        #  align_good <- dum1
        #  if(dum1$qstart<(dum1$qlen-dum1$qend)) {
        #    start <- as.numeric(dum1$send)
        #    len <- dum1$qlen-dum1$qend
        #  } else {
        #    start <- as.numeric(dum1$sstart)
        #    len <- dum1$qstart
        #  }
        #  seqname <- unique(dum1$sseqid)
        #  read <- unique(dum1$qseqid)
        #  report_good <- c(seqname, start, len, read, "BP_other")
        #  return(list(align_good, report_good))
        #}

        #ind <- match(dum1$qseqid, indels[,4])
        #if(!is.na(ind)) {
        #  dum1$SV <- indels[ind,5]
        #  return(list(dum1, as.character(indels[ind,])))
        #}
        
        #note that some reads have multiple entries reported from parsing pbmm2
        ind <- which(indels[,4]==dum1$qseqid)
        if(length(ind)>0) {
          dum2 <- list()
          for(i in 1:length(ind)) {
            dum3 <- dum1
            dum3$SV <- indels[ind[i],5]
            dum2[[i]] <- dum3
          }
          return(list(do.call("rbind", dum2), indels[ind,]))
        }
        #} else if  (nrow(dum1)>2 &nrow(dum) <=4 & length(unique(dum1$sstrand))==2) {
      } else if  (nrow(dum1)==3 & (sum(dum1$sstrand==c("pos", "neg", "pos"))==3 | sum(dum1$sstrand==c("neg", "pos", "neg"))==3)) {
        temp <- convert_alignment_new(dum1)
        a1 <- temp[[1]] #reference
        b1 <- temp[[2]] #query
        
        if (sum(b1==0)<=min_sv_size & sum(a1==0)<=min_sv_size & sum(b1>1)<=min_sv_size & sum(a1>1)>=min_sv_size) {
          #print("IVT")
          dum1$SV <- "IVT"
          align_good <- dum1
          
          if(dum1$sstrand[2]=="pos") {
            start <- dum1$sstart[2]
            end <- dum1$send[2]
          } else {
            end <- dum1$sstart[2]
            start <- dum1$send[2]
          }
          seqname <- unique(dum1$sseqid)
          read <- unique(dum1$qseqid)
          report_good <- c(seqname, start, abs(start-end), read, "IVT")
          
          return(list(align_good, report_good))
        } else {
          dum1$SV <- "inversion"
          align_good <- dum1
          
          if(dum1$sstrand[2]=="pos") {
            start <- dum1$sstart[2]
            end <- dum1$send[2]
          } else {
            end <- dum1$sstart[2]
            start <- dum1$send[2]
          }
          seqname <- unique(dum1$sseqid)
          read <- unique(dum1$qseqid)
          report_good <- c(seqname, start, abs(start-end), read, "inversion")
          
          return(list(align_good, report_good))
        }
      } else if  (nrow(dum1)==2 & length(unique(dum1$sstrand))==2) {
        temp <- convert_alignment_new(dum1)
        a1 <- temp[[1]] #reference
        b1 <- temp[[2]] #query
        
        if (sum(b1==0)<=min_sv_size & sum(a1==0)<=min_sv_size & sum(b1>1)<=min_sv_size) {
          dum1$SV <- "inversion_other"
          align_good <- dum1
          
          if(dum$sstrand[1]=="pos") {
            start <- dum1$send[1]
          } else {
            start <- dum1$sstart[1]
          }  
          
          if(dum1$sstrand[2]=="pos") {
            end <- dum1$sstart[2]
          } else {
            end <- dum1$sstart[2]
          }
          
          seqname <- unique(dum1$sseqid)
          read <- unique(dum1$qseqid)
          report_good <- c(seqname, start, abs(start-end), read, "inversion_other")
          #return(NULL)
          return(list(align_good, report_good))
        }
      } else if(nrow(dum1)==2 & length(unique(dum1$sstrand))==1) {
        temp <- convert_alignment_new(dum1)
        a1 <- temp[[1]] #reference
        b1 <- temp[[2]] #query
        
        #if(sum(a1!=1)<min_sv_size & sum(b1!=1)<min_sv_size) {
        #  return(NULL)
        #}
        
        if (sum(a1==0)<=min_sv_size & sum(b1==0)<=min_sv_size & sum(a1>1)>=min_sv_size &sum(b1>1)<=min_sv_size) {
        #  if (sum(a1==0)<=min_sv_size & sum(b1==0)<=min_sv_size & sum(a1>1)>=min_sv_size) {
            #likely miss very large duplications, due to read length
          dum1$SV <- "duplication"
          align_good <- dum1
          foo1 <- range(which(a1>1))
          start <- min(dum1[,c("sstart", "send")])+foo1[1]
          end <- min(dum1[,c("sstart", "send")])+foo1[2]
          seqname <- unique(dum1$sseqid)
          read <- unique(dum1$qseqid)
          report_good <- c(seqname, start, abs(start-end), read, "duplication")
          
          return(list(align_good, report_good))
        } else if (sum(a1==0)<=min_sv_size & sum(a1>1) <= min_sv_size & sum(b1==0)>=min_sv_size&sum(b1>1)<=min_sv_size) {
        #} else if (sum(a1==0)<sum(b1==0) & sum(b1==0)>=min_sv_size) {
          dum1$SV <- "insertion"
          align_good <- dum1
          
          start <- round(mean(as.numeric(dum1$send[1]), as.numeric(dum1$sstart[2])))
          seqname <- unique(dum1$sseqid)
          read <- unique(dum1$qseqid)
          report_good <- c(seqname, start, abs(dum1$qstart[2]-dum1$qend[1]), read, "insertion")
          return(list(align_good, report_good))
        } else if(sum(b1==0)<=min_sv_size & sum(b1>1) <=min_sv_size & sum(a1==0)>=min_sv_size & sum(a1>1)<=min_sv_size) {
        #} else if(sum(b1==0)<sum(a1==0) & sum(a1==0)>=min_sv_size) {
          dum1$SV <- "deletion"
          align_good <- dum1
          if(unique(dum1$sstrand)=="pos") {
            start <- dum1$send[1]
            end <- dum1$sstart[2]
          } else {
            end <- dum1$send[1]
            start <- dum1$sstart[2]
          }
          seqname <- unique(dum1$sseqid)
          read <- unique(dum1$qseqid)
          report_good <- c(seqname, start, abs(start-end), read, "deletion")
          return(list(align_good, report_good))
        } else if (sum(a1==0)>=min_sv_size & sum(a1>1) <= min_sv_size & sum(b1==0)>=min_sv_size&sum(b1>1)<=min_sv_size) {
          dum1$SV <- "CSUB"
          align_good <- dum1
          if(unique(dum1$sstrand)=="pos") {
            start <- dum1$send[1]
            end <- dum1$sstart[2]
          } else {
            end <- dum1$send[1]
            start <- dum1$sstart[2]
          }
          seqname <- unique(dum1$sseqid)
          read <- unique(dum1$qseqid)
          report_good <- c(seqname, start, abs(start-end), read, "CSUB")
          return(list(align_good, report_good))
        } else {
          #dum1$SV <- "BP_other"
          #return(NULL)
          ind <- which(indels[,4]==unique(dum1$qseqid))
          if(length(ind)>0) {
            dum2 <- list()
            for(i in 1:length(ind)) {
              dum3 <- dum1
              dum3$SV <- indels[ind[i],5]
              dum2[[i]] <- dum3
            }
            #print(indels[good,])
            return(list(do.call("rbind", dum2), indels[ind,]))
          }
        }
        
      } else if(nrow(dum1) >=3 & length(unique(dum1$sstrand))==1) {
        temp <- convert_alignment_new(dum1)
        a1 <- temp[[1]] #reference
        b1 <- temp[[2]] #query
        
        if (sum(b1==0)<=min_sv_size & sum(a1>1)>=min_sv_size & sum(a1==0)<=min_sv_size) {
          dum1$SV <- "duplication"
          align_good <- dum1
          foo1 <- range(which(a1>1))
          start <- min(dum1[,c("sstart", "send")])+foo1[1]
          end <- min(dum1[,c("sstart", "send")])+foo1[2]
          seqname <- unique(dum1$sseqid)
          read <- unique(dum1$qseqid)
          report_good <- c(seqname, start, abs(start-end), read, "duplication")
          
          return(list(align_good, report_good))
        } else if (sum(b1==0)<=min_sv_size & sum(a1==0)>=min_sv_size) {
          dum1$SV <- "deletion"
          align_good <- dum1
          if(unique(dum1$sstrand)=="pos") {
            start <- dum1$send[1]
            end <- dum1$sstart[2]
          } else {
            end <- dum1$send[1]
            start <- dum1$sstart[2]
          }
          seqname <- unique(dum1$sseqid)
          read <- unique(dum1$qseqid)
          report_good <- c(seqname, start, abs(start-end), read, "deletion")
          return(list(align_good, report_good))
        } 
        #else if (sum(a1==0)<=min_sv_size & sum(b1>=2)>=min_sv_size) {
          #repeats in reference
        #  return(NULL)
        #} 
      else {
          ind <- which(indels[,4]==unique(dum1$qseqid))
          if(length(ind)>0) {
            dum2 <- list()
            for(i in 1:length(ind)) {
              dum3 <- dum1
              dum3$SV <- indels[ind[i],5]
              dum2[[i]] <- dum3
            }
            #print(indels[good,])
            return(list(do.call("rbind", dum2), indels[ind,]))
          }
        }
      } else {
        ind <- which(indels[,4]==unique(dum1$qseqid))
        if(length(ind)>0) {
          dum2 <- list()
          for(i in 1:length(ind)) {
            dum3 <- dum1
            dum3$SV <- indels[ind[i],5]
            dum2[[i]] <- dum3
          }
          #print(indels[good,])
          return(list(do.call("rbind", dum2), indels[ind,]))
        }
      }
      #result <- list(insertion=insertion, deletion=deletion,inversion=inversion, duplication=duplication, tandem=tandem_repeats, translocation=translocation, other=other)
      #return(dum1)
    }, mc.cores=node)
    
    good_num <- sum(sapply(temp_result,function(X) sum(sapply(X, function(Y) length(Y[[1]]))!=0)))
    print(paste("Found", good_num, "SV reads so far."))
    gc()
  }
  #toc()
  close(fs)
  saveRDS(temp_result, file="temp_result.rds")
  #browser()
  good <- lapply(temp_result, function(X) X[sapply(X, function(Y) !is.null(Y)&sum(is.na(Y))==0&class(Y)!="try-error")])
  result <- do.call("rbind", lapply(good, function(X) do.call("rbind", sapply(X, function(Y) Y[1]))))
  a <- split(result, f=result$SV)
  result_list <- lapply(a, function(X) split(X, f=X$qseqid))
  print(sapply(result_list, length))
  #saveRDS(result_list, file="sv_all.rds")
  #report <- do.call("rbind", lapply(good, function(X) do.call("rbind", sapply(X, function(Y) Y[2]))))
  report <- lapply(good, function(X) do.call("rbind", sapply(X, function(Y) Y[2])))
  report_filtered <- report[!sapply(report, is.null)]
  for(z in 1:length(report_filtered)){
    colnames(report_filtered[[z]]) <- colnames(indels)
  }
  report <- do.call("rbind", report_filtered)
  #write.csv(report, file="sv_report.csv")
  return(list(result_list, report_filtered))
}

find_sv_hs_blastn_parallel_v3 <- function(FASTQ, candidates= NULL, min_sv_size = 30, max_dist=10000, blast_ref="Z:/genome/SacCer3/SacCer3_blast", masker_db=NULL, block_size=1000, platform="PacBio CCS", graph=T, node=8) {
  #install.packages('rBLAST', repos = 'https://mhahsler.r-universe.dev')
  library(ShortRead)
  library(parallel)

  #Sys.setenv(PATH=paste("/opt/hs-blastn/0.0.5/bin/", Sys.getenv("PATH"), sep=":"))
  #tic()
  indels <- read.csv("indels.csv")
  indels <- indels[,-1]
  #candidates=indels[,4]
  
  aligned_result <- list()
  count <- 0

  #print("Start detection process")
  fq <- FASTQ[width(FASTQ)>=2500]
    
    ind <- which(is.element(names(fq), candidates)==T)
    if(length(ind)==0) {
      return(NULL)
    }
    count <- count+length(ind)
    
    a <- fq[ind]
    #b <- a[width(a)>=2500]
    b <- a[match(unique(names(a)), names(a))]
    print(paste("Processing", length(b), "candidates..."))
    #browser()
    width <- width(b)
    names(width) <- names(b)
    writeXStringSet(b, file="temp.fa")
    #aligned <- predict(ref, a, BLAST_args = paste("-task megablast -window_masker_db /net/nfs-irwrsrchnas01/labs/xwu/genome/hg38/hg38.fa.counts.obinary -num_threads", node), custom_format = "qseqid qlen qstart qend sseqid sstart send sstrand length pident bitscore")
    if(is.null(masker_db)) {
      system(paste("hs-blastn align -query temp.fa -db", blast_ref, "-outfmt 6 -evalue 1e-10 -max_target_seqs 20 -num_threads", min(length(b), node), "> temp.out"))
    } else {
      system(paste("hs-blastn align -query temp.fa -db", blast_ref, "-window_masker_db", masker_db, "-max_target_seqs 20 -outfmt 6 -evalue 1e-10 -num_threads", min(length(b), node), "> temp.out"), ignore.stderr = T)
    }
    if(file.size("temp.out")) {
      aligned <- read.delim("temp.out", header=F)
    } else {
      return(NULL)
    }
    aligned$qlen <- width[match(aligned[,1], names(width))]
    colnames(aligned) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen")
    
    if(platform=="PacBio CCS") {
      aligned <- aligned[aligned$bitscore>=2*min_sv_size&aligned$pident>=94&aligned$sseqid!="chrM",]
    } else {
      aligned <- aligned[aligned$bitscore>=2*min_sv_size&aligned$pident>=87&aligned$sseqid!="chrM",]
    }
    aligned[,c("qstart", "qend")] <- t(apply(aligned[,c("qstart", "qend")], 1, function(X) sort(as.numeric(X), decreasing = F)))
    aligned <- aligned[order(aligned$qseqid, aligned$qstart),]
    aligned_list <- split(aligned, f=aligned$qseqid)
    print("Finished alignment with Blast.")
    #browser()
    temp_result <- mclapply(1:length(aligned_list), function(i) {
      #if(i %% 1000 ==0) {
      #  print(i)
      #}
      
      dum <- aligned_list[[i]]
      dum$sstrand <- apply(dum[,c("sstart", "send")], 1, function(X) ifelse(X[1]< X[2], "pos", "neg"))
      dum$cov <- (dum$qend - dum$qstart +1)/dum$qlen
      
      #Need to refine the subregion to make sure there is no redundant alignments
      if(nrow(dum)>2) {
        dum_temp <- dum
        dum_temp[,c("qstart", "qend")] <- t(apply(dum_temp[,c("qstart", "qend")], 1, function(X) sort(as.numeric(X), decreasing = F)))
        x <- GRanges(seqnames=unique(dum_temp$qseqid), ranges = IRanges(start=dum_temp$qstart, end=dum_temp$qend))
        y <- reduce(x-round(min_sv_size/2))
        #y <- reduce(x, min.gapwidth=round(min_sv_size/2))
        foo <- as.data.frame(findOverlaps(x,y, minoverlap = min_sv_size/2))
        #foo <- as.data.frame(findOverlaps(x,y))
        segmentation <- list()
        for(j in unique(foo$subjectHits)) {
          segmentation_temp <- dum_temp[foo$queryHits[which(foo$subjectHits==j)],]
          #segmentation[[j]] <- segmentation_temp
          segmentation[[j]] <- segmentation_temp[which(segmentation_temp$bitscore==max(segmentation_temp$bitscore)),]
        }
        dum <- do.call("rbind", segmentation)
      }
      
      best <- which(dum$bitscore==max(dum$bitscore))
      
      if((max(dum$cov[best])>=0.99 & max(dum$pident[best])>=97) | length(unique(dum$sseqid))>2) {
        #perfect align or too complicated align
        return(NULL)
      }
      
      if(length(unique(dum$sseqid))==2) {
        #browser()
        #1. Check if one chrom already covers entire read, go to following step for segmentation
        foo <- split(dum, f=dum$sseqid)
        each_list <- list()
        for(z in 1:length(foo)) {
          each <- array(0, dim=c(unique(foo[[z]]$qlen), nrow(foo[[z]])))
          for(j in 1:nrow(foo[[z]])) {
            each[foo[[z]][j, "qstart"]:foo[[z]][j, "qend"],j] <- 1
          }
          each_list[[z]] <- each
        }
        cov_list <- sapply(each_list, function(each) sum(apply(each, 1, sum)>0)/nrow(each))
        if(max(cov_list)>=0.99) {
          #one chrom covers entire read, so just use the chrom for next step and ignore the other chrom
          dum <- foo[[which(cov_list>= 0.99)]]
        } else {
          #Check if fragments on both chrom are contiguous and covers most of the read
          each_combined <- do.call("cbind", each_list)
          cov_combined <- sum(apply(each_combined, 1, sum)>0)/nrow(each_combined)
          if(cov_combined>= 0.99) {
            #translocation
            dum$SV <- "translocation"
            align_good <- dum
            dum <- dum[order(dum$sseqid),]
            count<- count+1
            if(dum$sstrand[1]=="pos") {
              start <- dum$send[1]
            } else {
              start <- dum$sstart[1]
            }  
            
            if(dum$sstrand[2]=="pos") {
              end <- dum$sstart[2]
            } else {
              end <- dum$sstart[2]
            }
            seqname <- paste(unique(dum$sseqid), collapse="-")
            read <- unique(dum$qseqid)
            report_good <- c(seqname, start, paste(start, end, sep=":"), read, "translocation")
            return(list(align_good, report_good))
          } else {
            return(NULL)
          }
        }
      }
      
      #add in segmentation step to resolve multiple reports involving distant regions
      dum_temp <- dum
      dum_temp[,c("sstart", "send")] <- t(apply(dum_temp[,c("sstart", "send")], 1, function(X) sort(as.numeric(X), decreasing = F)))
      x <- GRanges(seqnames=unique(dum_temp$sseqid), ranges = IRanges(start=dum_temp$sstart, end=dum_temp$send))
      y <- reduce(x+max_dist)
      foo <- as.data.frame(findOverlaps(x, y))
      segmentation <- array(dim=length(unique(foo$subjectHits)))
      for(j in unique(foo$subjectHits)) {
        segmentation_temp <- x[foo$queryHits[which(foo$subjectHits==j)]]
        segmentation[j] <- sum(width(segmentation_temp))
      }
      
      idx <- which(segmentation==max(segmentation))
      dum1 <- dum[foo$queryHits[which(foo$subjectHits==idx)],]
      
      #Now checking for SV 
      if(nrow(dum1)==1) {
        ind <- match(dum1$qseqid, indels[,4])
        if(!is.na(ind)) {
          #print("extra")
          dum1$SV <- indels[ind,5]
          return(list(dum1, as.character(indels[ind,])))
        }
        #normal
        #return(NULL)
        #} else if  (nrow(dum1)>2 &nrow(dum) <=4 & length(unique(dum1$sstrand))==2) {
      } else if  (nrow(dum1)==3 & (sum(dum1$sstrand==c("pos", "neg", "pos"))==3 | sum(dum1$sstrand==c("neg", "pos", "neg"))==3)) {
        temp <- convert_alignment_new(dum1)
        a1 <- temp[[1]] #reference
        b1 <- temp[[2]] #query
        
        if (sum(b1==0)<=min_sv_size & sum(a1==0)<=min_sv_size & sum(b1>1)<=min_sv_size) {
          dum1$SV <- "inversion"
          align_good <- dum1
          
          if(dum1$sstrand[2]=="pos") {
            start <- dum1$sstart[2]
            end <- dum1$send[2]
          } else {
            end <- dum1$sstart[2]
            start <- dum1$send[2]
          }
          seqname <- unique(dum1$sseqid)
          read <- unique(dum1$qseqid)
          report_good <- c(seqname, start, abs(start-end), read, "inversion")
          
          return(list(align_good, report_good))
        }
      } else if  (nrow(dum1)==2 & length(unique(dum1$sstrand))==2) {
        temp <- convert_alignment_new(dum1)
        a1 <- temp[[1]] #reference
        b1 <- temp[[2]] #query
        
        if (sum(b1==0)<=min_sv_size & sum(a1==0)<=min_sv_size & sum(b1>1)<=min_sv_size) {
          dum1$SV <- "inversion_other"
          align_good <- dum1
          
          if(dum$sstrand[1]=="pos") {
            start <- dum1$send[1]
          } else {
            start <- dum1$sstart[1]
          }  
          
          if(dum1$sstrand[2]=="pos") {
            end <- dum1$sstart[2]
          } else {
            end <- dum1$sstart[2]
          }
          
          seqname <- unique(dum1$sseqid)
          read <- unique(dum1$qseqid)
          report_good <- c(seqname, start, abs(start-end), read, "inversion_other")
          #return(NULL)
          return(list(align_good, report_good))
        }
      } else if(nrow(dum1)==2 & length(unique(dum1$sstrand))==1) {
        temp <- convert_alignment_new(dum1)
        a1 <- temp[[1]] #reference
        b1 <- temp[[2]] #query
        
        if(sum(a1!=1)<min_sv_size & sum(b1!=1)<min_sv_size) {
          return(NULL)
        }
        
        if (sum(a1==0)<=min_sv_size & sum(b1==0)<=min_sv_size & sum(a1>1)>=min_sv_size &sum(b1>1)<=min_sv_size) {
          dum1$SV <- "duplication"
          align_good <- dum1
          foo1 <- range(which(a1>1))
          start <- min(dum1[,c("sstart", "send")])+foo1[1]
          end <- min(dum1[,c("sstart", "send")])+foo1[2]
          seqname <- unique(dum1$sseqid)
          read <- unique(dum1$qseqid)
          report_good <- c(seqname, start, abs(start-end), read, "duplication")
          
          return(list(align_good, report_good))
        } else if (sum(a1==0)>=min_sv_size & sum(a1>1) <= min_sv_size & sum(b1==0)>=min_sv_size&sum(b1>1)<=min_sv_size) {
          dum1$SV <- "CSUB"
          align_good <- dum1
          if(unique(dum1$sstrand)=="pos") {
            start <- dum1$send[1]
            end <- dum1$sstart[2]
          } else {
            end <- dum1$send[1]
            start <- dum1$sstart[2]
          }
          seqname <- unique(dum1$sseqid)
          read <- unique(dum1$qseqid)
          report_good <- c(seqname, start, abs(start-end), read, "CSUB")
          return(list(align_good, report_good))
        } else if (sum(a1==0)<=min_sv_size/2 & sum(a1>1) <= min_sv_size & sum(b1==0)>=min_sv_size&sum(b1>1)<=min_sv_size) {
          dum1$SV <- "insertion"
          align_good <- dum1
          
          start <- round(mean(as.numeric(dum1$send[1]), as.numeric(dum1$sstart[2])))
          seqname <- unique(dum1$sseqid)
          read <- unique(dum1$qseqid)
          report_good <- c(seqname, start, abs(dum1$qstart[2]-dum1$qend[1]), read, "insertion")
          return(list(align_good, report_good))
        } else if(sum(b1==0)<=min_sv_size & sum(b1>1) <=min_sv_size & sum(a1==0)>=min_sv_size & sum(a1>1)<=min_sv_size) {
          dum1$SV <- "deletion"
          align_good <- dum1
          if(unique(dum1$sstrand)=="pos") {
            start <- dum1$send[1]
            end <- dum1$sstart[2]
          } else {
            end <- dum1$send[1]
            start <- dum1$sstart[2]
          }
          seqname <- unique(dum1$sseqid)
          read <- unique(dum1$qseqid)
          report_good <- c(seqname, start, abs(start-end), read, "deletion")
          return(list(align_good, report_good))
        } else {
          #dum1$SV <- "other"
          return(NULL)
        }
        
      } else if(nrow(dum1) >=3 & length(unique(dum1$sstrand))==1) {
        temp <- convert_alignment_new(dum1)
        a1 <- temp[[1]] #reference
        b1 <- temp[[2]] #query
        
        if (sum(b1>1)<=min_sv_size & sum(a1>2)>=min_sv_size & sum(a1==0)<=min_sv_size) {
          dum1$SV <- "duplication"
          align_good <- dum1
          foo1 <- range(which(a1>1))
          start <- min(dum1[,c("sstart", "send")])+foo1[1]
          end <- min(dum1[,c("sstart", "send")])+foo1[2]
          seqname <- unique(dum1$sseqid)
          read <- unique(dum1$qseqid)
          report_good <- c(seqname, start, abs(start-end), read, "duplication")
          
          return(list(align_good, report_good))
        } else if (sum(a1==0)<=min_sv_size & sum(b1>=2)>=min_sv_size) {
          #repeats in reference
          return(NULL)
        } else {
          #dum1$SV <- "other"
          return(NULL)
        }
      }
      #result <- list(insertion=insertion, deletion=deletion,inversion=inversion, duplication=duplication, tandem=tandem_repeats, translocation=translocation, other=other)
      #return(dum1)
    }, mc.cores=min(length(aligned_list), node))
    #browser() 
    good_num <- sum(sapply(temp_result, function(Y) length(Y[[1]]))!=0)
    print(paste("Found", good_num, "SV reads."))
    return(temp_result)

}



########################################################
#Take out function to call SV for two segment on one chr
########################################################

call_SV_internal <- function(dum1, min_sv_size, indels) {
  temp <- convert_alignment_new(dum1)
  a1 <- temp[[1]] #reference
  b1 <- temp[[2]] #query
  
  #if(sum(a1!=1)<min_sv_size & sum(b1!=1)<min_sv_size) {
  #  return(NULL)
  #}
  
  if (sum(a1==0)<=min_sv_size & sum(b1==0)<=min_sv_size & sum(a1>1)>=min_sv_size &sum(b1>1)<=min_sv_size) {
    #  if (sum(a1==0)<=min_sv_size & sum(b1==0)<=min_sv_size & sum(a1>1)>=min_sv_size) {
    #likely miss very large duplications, due to read length
    dum1$SV <- "duplication"
    align_good <- dum1
    foo1 <- range(which(a1>1))
    start <- min(dum1[,c("sstart", "send")])+foo1[1]
    end <- min(dum1[,c("sstart", "send")])+foo1[2]
    seqname <- unique(dum1$sseqid)
    read <- unique(dum1$qseqid)
    report_good <- c(seqname, start, abs(start-end), read, "duplication")
    
    return(list(align_good, t(report_good)))
  } else if (
            (sum(a1==0)<=min_sv_size & sum(a1>1) <= min_sv_size & sum(b1==0)>=min_sv_size&sum(b1>1)<=min_sv_size) |
            (sum(a1==0)>=min_sv_size & sum(a1>1) <= min_sv_size & sum(b1==0)>=min_sv_size&sum(b1>1)<=min_sv_size&sum(b1==0)/sum(a1==0)>2)
             ) {
    #} else if (sum(a1==0)<sum(b1==0) & sum(b1==0)>=min_sv_size) {
    dum1$SV <- "insertion"
    align_good <- dum1
    
    start <- round(mean(as.numeric(dum1$send[1]), as.numeric(dum1$sstart[2])))
    seqname <- unique(dum1$sseqid)
    read <- unique(dum1$qseqid)
    report_good <- c(seqname, start, abs(dum1$qstart[2]-dum1$qend[1]), read, "insertion")
    return(list(align_good, t(report_good)))
  } else if((sum(b1==0)<=min_sv_size & sum(b1>1) <=min_sv_size & sum(a1==0)>=min_sv_size & sum(a1>1)<=min_sv_size) |
            (sum(a1==0)>=min_sv_size & sum(a1>1) <= min_sv_size & sum(b1==0)>=min_sv_size&sum(b1>1)<=min_sv_size&sum(a1==0)/sum(b1==0)>2)) {
    #} else if(sum(b1==0)<sum(a1==0) & sum(a1==0)>=min_sv_size) {
    dum1$SV <- "deletion"
    align_good <- dum1
    if(unique(dum1$sstrand)=="pos") {
      start <- dum1$send[1]
      end <- dum1$sstart[2]
    } else {
      end <- dum1$send[1]
      start <- dum1$sstart[2]
    }
    seqname <- unique(dum1$sseqid)
    read <- unique(dum1$qseqid)
    report_good <- c(seqname, start, abs(start-end), read, "deletion")
    return(list(align_good, t(report_good)))
  } else if (sum(a1==0)>=min_sv_size & sum(a1>1) <= min_sv_size & sum(b1==0)>=min_sv_size&sum(b1>1)<=min_sv_size&sum(a1==0)/sum(b1==0)<2&sum(a1==0)/sum(b1==0)>0.5) {
    dum1$SV <- "CSUB"
    align_good <- dum1
    if(unique(dum1$sstrand)=="pos") {
      start <- dum1$send[1]
      end <- dum1$sstart[2]
    } else {
      end <- dum1$send[1]
      start <- dum1$sstart[2]
    }
    seqname <- unique(dum1$sseqid)
    read <- unique(dum1$qseqid)
    report_good <- c(seqname, start, abs(start-end), read, "CSUB")
    return(list(align_good, t(report_good)))
  } else {
    #dum1$SV <- "BP_other"
    return(NULL)
    #ind <- which(indels[,4]==unique(dum1$qseqid))
    #if(length(ind)>0) {
    #  dum2 <- list()
    #  for(i in 1:length(ind)) {
    #    dum3 <- dum1
    #    dum3$SV <- indels[ind[i],5]
    #    dum2[[i]] <- dum3
    #  }
      #print(indels[good,])
    #  return(list(do.call("rbind", dum2), as.character(indels[ind,])))
    #}
  }
}


#################################
#Take care of blast alignment more than 3 lines, iterate each pairs to get SV calls nearby
#################################
find_sv_hs_blastn_parallel_v4 <- function(FASTQ, candidates= NULL, min_sv_size = 30, max_dist=10000, blast_ref="Z:/genome/SacCer3/SacCer3_blast", masker_db=NULL, block_size=1000, platform="PacBio CCS", graph=T, node=8) {
  #install.packages('rBLAST', repos = 'https://mhahsler.r-universe.dev')
  library(ShortRead)
  library(parallel)
  #library(tictoc)
  #browser()
  #Sys.setenv(PATH=paste("/opt/hs-blastn/0.0.5/bin/", Sys.getenv("PATH"), sep=":"))
  #tic()
  indels <- read.csv("indels.csv")
  indels <- indels[,-1]
  #colnames(indels) <- ""
  #candidates=indels[,4]
  
  fs <- FastqStreamer(FASTQ, n=block_size)
  aligned_result <- list()
  count <- 0
  streamer_count <-  0
  temp_result <- list()
  print("Start detection process")
  while(length(fq <- yield(fs))) {
    #browser()
    streamer_count <- streamer_count+1
    print(paste("Imported", streamer_count*block_size, "reads."))
    temp_id <- sapply(strsplit((as.character(fq@id)), split=" "), function(X) X[[1]][1])
    ind <- which(is.element(temp_id, candidates)==T)
    if(length(ind)==0) {
      next
    }
    count <- count+length(ind)
    print(paste("Processing", length(ind), "candidates..."))
    a <- DNAStringSet(fq@sread[ind])
    names(a) <- temp_id[ind]
    
    width <- width(a)
    names(width) <- names(a)
    writeXStringSet(a, file="temp.fa")
    #aligned <- predict(ref, a, BLAST_args = paste("-task megablast -window_masker_db /net/nfs-irwrsrchnas01/labs/xwu/genome/hg38/hg38.fa.counts.obinary -num_threads", node), custom_format = "qseqid qlen qstart qend sseqid sstart send sstrand length pident bitscore")
    if(is.null(masker_db)) {
      #this will be too slow
      system(paste("hs-blastn align -query temp.fa -db", blast_ref, "-outfmt 6 -evalue 1e-10 -max_target_seqs 20 -num_threads", node, "> temp.out"))
    } else {
      #Fastest megablast
      system(paste("hs-blastn align -query temp.fa -db", blast_ref, "-window_masker_db", masker_db, "-max_target_seqs 20 -outfmt 6 -evalue 1e-10 -num_threads", node, "> temp.out"), ignore.stderr = T)
    }
    if(file.size("temp.out")) {
      aligned <- read.delim("temp.out", header=F)
    } else {
      next
    }
    aligned$qlen <- width[match(aligned[,1], names(width))]
    colnames(aligned) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen")
    write.csv(aligned, paste0("aligned_", streamer_count, ".csv"))
    
    if(platform=="PacBio CCS") {
      #use higher pindet for ccs data, remove chrM by default
      aligned <- aligned[(aligned$qend-aligned$qstart)>=min_sv_size&aligned$pident>=95&aligned$sseqid!="chrM",]
    } else {
      #allow more errors with non-CCS data, e.g. Nanopore or PacBio CLS
      aligned <- aligned[(aligned$qend-aligned$qstart)>=min_sv_size&aligned$pident>=87&aligned$sseqid!="chrM",]
    }
    aligned[,c("qstart", "qend")] <- t(apply(aligned[,c("qstart", "qend")], 1, function(X) sort(as.numeric(X), decreasing = F)))
    aligned <- aligned[order(aligned$qseqid, aligned$qstart),]
    aligned_list <- split(aligned, f=aligned$qseqid)
    print("Finished alignment with Blast.")
    
    temp_result[[streamer_count]] <- mclapply(1:length(aligned_list), function(i) {
      dum <- aligned_list[[i]]
      dum$sstrand <- apply(dum[,c("sstart", "send")], 1, function(X) ifelse(X[1]< X[2], "pos", "neg"))
      dum$cov <- (dum$qend - dum$qstart +1)/dum$qlen
      
      #Need to refine the subregion to make sure there is no redundant alignments
      if(nrow(dum)>2) {
        dum_temp <- dum
        dum_temp[,c("qstart", "qend")] <- t(apply(dum_temp[,c("qstart", "qend")], 1, function(X) sort(as.numeric(X), decreasing = F)))
        x <- GRanges(seqnames=unique(dum_temp$qseqid), ranges = IRanges(start=dum_temp$qstart, end=dum_temp$qend))
        y <- reduce(x-round(min_sv_size/2))
        foo <- as.data.frame(findOverlaps(x,y, minoverlap = min_sv_size/2))
        segmentation <- list()
        for(j in unique(foo$subjectHits)) {
          segmentation_temp <- dum_temp[foo$queryHits[which(foo$subjectHits==j)],]
          segmentation[[j]] <- segmentation_temp[which(segmentation_temp$bitscore==max(segmentation_temp$bitscore))[1],]
        }
        dum <- do.call("rbind", segmentation)
        #updated segmentation step to remove redundant, not as good
        #y <- findOverlaps(x,x+min_sv_size/2, type="within") 
        #a1 <- as.matrix(y)
        #b1 <- a1[which(a1[,1]!=a1[,2]),]
        #dum <- x[-unique(b1[,1]),]
      }
      best <- which(dum$bitscore==max(dum$bitscore))[1]
      if((max(dum$cov[best])>=0.99 & max(dum$pident[best])>=97) | length(unique(dum$sseqid))>2) {
        #perfect align or too complicated align
        return(NULL)
      }
      
      if(length(unique(dum$sseqid))==2) {
        #browser()
        #1. Check if one chrom already covers entire read, go to following step for segmentation
        foo <- split(dum, f=dum$sseqid)
        each_list <- list()
        for(z in 1:length(foo)) {
          each <- array(0, dim=c(unique(foo[[z]]$qlen), nrow(foo[[z]])))
          for(j in 1:nrow(foo[[z]])) {
            each[foo[[z]][j, "qstart"]:foo[[z]][j, "qend"],j] <- 1
          }
          each_list[[z]] <- each
        }
        cov_list <- sapply(each_list, function(each) sum(apply(each, 1, sum)>0)/nrow(each))
        if(max(cov_list)>=0.99) {
          #one chrom covers entire read, so just use the chrom for next step and ignore the other chrom
          dum <- foo[[which(cov_list>= 0.99)]]
        } else {
          #Check if fragments on both chrom are contiguous and covers most of the read
          each_combined <- do.call("cbind", each_list)
          cov_combined <- sum(apply(each_combined, 1, sum)>0)/nrow(each_combined)
          if(cov_combined>= 0.99) {
            #translocation
            dum$SV <- "translocation"
            align_good <- dum
            dum <- dum[order(dum$sseqid),]
            count<- count+1
            if(dum$sstrand[1]=="pos") {
              start <- dum$send[1]
            } else {
              start <- dum$sstart[1]
            }  
            
            if(dum$sstrand[2]=="pos") {
              end <- dum$sstart[2]
            } else {
              end <- dum$sstart[2]
            }
            seqname <- paste(unique(dum$sseqid), collapse="-")
            read <- unique(dum$qseqid)
            report_good <- c(seqname, start, paste(start, end, sep=":"), read, "translocation")
            return(list(align_good, t(report_good)))
          } else {
            #get the ones aligned to one chromosome with highest bitscore
            dum <- foo[[which(cov_list==max(cov_list))[1]]]
          }
        }
      }
      
      #add in segmentation step to resolve multiple reports involving distant regions
      dum_temp <- dum
      dum_temp[,c("sstart", "send")] <- t(apply(dum_temp[,c("sstart", "send")], 1, function(X) sort(as.numeric(X), decreasing = F)))
      x <- GRanges(seqnames=unique(dum_temp$sseqid), ranges = IRanges(start=dum_temp$sstart, end=dum_temp$send))
      y <- reduce(x+max_dist)
      foo <- as.data.frame(findOverlaps(x, y))
      segmentation <- array(dim=length(unique(foo$subjectHits)))
      for(j in unique(foo$subjectHits)) {
        segmentation_temp <- x[foo$queryHits[which(foo$subjectHits==j)]]
        segmentation[j] <- sum(width(segmentation_temp))
      }
      
      idx <- which(segmentation==max(segmentation)[1])
      dum1 <- dum[foo$queryHits[which(foo$subjectHits==idx)],]
      
      #Now checking for SV 
      if(nrow(dum1)==1) {
        #note that some reads have multiple entries reported from parsing pbmm2
        ind <- which(indels[,4]==dum1$qseqid)
        if(length(ind)>0) {
          dum2 <- list()
          for(z in 1:length(ind)) {
            dum3 <- dum1
            dum3$SV <- indels[ind[z],5]
            dum2[[z]] <- dum3
          }
          return(list(do.call("rbind", dum2), indels[ind,]))
        }
        #} else if  (nrow(dum1)>2 &nrow(dum) <=4 & length(unique(dum1$sstrand))==2) {
      } else if  (nrow(dum1)==3 & (sum(dum1$sstrand==c("pos", "neg", "pos"))==3 | sum(dum1$sstrand==c("neg", "pos", "neg"))==3)) {
        temp <- convert_alignment_new(dum1)
        a1 <- temp[[1]] #reference
        b1 <- temp[[2]] #query
        
        if (sum(b1==0)<=min_sv_size & sum(a1==0)<=min_sv_size & sum(b1>1)<=min_sv_size & sum(a1>1)>=min_sv_size) {
          #print("IVT")
          dum1$SV <- "IVT"
          align_good <- dum1
          
          if(dum1$sstrand[2]=="pos") {
            start <- dum1$sstart[2]
            end <- dum1$send[2]
          } else {
            end <- dum1$sstart[2]
            start <- dum1$send[2]
          }
          seqname <- unique(dum1$sseqid)
          read <- unique(dum1$qseqid)
          report_good <- c(seqname, start, abs(start-end), read, "IVT")
          
          return(list(align_good, t(report_good)))
        } else {
          dum1$SV <- "inversion"
          align_good <- dum1
          
          if(dum1$sstrand[2]=="pos") {
            start <- dum1$sstart[2]
            end <- dum1$send[2]
          } else {
            end <- dum1$sstart[2]
            start <- dum1$send[2]
          }
          seqname <- unique(dum1$sseqid)
          read <- unique(dum1$qseqid)
          report_good <- c(seqname, start, abs(start-end), read, "inversion")
          
          return(list(align_good, t(report_good)))
        }
      } else if  (nrow(dum1)==2 & length(unique(dum1$sstrand))==2) {
        temp <- convert_alignment_new(dum1)
        a1 <- temp[[1]] #reference
        b1 <- temp[[2]] #query
        
        if (sum(b1==0)<=min_sv_size & sum(a1==0)<=min_sv_size & sum(b1>1)<=min_sv_size) {
          dum1$SV <- "inversion_other"
          align_good <- dum1
          
          if(dum$sstrand[1]=="pos") {
            start <- dum1$send[1]
          } else {
            start <- dum1$sstart[1]
          }  
          
          if(dum1$sstrand[2]=="pos") {
            end <- dum1$sstart[2]
          } else {
            end <- dum1$sstart[2]
          }
          
          seqname <- unique(dum1$sseqid)
          read <- unique(dum1$qseqid)
          report_good <- c(seqname, start, abs(start-end), read, "inversion_other")
          #return(NULL)
          return(list(align_good, t(report_good)))
        }
      } else if(nrow(dum1)==2 & length(unique(dum1$sstrand))==1) {
          sv_internal <- call_SV_internal(dum1, min_sv_size, indels)
          if(is.null(sv_internal)) {
            ind <- which(indels[,4]==unique(dum1$qseqid))
            if(length(ind)>0) {
              dum2 <- list()
              for(z in 1:length(ind)) {
                dum3 <- dum1
                dum3$SV <- indels[ind[z],5]
                dum2[[z]] <- dum3
              }
              #print(indels[good,])
              return(list(do.call("rbind", dum2),indels[ind,]))
            }
          } else {
            return(sv_internal)
          }
          
      } else if(nrow(dum1) >=3 & length(unique(dum1$sstrand))==1) {
        temp <- convert_alignment_new(dum1)
        a1 <- temp[[1]] #reference
        b1 <- temp[[2]] #query
        
        if (sum(a1>1)>=min_sv_size & sum(a1==0)<=min_sv_size) {
        #if (sum(b1==0)<=min_sv_size & sum(a1>1)>=min_sv_size & sum(a1==0)<=min_sv_size) {
          dum1$SV <- "duplication"
          align_good <- dum1
          foo1 <- range(which(a1>1))
          start <- min(dum1[,c("sstart", "send")])+foo1[1]
          end <- min(dum1[,c("sstart", "send")])+foo1[2]
          seqname <- unique(dum1$sseqid)
          read <- unique(dum1$qseqid)
          report_good <- c(seqname, start, abs(start-end), read, "duplication")
          
          return(list(align_good, t(report_good)))
        }  else {
          #print("check")
          sv_internal <- list()
          for(z in 1:(nrow(dum1)-1)) {
            dum2 <- dum1[z:(z+1),]
            sv_internal[[z]] <- call_SV_internal(dum2, min_sv_size, indels)
          }
          #print(sv_internal)
          #return(list(lapply(sv_internal, function(X) X[[1]]), sapply(sv_internal, function(X) X[[2]])))
          if(sum(sapply(sv_internal, is.null))>0) {
            return(list(do.call("rbind", lapply(sv_internal, function(X) X[[1]])), 
                        do.call("rbind", lapply(sv_internal, function(X) X[[2]]))))
          } else {
            ind <- which(indels[,4]==unique(dum1$qseqid))
            if(length(ind)>0) {
              dum2 <- list()
              for(z in 1:length(ind)) {
                dum3 <- dum1
                dum3$SV <- indels[ind[z],5]
                dum2[[z]] <- dum3
              }
              #print(indels[good,])
              return(list(do.call("rbind", dum2), indels[ind,]))
            }
          }
          
        }
      } else {
        ind <- which(indels[,4]==unique(dum1$qseqid))
        if(length(ind)>0) {
          dum2 <- list()
          for(z in 1:length(ind)) {
            dum3 <- dum1
            dum3$SV <- indels[ind[z],5]
            dum2[[z]] <- dum3
          }
          #print(indels[good,])
          return(list(do.call("rbind", dum2), indels[ind,]))
        }
      }
      #result <- list(insertion=insertion, deletion=deletion,inversion=inversion, duplication=duplication, tandem=tandem_repeats, translocation=translocation, other=other)
      #return(dum1)
    }, mc.cores=node)
    
    good_num <- sum(sapply(temp_result,function(X) sum(sapply(X, function(Y) length(Y[[1]]))!=0)))
    print(paste("Found", good_num, "SV reads so far."))
    gc()
  }
  #toc()
  close(fs)
  saveRDS(temp_result, file="temp_result.rds")
  #browser()
  good <- lapply(temp_result, function(X) X[sapply(X, function(Y) !is.null(Y)&sum(is.na(Y))==0&class(Y)!="try-error")])
  result <- do.call("rbind", lapply(good, function(X) do.call("rbind", sapply(X, function(Y) Y[1]))))
  a <- split(result, f=result$SV)
  result_list <- lapply(a, function(X) split(X, f=X$qseqid))
  print(sapply(result_list, length))
  #saveRDS(result_list, file="sv_all.rds")
  #report <- do.call("rbind", lapply(good, function(X) do.call("rbind", sapply(X, function(Y) Y[2]))))
  report <- lapply(good, function(X) t(sapply(X, function(Y) as.data.frame(Y[[2]]))))
  report_filtered <- report[!sapply(report, is.null)]
  for(z in 1:length(report_filtered)){
    colnames(report_filtered[[z]]) <- colnames(indels)
  }
  report <- do.call("rbind", report_filtered)
  #write.csv(report, file="sv_report.csv")
  return(list(result_list, report_filtered))
}



#################################
#This is the main pipeline to for 
#long read SV detection
#################################
call_SV_pacbio <- function(BAM, fastq, candidates=NULL, ref=NULL, max_dist = 10000, min_sv_length=30, window_masker_db=NULL, block_size=1000, platform="PacBio CCS", graph=T, node=8, out="SV_result") {
  library(rtracklayer)
  library(tictoc)
  
  if(is.null(candidates)) {
    candidates <- get_candidate_v3(BAM=BAM, softclip_length = min_sv_length)
    writeLines(candidates, "candidates.txt")
  } else {
    candidates=readLines(candidates)
  }
  
  
  #start_time <- Sys.time()
  tic()
  sv_result <- find_sv_hs_blastn_parallel_v4(FASTQ = fastq,  blast_ref=ref, masker_db=window_masker_db, candidates=candidates,max_dist=max_dist, min_sv_size=min_sv_length, block_size=block_size, platform=platform, node=node)
  #end_time <- Sys.time()
  #print(paste(start_time, end_time))
  #print(paste("Total run time:", end_time - start_time))
  toc()
  
  saveRDS(sv_result, file="sv_all.rds")
  system(paste0("mkdir ", out))
  setwd(out)
  sv <- sv_result[[1]]
  for(i in names(sv)) {
    write.csv(do.call("rbind", sv[[i]]), file=paste0(i, ".csv"))
  }
  if(is.list(sv_result[[2]])) {
    dum <- lapply(sv_result[[2]], function(X) apply(X, 1, function(Y) do.call("cbind", Y)))
    report <- do.call("rbind", lapply(dum, function(X) do.call("rbind", X)))
  } else {
    report <- sv_result[[2]]
  }
  
  write.csv(report, "sv_report_all.csv")
  SV_collapsed <- Collapse_SV_v2(report)
  write.csv(SV_collapsed, "SV_report_summary.csv")
  
  if(graph) {
    print("Generating graphs...")
    for(i in names(sv)) {
      system(paste0("mkdir ", i))
      setwd(i)
      if(i=="translocation") {
        generate_image_translocation(sv[["translocation"]], output=T, cores=node)
      } else {
        generate_image(sv[[i]], output=T, cores=node)
      }
      setwd("../")
    }
  }  
  
  return(sv_result)
}



call_SV_pacbio_chunk <- function(BAM, fastq, candidates=NULL, ref=NULL, max_dist = 10000, min_sv_length=30, window_masker_db=NULL, block_size=10000, platform="PacBio CCS", graph=T, node=8, out="SV_result") {
  library(rtracklayer)
  library(Rsamtools)
  library(tictoc)
  
  print(paste("Start SV detection with", BAM))
  tic()
  sv_result_list <- list()
  candidates_list <- list()
  if(is.null(candidates)) {
    BAM_chunk <- BamFile(BAM, yieldSize=block_size)
    open(BAM_chunk)
    
    total_alignment <- 0
    reads <- scanBam(BAM_chunk)
    reads <- reads[[1]]
    #reads <- reads[which(reads$qwidth>=2500)]
    sequences <- reads$seq
    names(sequences) <- reads$qname
    count <- 1
    while(length(reads$qname)) {
      #browser()
      total_alignment <- total_alignment+length(reads$qname)
      print(paste("Imported", total_alignment, "alignments."))
      candidates_list[[count]] <- get_candidate_v5(read =reads, softclip_length = 30, chunk=block_size)
      sv_result_list[[count]] <- find_sv_hs_blastn_parallel_v3(FASTQ = sequences,  
                                                               blast_ref=ref, 
                                                               masker_db=window_masker_db, 
                                                               candidates=candidates_list[[count]], 
                                                               max_dist=max_dist, 
                                                               min_sv_size=min_sv_length,
                                                               block_size=block_size, 
                                                               platform=platform, 
                                                               node=node)
      
      reads <- scanBam(BAM_chunk)
      reads <- reads[[1]]
      sequences <- reads$seq
      names(sequences) <- reads$qname
      count <- count + 1
    }
    candidates <- unlist(candidates_list)
    writeLines(candidates, "candidates.txt")

    good <- lapply(sv_result_list, function(X) X[sapply(X, function(Y) !is.null(Y)&sum(is.na(Y))==0&class(Y)!="try-error")])
    result <- do.call("rbind", lapply(good, function(X) do.call("rbind", sapply(X, function(Y) Y[1]))))
    a <- split(result, f=result$SV)
    result_list <- lapply(a, function(X) split(X, f=X$qseqid))
    print(sapply(result_list, length))
    report <- do.call("rbind", lapply(good, function(X) do.call("rbind", sapply(X, function(Y) Y[2]))))
    sv_result <- list(result_list, report)
    
    
  } else {
    candidates=readLines(candidates)
    sv_result <- find_sv_hs_blastn_parallel_v2(FASTQ = fastq,  
                                               blast_ref=ref, 
                                               masker_db=window_masker_db, 
                                               candidates=candidates,
                                               max_dist=max_dist, 
                                               block_size=block_size, 
                                               platform=platform, 
                                               node=node)
  }
  

  close(BAM_chunk)
  toc()
  #browser()
  
  system(paste0("mkdir ", out))
  setwd(out)
  saveRDS(sv_result, file="sv_all.rds")
  sv <- sv_result[[1]]
  for(i in names(sv)) {
    write.csv(do.call("rbind", sv[[i]]), file=paste0(i, ".csv"))
  }
  
  report <- sv_result[[2]]
  write.csv(report, "sv_report_all.csv")
  SV_collapsed <- Collapse_SV_v2(report)
  write.csv(SV_collapsed, "SV_report_summary.csv")
  
  if(graph) {
    print("Generating graphs...")
    for(i in names(sv)) {
      system(paste0("mkdir ", i))
      setwd(i)
      if(i=="translocation") {
        generate_image_translocation(sv[["translocation"]], output=T, cores=node)
      } else {
        generate_image(sv[[i]], output=T, cores=node)
      }
      setwd("../")
    }
  }  
  
  return(sv_result)
}


############################
#It is a helper function to check 
#any given read 
#for its alignment and graph 
#generation
############################
check_read <- function(fastq, ID, blast_ref, min_sv_size=50) {
  if(is.null(ID)) {
    foo <- readFastq(fastq)
  } else {
    system(paste("grep -A 3", ID, fastq, "> test.fastq"))
    foo <- readFastq("test.fastq")
  }
  print("found read")
  a <- foo@sread
  names(a) <- foo@id[1]
  #ref <- blast(blast_ref)
  #aligned <- predict(blast_ref, a, custom_format = "qseqid qlen qstart qend sseqid sstart send sstrand length pident bitscore")
  system(paste("hs-blastn align -query test.fastq -db", blast_ref, "-window_masker_db", paste0(blast_ref, ".counts.obinary"), "-outfmt 6 -evalue 1e-10 -max_target_seqs 20 -gapopen 3 -gapextend 3 -reward 1 -penalty -5 > temp.out"))
  aligned <- read.delim("temp.out", header=F)
  aligned$qlen <- width(a)
  colnames(aligned) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen")
  print(aligned)
  aligned <- aligned[(aligned$qend-aligned$qstart)>=50&aligned$pident>=94,]
  aligned[,c("qstart", "qend")] <- t(apply(aligned[,c("qstart", "qend")], 1, function(X) sort(as.numeric(X), decreasing = F)))
  aligned <- aligned[order(aligned$qseqid, aligned$qstart),]
  
  dum <- aligned
  dum$sstrand <- apply(dum[,c("sstart", "send")], 1, function(X) ifelse(X[1]< X[2], "pos", "neg"))
  dum$cov <- (dum$qend - dum$qstart +1)/dum$qlen
  print(dum)
  #generate_image(list(dum), output=T)
  #Need to refine the subregion to make sure there is no redundant alignments
  if(nrow(dum)>2) {
    dum_temp <- dum
    dum_temp[,c("qstart", "qend")] <- t(apply(dum_temp[,c("qstart", "qend")], 1, function(X) sort(as.numeric(X), decreasing = F)))
    x <- GRanges(seqnames=unique(dum_temp$qseqid), ranges = IRanges(start=dum_temp$qstart, end=dum_temp$qend))
    y <- reduce(x-round(min_sv_size/2))
    #y <- reduce(x, min.gapwidth=round(min_sv_size/2))
    browser()
    foo <- as.data.frame(findOverlaps(x,y, minoverlap = min_sv_size/2))
    segmentation <- list()
    for(j in unique(foo$subjectHits)) {
      segmentation_temp <- dum_temp[foo$queryHits[which(foo$subjectHits==j)],]
      segmentation[[j]] <- segmentation_temp[which(segmentation_temp$bitscore==max(segmentation_temp$bitscore)),]
    }
    dum1 <- do.call("rbind", segmentation)
  }
  print(dum1)
  if(length(unique(dum1$sseqid))==1) {
    generate_image(list(dum1), output=T)
  } else {
    generate_image_translocation(list(dum1), output=T)
  }
  return(dum1)
}


eval_sim <- function(truth, result, gap=500) {
  colnames(truth)[1:2] <- c("seqname", "start")
  truth[,1] <- paste0("chr", truth[,1])
  truth$SVLENGTH <- gsub("x", "*",truth$SVLENGTH)
  for(i in 1:nrow(truth)) {
    truth$SVLENGTH[i] <- eval(parse(text=truth$SVLENGTH[i]))
  }

  truth$end <- truth$start+as.numeric(as.character(truth$SVLENGTH))
  truth_range <- GRanges(truth)
  
  colnames(result)[2:3] <- c("seqname", "start")
  result$end <- result$start+result$Length
  result_range <- GRanges(result)
  
  overlaps <- array(dim=c(length(unique(truth$TYPE)), length(unique(result$sv_type))))
  missed <- list()
  for(i in 1:length(unique(truth$TYPE))) {
    for(j in 1:length(unique(result$sv_type))) {
      #print(paste(i,j))
      dum <- findOverlaps(result_range[result_range$sv_type==unique(result$sv_type)[j]], 
                          truth_range[truth_range$TYPE==unique(truth$TYPE)[i]], maxgap = gap)
      overlaps[i,j] <- length(unique(as.matrix(dum)[,2]))
      if(i==4) {
        #missed[[j]] <- result_range[result_range$sv_type==unique(result$sv_type)[j]][unique(as.matrix(dum)[,1])]
        missed[[j]] <- truth_range[truth_range$TYPE==unique(truth$TYPE)[i]][unique(as.matrix(dum)[,2])]
      }
       
    }
  }
  colnames(overlaps) <- unique(result$sv_type)
  row.names(overlaps) <- unique(truth$TYPE)
  print(overlaps)
  return(missed)
}
