get_candidate <- function(reads=NULL, BAM, softclip_length) {
  library(ShortRead)
  if(is.null(reads)){
    print("Importing bam file...")
    read <- readGAlignments(BAM, use.names = T)
    print("Done importing bam.")
  } else {
    read <- reads
  }

  #Reads with large insertion or softclip
  dum <- read@cigar
  a <- grep("^[0-9]*S", dum)
  bad1 <- a[which(sapply(strsplit(dum[a], "S"), function(X) as.integer(X[1]))>=softclip_length)]
  b <- grep("[0-9]*S$", dum)
  b1 <- gsub("[0-9]*[IDX=]", "", dum[b])
  bad2 <- b[which(sapply(strsplit(b1, split="S"), function(X) {
    temp <- as.integer(X)
    temp[length(temp)]
  })>=softclip_length)]
  #Also include large insertion in the middle
  c <- grep("I", dum)
  c1 <- gsub("[0-9]*[SDX=]", "", dum[c])
  bad3 <- c[sapply(strsplit(c1, split="I"), function(X) {
    temp <- as.integer(X)
    good <- FALSE
    if(max(temp)>=softclip_length) {
      good <- TRUE
    }
    return(good)
  })]
  
  #Also include large deletion in the middle
  d <- grep("D", dum)
  d1 <- gsub("[0-9]*[SIX=]", "", dum[d])
  bad4 <- d[sapply(strsplit(d1, split="D"), function(X) {
    temp <- as.integer(X)
    good <- FALSE
    if(max(temp)>=softclip_length) {
      good <- TRUE
    }
    return(good)
  })]
  #browser()
  bad <- unique(c(bad1, bad2, bad3, bad4))
  #browser()
  print(paste("There are", length(bad), "reads with long softclips or large indels."))
  return(names(read)[bad])
}



find_sv_blast <- function(FASTQ, candidates= NULL, min_sv_size = 50, max_dist=10000, blast_ref="Z:/genome/SacCer3/SacCer3_blast", block_size=1000, exclude_chrM=T, graph=T) {
  #install.packages('rBLAST', repos = 'https://mhahsler.r-universe.dev')
  library(ShortRead)
  library(rBLAST)
  library(parallel)
  #Please make sure blastn is in PATH. You can set the parameters with the following command
  dum <- Sys.info()
  if(dum[1] == "Windows") {
    Sys.setenv(PATH=paste(Sys.which("blastn"), Sys.getenv("PATH"), sep=";"))
  } else if(dum[1]== "Linux") {
    Sys.setenv(PATH=paste(Sys.which("blastn"), Sys.getenv("PATH"), sep=":"))
  }
  ref <- blast(blast_ref)
  fs <- FastqStreamer(FASTQ, n=block_size)
  aligned_result <- list()
  count <- 0
  streamer_count <-  0
  insertion <- list()
  deletion <- list()
  inversion <- list()
  duplication <- list()
  tandem_repeats <- list()
  other <- list()
  print("Start detection process")
  while(length(fq <- yield(fs))) {
    streamer_count <- streamer_count+1
    print(paste("Imported", streamer_count*block_size, "reads."))
    ind <- which(is.element(fq@id, candidates)==T)
    count <- count+length(ind)
    print(paste("Processing", count, "candidates..."))
    a <- DNAStringSet(fq@sread[ind])
    names(a) <- fq@id[ind]
    aligned <- predict(ref, a, custom_format = "qseqid qlen qstart qend sseqid sstart send sstrand length pident bitscore")
    if(exclude_chrM) {
      aligned <- aligned[aligned$bitscore>=500&aligned$pident>=95&aligned$sseqid!="chrM",]
    }
    aligned[,c("qstart", "qend")] <- t(apply(aligned[,c("qstart", "qend")], 1, function(X) sort(as.numeric(X), decreasing = F)))
    aligned <- aligned[order(aligned$qseqid, aligned$qstart),]
    aligned_list <- split(aligned, f=aligned$qseqid)
    print("Finished alignment with Blast.")
    for(i in 1:length(aligned_list)) {
      #
      if(i %% 100 ==0) {
        print(i)
      }

      dum <- aligned_list[[i]]
      dum$cov <- (dum$qend - dum$qstart +1)/dum$qlen
      if(max(dum$cov[which(dum$bitscore==max(dum$bitscore))])>=0.95 |
         length(unique(dum$sseqid))>1) {
        next
      }

      #add in segmenation step to resolve multiple reports involving distant regions
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

      if(nrow(dum1)==1) {
        #normal
        next
      } else if(nrow(dum1)==3) {
        if(length(unique(dum1$sstrand))==2) {
          inversion[[names(aligned_list)[i]]] <- dum1
          next
        } else {
          other[[names(aligned_list)[i]]] <- dum1
          next
        }
      } else if  (nrow(dum1)==2 & length(unique(dum1$sstrand))==2) {
        inversion[[names(aligned_list)[i]]] <- dum1
        next
      } else if(nrow(dum1)==2 & length(unique(dum1$sstrand))==1) {
        temp <- convert_alignment(dum1)

        #reference
        #a1 <- rowSums(temp)
        a1 <- rowSums(abs(temp))
        #query
        #b1 <- colSums(temp)
        b1 <- colSums(abs(temp))

        if(sum(a1==0)==0 & sum(b1==0)==0) {
          next
        }

        if (sum(b1==0)<=min_sv_size & sum(a1>1)>=min_sv_size) {
          duplication[[names(aligned_list)[i]]] <- dum1
        } else if (sum(a1==0)<=min_sv_size & sum(a1>1) <= min_sv_size & sum(b1==0)>=min_sv_size) {
          insertion[[names(aligned_list)[i]]] <- dum1
        } else if(sum(b1==0)<=min_sv_size & sum(b1>1) <=min_sv_size & sum(a1==0)>=min_sv_size) {
          deletion[[names(aligned_list)[i]]] <- dum1
        } else {
          other[[names(aligned_list)[i]]] <- dum1
        }

      } else if(nrow(dum1) >3) {
        temp <- convert_alignment(dum1)

        #reference
        a1 <- rowSums(abs(temp))
        #query
        b1 <- colSums(abs(temp))
        if (sum(b1==0)<=min_sv_size & sum(a1>2)>=min_sv_size) {
          tandem_repeats[[names(aligned_list)[i]]] <- dum1
          next
        } else if (sum(a1==0)<=min_sv_size & sum(b1>=2)>=min_sv_size) {
          #normal
          next
        } else {
          other[[names(aligned_list)[i]]] <- dum1
          next
        }
      }
    }
    good_num <- length(insertion) + length(deletion) + length(inversion) + length(duplication) + length(tandem_repeats) + length(other)
    print(paste("Found", good_num, "SV events so far."))
    gc()
  }
  close(fs)
  write.csv(do.call("rbind", insertion), file="SV_result/insertion.csv")
  write.csv(do.call("rbind", deletion), file="SV_result/deletion.csv")
  write.csv(do.call("rbind", inversion), file="SV_result/inversion.csv")
  write.csv(do.call("rbind", duplication), file="SV_result/duplication.csv")
  write.csv(do.call("rbind", tandem_repeats), file="SV_result/tandem_repeats.csv")
  write.csv(do.call("rbind", other), file="SV_result/other.csv")

  result <- list(insertion=insertion, deletion=deletion, duplication=duplication, tandem=tandem_repeats, inversion=inversion, other=other)
  print("summary of findings:")
  print(paste("Insertion:", length(insertion)))
  print(paste("Deletion:", length(deletion)))
  print(paste("Duplication:", length(duplication)))
  print(paste("Tandem repeats:", length(tandem_repeats)))
  print(paste("Inversion:", length(inversion)))
  print(paste("Other:", length(other)))

  return(result)
}
###########################
#Now Detect translocation
#06/08/2023
###########################
find_sv_blast_1 <- function(FASTQ, candidates= NULL, min_sv_size = 50, max_dist=10000, blast_ref="Z:/genome/SacCer3/SacCer3_blast", block_size=1000, exclude_chrM=T, graph=T) {
  #install.packages('rBLAST', repos = 'https://mhahsler.r-universe.dev')
  library(ShortRead)
  library(rBLAST)
  library(parallel)
  #Please make sure blastn is in PATH. You can set the parameters with the following command
  dum <- Sys.info()
  if(dum[1] == "Windows") {
    Sys.setenv(PATH=paste(Sys.which("blastn"), Sys.getenv("PATH"), sep=";"))
  } else if(dum[1]== "Linux") {
    Sys.setenv(PATH=paste(Sys.which("blastn"), Sys.getenv("PATH"), sep=":"))
  }
  ref <- blast(blast_ref)
  fs <- FastqStreamer(FASTQ, n=block_size)
  aligned_result <- list()
  count <- 0
  streamer_count <-  0
  insertion <- list()
  deletion <- list()
  inversion <- list()
  duplication <- list()
  tandem_repeats <- list()
  translocation <- list()
  other <- list()
  print("Start detection process")
  while(length(fq <- yield(fs))) {
    streamer_count <- streamer_count+1
    print(paste("Imported", streamer_count*block_size, "reads."))
    ind <- which(is.element(fq@id, candidates)==T)
    count <- count+length(ind)
    print(paste("Processing", count, "candidates..."))
    a <- DNAStringSet(fq@sread[ind])
    names(a) <- fq@id[ind]
    aligned <- predict(ref, a, custom_format = "qseqid qlen qstart qend sseqid sstart send sstrand length pident bitscore")
    if(exclude_chrM) {
      aligned <- aligned[aligned$bitscore>=500&aligned$pident>=95&aligned$sseqid!="chrM",]
    }
    aligned[,c("qstart", "qend")] <- t(apply(aligned[,c("qstart", "qend")], 1, function(X) sort(as.numeric(X), decreasing = F)))
    aligned <- aligned[order(aligned$qseqid, aligned$qstart),]
    aligned_list <- split(aligned, f=aligned$qseqid)
    print("Finished alignment with Blast.")
    #browser()
    for(i in 1:length(aligned_list)) {
      #
      if(i %% 100 ==0) {
        print(i)
      }

      dum <- aligned_list[[i]]
      dum$cov <- (dum$qend - dum$qstart +1)/dum$qlen
      best <- which(dum$bitscore==max(dum$bitscore))

      if((max(dum$cov[best])>=0.99 & max(dum$pident[best])>=0.98) | length(unique(dum$sseqid))>2) {
      #perfect align
          next
      }

      if(length(unique(dum$sseqid))==2) {
      #translocation
        translocation[[names(aligned_list)[i]]] <- dum
        next
      }
      #add in segmenation step to resolve multiple reports involving distant regions
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

      if(nrow(dum1)==1) {
        #normal
        next
      } else if(nrow(dum1)==3) {
        if(length(unique(dum1$sstrand))==2) {
          inversion[[names(aligned_list)[i]]] <- dum1
          next
        } else {
          other[[names(aligned_list)[i]]] <- dum1
          next
        }
      } else if  (nrow(dum1)==2 & length(unique(dum1$sstrand))==2) {
        inversion[[names(aligned_list)[i]]] <- dum1
        next
      } else if(nrow(dum1)==2 & length(unique(dum1$sstrand))==1) {
        temp <- convert_alignment(dum1)

        #reference
        #a1 <- rowSums(temp)
        a1 <- rowSums(abs(temp))
        #query
        #b1 <- colSums(temp)
        b1 <- colSums(abs(temp))

        if(sum(a1==0)==0 & sum(b1==0)==0) {
          next
        }

        if (sum(b1==0)<=min_sv_size & sum(a1>1)>=min_sv_size) {
          duplication[[names(aligned_list)[i]]] <- dum1
        } else if (sum(a1==0)<=min_sv_size & sum(a1>1) <= min_sv_size & sum(b1==0)>=min_sv_size) {
          insertion[[names(aligned_list)[i]]] <- dum1
        } else if(sum(b1==0)<=min_sv_size & sum(b1>1) <=min_sv_size & sum(a1==0)>=min_sv_size) {
          deletion[[names(aligned_list)[i]]] <- dum1
        } else {
          other[[names(aligned_list)[i]]] <- dum1
        }

      } else if(nrow(dum1) >3) {
        temp <- convert_alignment(dum1)

        #reference
        a1 <- rowSums(abs(temp))
        #query
        b1 <- colSums(abs(temp))
        if (sum(b1==0)<=min_sv_size & sum(a1>2)>=min_sv_size) {
          tandem_repeats[[names(aligned_list)[i]]] <- dum1
          next
        } else if (sum(a1==0)<=min_sv_size & sum(b1>=2)>=min_sv_size) {
          #normal
          next
        } else {
          other[[names(aligned_list)[i]]] <- dum1
          next
        }
      }
    }
    good_num <- length(insertion) + length(deletion) + length(inversion) + length(duplication) + length(tandem_repeats) + length(translocation) + length(other)
    print(paste("Found", good_num, "SV events so far."))
    gc()
  }
  close(fs)
  write.csv(do.call("rbind", insertion), file="SV_result/insertion.csv")
  write.csv(do.call("rbind", deletion), file="SV_result/deletion.csv")
  write.csv(do.call("rbind", inversion), file="SV_result/inversion.csv")
  write.csv(do.call("rbind", duplication), file="SV_result/duplication.csv")
  write.csv(do.call("rbind", tandem_repeats), file="SV_result/tandem_repeats.csv")
  write.csv(do.call("rbind", translocation), file="SV_result/translocation.csv")
  write.csv(do.call("rbind", other), file="SV_result/other.csv")

  result <- list(insertion=insertion, deletion=deletion, duplication=duplication, tandem=tandem_repeats, inversion=inversion, other=other)
  print("summary of findings:")
  print(paste("Insertion:", length(insertion)))
  print(paste("Deletion:", length(deletion)))
  print(paste("Duplication:", length(duplication)))
  print(paste("Tandem repeats:", length(tandem_repeats)))
  print(paste("Inversion:", length(inversion)))
  print(paste("translocation:", length(translocation)))
  print(paste("Other:", length(other)))

  return(result)
}

#########################
#06/09/2023
#parallel version
#########################
find_sv_blast_parallel <- function(FASTQ, candidates= NULL, min_sv_size = 50, max_dist=10000, blast_ref="Z:/genome/SacCer3/SacCer3_blast", block_size=1000, exclude_chrM=T, graph=T, node=8) {
  #install.packages('rBLAST', repos = 'https://mhahsler.r-universe.dev')
  library(ShortRead)
  library(rBLAST)
  library(parallel)
  #Please make sure blastn is in PATH. You can set the parameters with the following command
  dum <- Sys.info()
  if(dum[1] == "Windows") {
    Sys.setenv(PATH=paste(Sys.which("blastn"), Sys.getenv("PATH"), sep=";"))
  } else if(dum[1]== "Linux") {
    Sys.setenv(PATH=paste(Sys.which("blastn"), Sys.getenv("PATH"), sep=":"))
  }
  ref <- blast(blast_ref)
  fs <- FastqStreamer(FASTQ, n=block_size)
  aligned_result <- list()
  count <- 0
  streamer_count <-  0
  temp_result <- list()
  print("Start detection process")
  while(length(fq <- yield(fs))) {
    streamer_count <- streamer_count+1
    print(paste("Imported", streamer_count*block_size, "reads."))
    ind <- which(is.element(fq@id, candidates)==T)
	if(length(ind)==0) {
		next
	}
    count <- count+length(ind)
    print(paste("Processing", count, "candidates..."))
    a <- DNAStringSet(fq@sread[ind])
    names(a) <- fq@id[ind]
    aligned <- predict(ref, a, BLAST_args = paste("-task megablast -window_masker_db /net/nfs-irwrsrchnas01/labs/xwu/genome/hg38/hg38.fa.counts.obinary -num_threads", node), custom_format = "qseqid qlen qstart qend sseqid sstart send sstrand length pident bitscore")
    if(exclude_chrM) {
      aligned <- aligned[aligned$bitscore>=200&aligned$pident>=95&aligned$sseqid!="chrM",]
    }
    aligned[,c("qstart", "qend")] <- t(apply(aligned[,c("qstart", "qend")], 1, function(X) sort(as.numeric(X), decreasing = F)))
    aligned <- aligned[order(aligned$qseqid, aligned$qstart),]
    aligned_list <- split(aligned, f=aligned$qseqid)
    print("Finished alignment with Blast.")

    temp_result[[streamer_count]] <- mclapply(1:length(aligned_list), function(i) {


		dum <- aligned_list[[i]]
        dum$cov <- (dum$qend - dum$qstart +1)/dum$qlen
        best <- which(dum$bitscore==max(dum$bitscore))

        if((max(dum$cov[best])>=0.99 & max(dum$pident[best])>=0.98) | length(unique(dum$sseqid))>2) {
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
				#Check if fragments on both chom are contiguous and covers most of the read
				each_combined <- do.call("cbind", each_list)
				cov_combined <- sum(apply(each_combined, 1, sum)>0)/nrow(each_combined)
				if(cov_combined>= 0.99) {
				    #translocation
					dum$SV <- "translocation"
					return(dum)
				} else {
					return(NULL)
				}
			}
        }
        #add in segmenation step to resolve multiple reports involving distant regions
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
          #normal
          return(NULL)
        } else if  (nrow(dum1)<=3 & length(unique(dum1$sstrand))==2) {
          dum1$SV <- "inversion"
		  return(dum1)
        } else if(nrow(dum1)<=3 & length(unique(dum1$sstrand))==1) {
          temp <- convert_alignment(dum1)
          a1 <- rowSums(abs(temp)) #reference
          b1 <- colSums(abs(temp)) #query

          if(sum(a1==0)==0 & sum(b1==0)==0) {
            return(NULL)
          }

          if (sum(b1==0)<=min_sv_size & sum(a1>1)>=min_sv_size &sum(b1>1)<=min_sv_size) {
            dum1$SV <- "duplication"
			return(dum1)
          } else if (sum(a1==0)<=min_sv_size & sum(a1>1) <= min_sv_size & sum(b1==0)>=min_sv_size&sum(b1>1)<=min_sv_size) {
            dum1$SV <- "insertion"
			return(dum1)
          } else if(sum(b1==0)<=min_sv_size & sum(b1>1) <=min_sv_size & sum(a1==0)>=min_sv_size & sum(a1>1)<=min_sv_size) {
            dum1$SV <- "deletion"
			return(dum1)
          } else {
            dum1$SV <- "other"
			return(dum1)
          }

        } else if(nrow(dum1) >3 & length(unique(dum1$sstrand))==1) {
          temp <- convert_alignment(dum1)
          #reference
          a1 <- rowSums(abs(temp))
          #query
          b1 <- colSums(abs(temp))
          if (sum(b1==0)<=min_sv_size & sum(a1>2)>=min_sv_size & sum(a1==0)<=min_sv_size) {
            dum1$SV <- "tandem_repeats"
			#print(dum1)
			return(dum1)
          } else if (sum(a1==0)<=min_sv_size & sum(b1>=2)>=min_sv_size) {
            #repeats in reference
            return(NULL)
          } else {
            dum1$SV <- "other"
			return(dum1)
          }
        }
		#result <- list(insertion=insertion, deletion=deletion,inversion=inversion, duplication=duplication, tandem=tandem_repeats, translocation=translocation, other=other)
		#return(dum1)
    }, mc.cores=node)
	#browser()
	good_num <- sum(sapply(temp_result,function(X) sum(sapply(X, length)!=0)))
    print(paste("Found", good_num, "SV events so far."))
    gc()
  }
  close(fs)
  #browser()
  saveRDS(temp_result, file="temp_result.rds")
  result <- do.call("rbind", lapply(temp_result, function(X) do.call("rbind", X)))
  a <- split(result, f=result$SV)
	
  result_list <- lapply(a, function(X) split(X, f=X$qseqid))
  print(sapply(result_list, length))
  
  #insertion <- split(a[["insertion"]], f=a[["insertion"]]$qseqid)
  #deletion <- split(a[["deletion"]], f=a[["deletion"]]$qseqid)
  #inversion <- split(a[["inversion"]], f=a[["inversion"]]$qseqid)
  #duplication  <- split(a[["duplication"]], f=a[["duplication"]]$qseqid)
  #tandem_repeats <- split(a[["tandem_repeats"]], f=a[["tandem_repeats"]]$qseqid)
  #translocation <- split(a[["translocation"]], f=a[["translocation"]]$qseqid)
  #other <- split(a[["other"]], f=a[["other"]]$qseqid)

  #result_list <- list(insertion=insertion, deletion=deletion, duplication=duplication, tandem=tandem_repeats, translocation=translocation, inversion=inversion, other=other)
  #print("summary of findings:")
  #print(paste("Insertion:", length(insertion)))
  #print(paste("Deletion:", length(deletion)))
  #print(paste("Duplication:", length(duplication)))
  #print(paste("Tandem repeats:", length(tandem_repeats)))
  #print(paste("Inversion:", length(inversion)))
  #print(paste("translocation:", length(translocation)))
  #print(paste("Other:", length(other)))

  return(result_list)
}



generate_image_translocation <- function(aligned_list, output=T) {
  #library(rBLAST)
  if(length(aligned_list)==0) {
    print("Nothing to plot")
    return(NULL)
  }
  for(j in 1:length(aligned_list)) {
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
  }
}


generate_image <- function(aligned_list, output=T) {
  #library(rBLAST)
  if(length(aligned_list)==0) {
    print("Nothing to plot")
    return(NULL)
  }
  #browser()
  for(j in 1:length(aligned_list)) {
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
  }
}

convert_alignment <- function(alignment, scale.down=1) {

    if(length(unique(alignment$sseqid))>1) {
      print("More than one chromosome.")
      next
    }
    min_Y <- min(c(alignment$qstart, alignment$qend))
    max_Y <- max(c(alignment$qstart, alignment$qend))
    min_X <- min(c(alignment$sstart, alignment$send))
    max_X <- max(c(alignment$sstart, alignment$send))
    temp <- array(0, dim=c(max_X-min_X+1,
                           max_Y-min_Y+1))
    for(i in 1:nrow(alignment)) {
      #print(i)
      length <- min(abs(alignment$send[i] - alignment$sstart[i]), abs(alignment$qend[i] - alignment$qstart[i]))
      y_start <- alignment$qstart[i]
      y_end <- y_start+length
      #print(paste(i, start, end))
      x_range <- alignment$sstart[i]:alignment$send[i]
      y_range <- y_start:y_end
      for(k in 1:length(x_range)) {
        if(alignment$sstrand[i]=="plus"){
          temp[x_range[k]-min_X, y_range[k]-min_Y] <- 1
        } else {
          temp[x_range[k]-min_X, y_range[k]-min_Y] <- -1
        }
      }
    }
  return(temp)
}


call_SV_pacbio <- function(BAM, fastq, candidates=NULL, ref=NULL, max_dist = 1000, window_masker_db=NULL, block_size=1000, graph=T, node=8, out="SV_result") {
  library(rtracklayer)
  if(is.null(candidates)) {
    candidates <- get_candidate(BAM=BAM, softclip_length = 50)
    writeLines(candidates, "candidates.txt")
  } else {
    candidates=readLines(candidates)
  }


  start_time <- Sys.time()
  sv <- find_sv_hs_blastn_parallel(FASTQ = fastq,  blast_ref=ref, masker_db=window_masker_db, candidates=candidates,max_dist=max_dist, block_size=block_size, node=node)
  end_time <- Sys.time()
  print(paste(start_time, end_time))
  print(paste("Total run time:", end_time - start_time))
  
  system(paste0("mkdir ", out))
  setwd(out)
  saveRDS(sv, file="sv_all.rds")
  
  for(i in names(sv)) {
	  write.csv(do.call("rbind", sv[[i]]), file=paste0(i, ".csv"))
  }

  result <- SV_report(sv)
  write.csv(result, "result_summary.csv")
  SV_collapsed <- Collapse_SV(result)
  export(SV_collapsed, "SV.gtf")
  
  if(graph) {
    system("mkdir insertion")
	system("mkdir deletion")
	system("mkdir inversion")
	system("mkdir duplication")
	system("mkdir CSUB")
	system("mkdir translocation")
	system("mkdir other")
 
    print("Generating graphs...")
    setwd("insertion")
    generate_image(sv$insertion, output=T)
    setwd("../deletion")
    generate_image(sv$deletion, output=T)
    setwd("../inversion")
    generate_image(sv[["inversion"]], output=T)
    setwd("../duplication")
    generate_image(sv[["duplication"]], output=T)
    setwd("../CSUB")
    generate_image(sv[["CSUB"]], output=T)
	  setwd("../translocation")
    generate_image_translocation(sv[["translocation"]], output=T)
    setwd("../other")
    generate_image(sv[["other"]], output=T)
    setwd("../../")
  }

  return(sv)
}

SV_report <- function(SV) {
  result <- list()
  count <- 0
	temp <- SV[["deletion"]]
	if(!is.null(temp)) {	
		for(j in 1:length(temp)) {
		  
		  dum <- temp[[j]]
		  if(nrow(dum)!=2) {
		    next
		  }
		  count <- count+1
			if(unique(dum$sstrand)=="pos") {
			  start <- dum$send[1]
        end <- dum$sstart[2]
		  } else {
		    end <- dum$send[1]
		    start <- dum$sstart[2]
		  }
		  seqname <- unique(dum$sseqid)
		  read <- unique(dum$qseqid)
		  result[[count]] <- c(seqname, start, end, read, "deletion")
		}
	}
	
	temp <- SV[["insertion"]]
	if(!is.null(temp)) {	
	  for(j in 1:length(temp)) {
	    
	    dum <- temp[[j]]
	    if(nrow(dum)!=2) {
	      next
	    }
	    count <- count+1
	    start <- mean(as.numeric(dum$send[1]), as.numeric(dum$sstart[2]))
      end <- start
	    seqname <- unique(dum$sseqid)
	    read <- unique(dum$qseqid)
	    result[[count]] <- c(seqname, start, end, read, "insertion")
	  }
	}
	
	temp <- SV[["duplication"]]
	if(!is.null(temp)) {	
	  for(j in 1:length(temp)) {

	    dum <- temp[[j]]
	    if(nrow(dum)!=2) {
	      next
	    }
	    count<- count+1
	    if(unique(dum$sstrand)=="pos") {
	      start <- dum$send[1]
	      end <- dum$sstart[2]
	    } else {
	      end <- dum$send[1]
	      start <- dum$sstart[2]
	    }
	    seqname <- unique(dum$sseqid)
	    read <- unique(dum$qseqid)
	    result[[count]] <- c(seqname, start, end, read, "duplication")
	  }
	}
	
	temp <- SV[["inversion"]]
	if(!is.null(temp)) {	
	  for(j in 1:length(temp)) {
	    dum <- temp[[j]]
	    if(nrow(dum)>3) {
	      next
	    } 
	    count <- count+1
	    if(nrow(dum)==2) {
        ind <- which(dum$cov==min(dum$cov))[1]
	    } else if(nrow(dum)==3) {
	      ind <- 2
	    }
	    foo <- table(dum$sstrand)
      if(names(sort(foo)[2])=="neg") {
        start <- dum$send[ind]
        end <- dum$sstart[ind]
      } else {
        end <- dum$send[ind]
        start <- dum$sstart[ind]
      }
	    seqname <- unique(dum$sseqid)
	    read <- unique(dum$qseqid)
	    result[[count]] <- c(seqname, start, end, read, "inversion")
	   }
	}
	
	  temp <- SV[["translocation"]]
	  if(!is.null(temp)) {	
	    for(j in 1:length(temp)) {
	      dum <- temp[[j]]
	      dum_list <- split(dum, f=dum$sseqid)
	      for(i in length(dum_list)) {
	        dum_temp <- dum_list[[i]]
	        if(nrow(dum_temp)==1) {
	          dum_list[[i]] <- dum_temp
	        } else {
	          dum_temp[,c("sstart", "send")] <- t(apply(dum_temp[,c("sstart", "send")], 1, function(X) sort(as.numeric(X), decreasing = F)))
	          #x <- GRanges(seqnames=unique(dum_temp$sseqid), ranges = IRanges(start=dum_temp$sstart, end=dum_temp$send))
	          cov_1 <- as.numeric(dum_temp$length)/(max(as.numeric(dum_temp$qend))-min(as.numeric(dum_temp$qstart)))
	          if(max(cov_1)>= 0.99) {
	            ind <- which(cov_1 == max(cov_1))
	            dum_list[[i]] <- dum_temp[ind,]
	          } else if (nrow(dum_temp)>2) {
	            break
	          } else {
	            dum_list[[i]] <- dum_temp[1,]
	          }
	        }
	      }
        dum <- do.call("rbind", dum_list)
        dum <- dum[order(dum$qstart),]
        if(nrow(dum)!=2) {
          next
        } else {
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
          result[[count]] <- c(seqname, start, end, read, "translocation")
        }
	    }
	  }
	  #browser()
	  result <- do.call("rbind", result)

	  return(result)
}

Collapse_SV <- function(result_table) {
  #result_table <- result[result[,5]!="translocation",]
  result_table[,2:3] <- t(apply(result_table[,2:3], 1, function(X) sort(as.numeric(X))))
  #result_table <- rbind(result_table, result[result[,5]=="translocation",])
  colnames(result_table) <- c("seqname", "start", "end", "read_id", "sv_type")
  result_table <- as.data.frame(result_table)
  result_table <- result_table[order(result_table$seqname, result_table$start),]
  result_range <- GRanges(result_table)
  combined <- reduce(result_range, min.gapwidth=25)
  foo <- findOverlaps(result_range, combined)
  combined$support_reads <- table(subjectHits(foo))
  combined$sv_type <- result_table$sv_type[queryHits(foo)[match(sort(unique(subjectHits(foo))), as.integer(subjectHits(foo)))]]
  return(combined)
}


##################################
#6/21/23
#Using hs-blastn to improve on speed
##################################
find_sv_hs_blastn_parallel <- function(FASTQ, candidates= NULL, min_sv_size = 50, max_dist=10000, blast_ref="Z:/genome/SacCer3/SacCer3_blast", masker_db=NULL, block_size=1000, exclude_chrM=T, graph=T, node=8) {
  #install.packages('rBLAST', repos = 'https://mhahsler.r-universe.dev')
  library(ShortRead)
  library(parallel)
  #browser()
  Sys.setenv(PATH=paste("/opt/hs-blastn/0.0.5/bin/", Sys.getenv("PATH"), sep=":"))
  
  fs <- FastqStreamer(FASTQ, n=block_size)
  aligned_result <- list()
  count <- 0
  streamer_count <-  0
  temp_result <- list()
  print("Start detection process")
  while(length(fq <- yield(fs))) {
    streamer_count <- streamer_count+1
    print(paste("Imported", streamer_count*block_size, "reads."))
    ind <- which(is.element(fq@id, candidates)==T)
	if(length(ind)==0) {
		next
	}
    count <- count+length(ind)
    print(paste("Processing", count, "candidates..."))
    a <- DNAStringSet(fq@sread[ind])
    names(a) <- fq@id[ind]
	#browser()
	width <- width(a)
	names(width) <- names(a)
	writeXStringSet(a, file="temp.fa")
    #aligned <- predict(ref, a, BLAST_args = paste("-task megablast -window_masker_db /net/nfs-irwrsrchnas01/labs/xwu/genome/hg38/hg38.fa.counts.obinary -num_threads", node), custom_format = "qseqid qlen qstart qend sseqid sstart send sstrand length pident bitscore")
    if(is.null(masker_db)) {
		system(paste("hs-blastn align -query temp.fa -db", blast_ref, "-outfmt 6 -evalue 1e-10 -max_target_seqs 20 -num_threads", node, "> temp.out"))
	} else {
		system(paste("hs-blastn align -query temp.fa -db", blast_ref, "-window_masker_db", masker_db, "-max_target_seqs 20 -outfmt 6 -evalue 1e-10 -num_threads", node, "> temp.out"))
	}
	if(file.size("temp.out")) {
		aligned <- read.delim("temp.out", header=F)
	} else {
		next
	}
	aligned$qlen <- width[match(aligned[,1], names(width))]
	colnames(aligned) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen")

	if(exclude_chrM) {
      aligned <- aligned[aligned$bitscore>=200&aligned$pident>=95&aligned$sseqid!="chrM",]
    }
    aligned[,c("qstart", "qend")] <- t(apply(aligned[,c("qstart", "qend")], 1, function(X) sort(as.numeric(X), decreasing = F)))
    aligned <- aligned[order(aligned$qseqid, aligned$qstart),]
    aligned_list <- split(aligned, f=aligned$qseqid)
    print("Finished alignment with Blast.")
    #browser()
    temp_result[[streamer_count]] <- mclapply(1:length(aligned_list), function(i) {
		if(i %% 1000 ==0) {
			print(i)
		}

		dum <- aligned_list[[i]]
		dum$sstrand <- apply(dum[,c("sstart", "send")], 1, function(X) ifelse(X[1]< X[2], "pos", "neg"))
    dum$cov <- (dum$qend - dum$qstart +1)/dum$qlen
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
					#print(paste("translocation",ncol(dum)))
					return(dum)
				} else {
					return(NULL)
				}
			}
    }
    
    #add in segmenation step to resolve multiple reports involving distant regions
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
          #normal
          return(NULL)
        } else if  (nrow(dum1)<=3 & length(unique(dum1$sstrand))==2) {
          temp <- convert_alignment(dum1)
          a1 <- rowSums(abs(temp)) #reference
          b1 <- colSums(abs(temp)) #query
          
          if (sum(b1==0)<=min_sv_size & sum(a1==0)<=min_sv_size & sum(a1>1)<=min_sv_size & sum(b1>1)<=min_sv_size) {
            dum1$SV <- "inversion"
		        return(dum1)
          }
        } else if(nrow(dum1)<=3 & length(unique(dum1$sstrand))==1) {
          temp <- convert_alignment(dum1)
          a1 <- rowSums(abs(temp)) #reference
          b1 <- colSums(abs(temp)) #query

          if(sum(a1==0)==0 & sum(b1==0)==0) {
            return(NULL)
          }

          if (sum(a1==0)<=min_sv_size & sum(a1>1)>=min_sv_size &sum(b1>1)<=min_sv_size) {
            dum1$SV <- "duplication"
			      return(dum1)
          } else if (sum(a1==0)>=min_sv_size & sum(a1>1) <= min_sv_size & sum(b1==0)>=min_sv_size&sum(b1>1)<=min_sv_size) {
            dum1$SV <- "CSUB"
          } else if (sum(a1==0)<=min_sv_size & sum(a1>1) <= min_sv_size & sum(b1==0)>=min_sv_size&sum(b1>1)<=min_sv_size) {
            dum1$SV <- "insertion"
			      return(dum1)
          } else if(sum(b1==0)<=min_sv_size & sum(b1>1) <=min_sv_size & sum(a1==0)>=min_sv_size & sum(a1>1)<=min_sv_size) {
            dum1$SV <- "deletion"
			      return(dum1)
          } else {
            dum1$SV <- "other"
			      return(dum1)
          }

        } else if(nrow(dum1) >3 & length(unique(dum1$sstrand))==1) {
          temp <- convert_alignment(dum1)
          #reference
          a1 <- rowSums(abs(temp))
          #query
          b1 <- colSums(abs(temp))
          if (sum(b1>1)<=min_sv_size & sum(a1>2)>=min_sv_size & sum(a1==0)<=min_sv_size) {
            dum1$SV <- "duplication"
			      #print(dum1)
			      return(dum1)
          } else if (sum(a1==0)<=min_sv_size & sum(b1>=2)>=min_sv_size) {
            #repeats in reference
            return(NULL)
          } else {
            dum1$SV <- "other"
			      return(dum1)
          }
        }
		#result <- list(insertion=insertion, deletion=deletion,inversion=inversion, duplication=duplication, tandem=tandem_repeats, translocation=translocation, other=other)
		#return(dum1)
    }, mc.cores=node)

	good_num <- sum(sapply(temp_result,function(X) sum(sapply(X, length)!=0)))
    print(paste("Found", good_num, "SV events so far."))
    gc()
  }
  close(fs)
  #browser()
  good_result <- temp_result[!sapply(temp_result, is.null)]
  saveRDS(good_result, file="temp_result.rds")
  result <- do.call("rbind", lapply(good_result, function(X) do.call("rbind", X)))
  a <- split(result, f=result$SV)
  result_list <- lapply(a, function(X) split(X, f=X$qseqid))
  print(sapply(result_list, length))
return(result_list)
}
