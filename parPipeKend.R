#!/usr/bin/env Rscript
library(methods)
args = commandArgs(trailingOnly=TRUE)

samp_dir <- args[1]
out_dir <- args[2]
library(VennDiagram)
library(pracma)
library(biomaRt)
library(ggplot2)
library(tidyr)
library(dplyr)
library(GO.db)
library(biomaRt)
library(VennDiagram)
library(ggplot2)
library(doParallel)
library(hypeR)
library(purrr)
library(tidyverse)
library(bedr)
library(beepr)


makeConcatHIT <- function(sampleNames, samp_dir) {
  totAFE <- read.delim(paste(samp_dir, sampleNames[1], "/hit/", sampleNames[1],".AFEPSI", sep = ""))
  totAFE$sampleNum <- 1
  for (i in 2:length(sampleNames)) {
    tempAFE <- read.delim(paste(samp_dir, sampleNames[i], "/hit/", sampleNames[i],".AFEPSI", sep = ""))
    tempAFE$sampleNum <- i
    totAFE <- rbind(totAFE, tempAFE)
  }
  totAFE <- totAFE %>% dplyr::select(gene, sampleNum, exon, strand, AFEPSI) %>% arrange(gene, exon)
  totAFE$bedr.loc <- unlist(lapply(totAFE$exon, function(x) paste("chr", x,sep="")))
  totAFE <- totAFE %>% separate(bedr.loc, c("chr", "loc"), ":", remove = FALSE) %>% dplyr::select(gene, sampleNum, strand, bedr.loc, chr, loc, AFEPSI)
  # totAFE <- totAFE[(totAFE$AFEPSI > 0 & totAFE$AFEPSI < 1.0),]
  return(totAFE = totAFE)
}

makeConcatrMATS <- function(sampleNames, samp_dir) {
  totSE <- read.delim(paste(samp_dir, sampleNames[1], "/rmats/SE.MATS.JC.txt", sep = ""))
  totSE$sampleNum <- 1
  for (i in 2:length(sampleNames)) {
    tempSE <- read.delim(paste(samp_dir, sampleNames[i], "/rmats/SE.MATS.JC.txt", sep = ""))
    tempSE$sampleNum <- i
    totSE <- rbind(totSE, tempSE)
  }
  totSE <- totSE[(!is.na(totSE$IncLevel1)),] %>% dplyr::select(GeneID, strand, sampleNum, exonStart_0base,
                                                               exonEnd, IncLevel1)


  #totSE$bedr.loc <- unlist(lapply(totSE$exon, function(x) paste("chr", x,sep="")))
  return(totSE)
}

mergeProm.par <- function(concatSamp) {
  finAFE <- data.frame(
    gene = c("place"),
    strand = c("place"),
    bedr.loc = c("place")
  )
  cl <- makePSOCKcluster(detectCores() - 4)
  registerDoParallel(cl)
  finAFE <- rbind(finAFE, foreach (i = 1:length(unique(concatSamp$gene)), .combine = rbind) %dopar% {
    library(foreach)
    library(doParallel)
    itAFE <- concatSamp[concatSamp$gene %in% unique(concatSamp$gene)[i],]

    foreach (j = 1:length(unique(itAFE$strand)), .combine = rbind) %dopar% {
      library(bedr)
      library(dplyr)
      s.itAFE <- itAFE[itAFE$strand %in% unique(itAFE$strand)[j],]
      merged.s.itAFE <- bedr.merge.region(bedr.sort.region(s.itAFE$bedr.loc))

      data.frame(
        gene = rep(unique(concatSamp$gene)[i], length(merged.s.itAFE)),
        strand = rep(unique(itAFE$strand)[j], length(merged.s.itAFE)),
        bedr.loc = merged.s.itAFE
      )
    }
  })
  finAFE <- finAFE[!(finAFE$gene %in% "place"),]
  finAFE <- finAFE %>% separate(bedr.loc, c("chr", "loc"), ":")
  return(finAFE = finAFE)
}


mergeSE.par <- function(concatSamp) {
  cl <- makePSOCKcluster(detectCores() - 4)
  registerDoParallel(cl)
  finSE <- data.frame(
    gene = c("place"),
    strand = c("place"),
    exon_start = c("place"),
    exon_end = c("place"),
    exon = c("place")
  )
  finSE <- rbind(finSE, foreach (i = 1:length(unique(concatSamp$GeneID)), .combine = rbind) %dopar% {
    library(foreach)
    library(doParallel)
    itSE <- concatSamp[concatSamp$GeneID %in% unique(concatSamp$GeneID)[i],]
    foreach (j = 1:length(unique(itSE$strand)), .combine = rbind) %dopar% {
      library(bedr)
      library(dplyr)
      s.itSE <- itSE[itSE$strand %in% unique(itSE$strand)[j],]

      s.o <- paste(unlist(lapply(s.itSE$exonStart_0base, function(x) paste("chr1:", x,sep=""))), "-", s.itSE$exonEnd, sep="")
      s.t <- bedr::bedr.merge.region(bedr::bedr.sort.region(s.o))
      merged.s.itSE <- strsplit(s.t, ":")
      merged.s.itSE <- sapply(merged.s.itSE, "[[", 2)
      merged.s.itSE <- t(data.frame(strsplit(merged.s.itSE, "-")))
      rownames(merged.s.itSE) <- NULL
      data.frame(
        gene = rep(unique(concatSamp$GeneID)[i], length(merged.s.itSE)/2),
        strand = rep(unique(itSE$strand)[j], length(merged.s.itSE)/2),
        exon_start = merged.s.itSE[,1],
        exon_end = merged.s.itSE[,2],
        exon = paste(merged.s.itSE[,1], "-", merged.s.itSE[,2], sep="")
      )
    }
  }
  )
  finSE <- finSE[2:length(finSE$gene),]
  return(finSE = finSE)

}

twoInt <- function(x1, x2, y) {
  y.s <- as.integer(strsplit(y, "-")[[1]])
  inter <- length(intersect(seq(x1, x2), seq(y.s[1], y.s[2]))) > 0
  return(inter)
}


makeDFperSE.par <- function(cat.SE, cat.AFE, merge.SE, merge.AFE) {


  cl <- makePSOCKcluster(detectCores() - 4)
  registerDoParallel(cl)
  indSElistT <- foreach (i = 1:length(unique(merge.SE$gene))) %dopar% {
    indSElist <- list()

    twoInt <- function(x1, x2, y) {
      y.s <- as.integer(strsplit(y, "-")[[1]])
      inter <- length(intersect(seq(x1, x2), seq(y.s[1], y.s[2]))) > 0
      return(inter)
    }
    library(dplyr)
    g <- unique(merge.SE$gene)[i]
    g.AFE <- merge.AFE[merge.AFE$gene %in% g,]
    g.SE <- merge.SE[merge.SE$gene %in% g,]


    for (j in 1:length(unique(g.SE$strand))) {

      s <- unique(g.SE$strand)[j]
      gs.SE <- g.SE[g.SE$strand %in% s,]
      gs.AFE <- g.AFE[g.AFE$strand %in% s,]

      for (k in 1:length(unique(gs.SE$exon))) {

        e <- unique(gs.SE$exon)[k]
        gse.SE <- gs.SE[gs.SE$exon %in% e,]



        tempRmats <- cat.SE[(cat.SE$GeneID %in% g &
                               cat.SE$strand %in% s),]

        locs <- as.list(as.data.frame(t(tempRmats[,4:5])))

        inEx <- unlist(lapply(1:length(locs), function(x) twoInt(locs[[x]][1], locs[[x]][2], e)))

        indSElist[[(length(indSElist)+1)]] <- data.frame(
          tempRmats[inEx, c(1:6)]
        )

        promCols <- length(gs.AFE$loc)
        if (promCols == 0) {
          indSElist[[(length(indSElist))]] <- cbind(indSElist[[(length(indSElist))]], t(rep(NA,2)))
          colnames(indSElist[[(length(indSElist))]]) <- c("gene", "strand", "sampleNum", "seStart", "seEnd", "psiSE", "noProms", "noProms")

        } else {
          indSElist[[(length(indSElist))]] <- cbind(indSElist[[(length(indSElist))]], t(rep(NA,promCols)))
          colnames(indSElist[[(length(indSElist))]]) <- c("gene", "strand", "sampleNum", "seStart", "seEnd", "psiSE", gs.AFE$loc)

        }

        tempHIT <- cat.AFE[cat.AFE$gene %in% g & cat.AFE$strand %in% s & cat.AFE$sampleNum %in% unique(indSElist[[(length(indSElist))]]$sampleNum),]

        for (l in 1:length(unique(tempHIT$sampleNum))) {
          furHIT <- tempHIT[tempHIT$sampleNum %in% unique(tempHIT$sampleNum)[l],]

          indPSI <- unlist(lapply(furHIT$loc, function(y)
            which(unlist(lapply(gs.AFE$loc, function(x)

              twoInt(as.integer(unlist(strsplit(y, "-")))[1], as.integer(unlist(strsplit(y, "-")))[2], x))))))


          indSElist[[length(indSElist)]][indSElist[[length(indSElist)]]$sampleNum %in% unique(tempHIT$sampleNum)[l], 6+indPSI] <-
            rep(furHIT$AFEPSI,
                each =    dim(indSElist[[length(indSElist)]][indSElist[[length(indSElist)]]$sampleNum %in% unique(tempHIT$sampleNum)[l],])[1])

          indSElist[[length(indSElist)]][is.na(indSElist[[length(indSElist)]])] <- 0
        }
      }

    }
    indSElist
  }
  return(unlist(indSElistT,recursive=FALSE))
}

calcRho_noO.par <- function(SEdf) {
  rhoDF <- data.frame(matrix(NA, nrow = 1, ncol = 6))
  cl <- makePSOCKcluster(detectCores() - 4)
  registerDoParallel(cl)
  rhoDF <- foreach (t = 1:length(SEdf), .combine = rbind) %dopar% {
    library(foreach)
    library(doParallel)
    if (dim(SEdf[[t]])[2] > 6 & sum("noProms" %in% colnames(SEdf[[t]])) == 0) {
      # rhoTM <- data.frame(matrix(NA, nrow = 1, ncol = 6))
      rhod <- foreach (i = 7:dim(SEdf[[t]])[2], .combine = rbind) %dopar% {

        proPSI <- SEdf[[t]][,i]
        sePSI <- SEdf[[t]][,6]
        proInd <- proPSI != 1
        seInd <- sePSI != 1
        combInd <- proInd & seInd
        pro = proPSI[combInd]
        se = sePSI[combInd]
        # pro <- SEdf[[t]][,i]
        # se <- SEdf[[t]][,6]
        if (length(se) > 0 & length(pro) > 0) {
          rho <- suppressWarnings(cor(se, pro, method = "kendall"))
        } else {
          rho <- NA
        }
        # rhoTM <- rbind(rhoTM, c(unique(SEdf[[t]][,1]), unique(SEdf[[t]][,2]), min(SEdf[[t]][,4]), max(SEdf[[t]][,5]), colnames(SEdf[[t]])[i], rho))
        # rhoTM <- rhoTM[!is.na(rhoTM[,1]),]
        rhoTM <- c(unique(SEdf[[t]][,1]), unique(SEdf[[t]][,2]), min(SEdf[[t]][,4]), max(SEdf[[t]][,5]), colnames(SEdf[[t]])[i], rho)
        rhoTM
      }
      rhod
    }

  }
  return(data.frame(rhoDF))
}

calcRho_noZO <- function(SEdf) {
  rhoDF <- data.frame(matrix(NA, nrow = 1, ncol = 6))
  for (t in 1:length(SEdf)) {
    if (dim(SEdf[[t]])[2] > 6 & sum(grepl("noProms", colnames(SEdf[[t]]))) == 0) {
      for (i in 7:dim(SEdf[[t]])[2]) {

        proPSI <- SEdf[[t]][,i]
        sePSI <- SEdf[[t]][,6]
        proInd <- proPSI != 1 & proPSI != 0
        seInd <- sePSI != 1 & sePSI != 0
        combInd <- proInd & seInd
        pro = proPSI[combInd]
        se = sePSI[combInd]
        # pro <- SEdf[[t]][,i]
        # se <- SEdf[[t]][,6]
        if (length(se) > 0 & length(pro) > 0) {
          rho <- suppressWarnings(cor(se, pro, method = "spearman"))
        } else {
          rho <- NA
        }
        rhoDF <- rbind(rhoDF, c(unique(SEdf[[t]][,1]), unique(SEdf[[t]][,2]), min(SEdf[[t]][,4]), max(SEdf[[t]][,5]), colnames(SEdf[[t]])[i], rho))
      }

    }

  }
  return(rhoDF)
}

calcRho_noZ <- function(SEdf) {
  rhoDF <- data.frame(matrix(NA, nrow = 1, ncol = 6))
  for (t in 1:length(SEdf)) {
    if (dim(SEdf[[t]])[2] > 6 & sum(grepl("noProms", colnames(SEdf[[t]]))) == 0) {
      for (i in 7:dim(SEdf[[t]])[2]) {

        proPSI <- SEdf[[t]][,i]
        sePSI <- SEdf[[t]][,6]
        proInd <- proPSI != 0
        seInd <- sePSI != 0
        combInd <- proInd & seInd
        pro = proPSI[combInd]
        se = sePSI[combInd]
        # pro <- SEdf[[t]][,i]
        # se <- SEdf[[t]][,6]
        if (length(se) > 0 & length(pro) > 0) {
          rho <- suppressWarnings(cor(se, pro, method = "spearman"))
        } else {
          rho <- NA
        }
        rhoDF <- rbind(rhoDF, c(unique(SEdf[[t]][,1]), unique(SEdf[[t]][,2]), min(SEdf[[t]][,4]), max(SEdf[[t]][,5]), colnames(SEdf[[t]])[i], rho))
      }

    }

  }
  return(rhoDF)
}

calcRho <- function(SEdf) {
  rhoDF <- data.frame(matrix(NA, nrow = 1, ncol = 6))
  for (t in 1:length(SEdf)) {
    if (dim(SEdf[[t]])[2] > 6 & sum(grepl("noProms", colnames(SEdf[[t]]))) == 0) {
      for (i in 7:dim(SEdf[[t]])[2]) {
        pro <- SEdf[[t]][,i]
        se <- SEdf[[t]][,6]
        if (length(se) > 0 & length(pro) > 0) {
          rho <- suppressWarnings(cor(se, pro, method = "spearman"))
        } else {
          rho <- NA
        }
        rhoDF <- rbind(rhoDF, c(unique(SEdf[[t]][,1]), unique(SEdf[[t]][,2]), min(SEdf[[t]][,4]), max(SEdf[[t]][,5]), colnames(SEdf[[t]])[i], rho))
      }

    }

  }
  return(rhoDF)
}

se.classifyDF.par <- function(rhoDF) {
  posL <- list()
  disL <- list()
  kbL <- list()
  disAltL <- list()
  disClass <- c(-20011, seq(-20000, 19000, 1000), 20001, 20011)
  seqs <- lapply(seq(1:(length(disClass)-1)), function(x) c(disClass[x], disClass[x+1]))
  prom <- lapply(rhoDF$X5, function(x) as.integer(unlist(strsplit(x, "-"))))
  statDF <- foreach (i = 1:length(rhoDF[,1]), .combine = rbind) %dopar% {
    if (rhoDF$X2[i] == "+") {
      ##upstream
      if (as.integer(prom[[i]][2]) <= rhoDF$X3[i] & prom[[i]][1] <= rhoDF$X3[i]) {
        posL <- "up"
        disL <- prom[[i]][2]-rhoDF$X3[i]
        if (disL < -20000) {
          disAltL <- disClass[1]+2
        } else if (disL > 20000) {
          disAltL <- disClass[length(disClass)]-2
        } else {
          disAltL <- disL
        }
        kbL <- which(unlist(lapply(seqs, function(x) disAltL %in% seq(x[1], x[2]-1))))
      } else {#if (as.integer(prom[[i]][2]) >= rhoDF$X4[i] & prom[[i]][1] >= rhoDF$X4[i]) {
        posL <- "down"
        disL <- prom[[i]][1]-rhoDF$X4[i]
        if (disL < -20000) {
          disAltL <- disClass[1]+2
        } else if (disL > 20000) {
          disAltL <- disClass[length(disClass)]-2
        } else {
          disAltL <- disL
        }
        kbL <- which(unlist(lapply(seqs, function(x) disAltL %in% seq(x[1], x[2]-1))))
      }

    } else {
      ##upstream
      if (as.integer(prom[[i]][1]) >= rhoDF$X4[i] & as.integer(prom[[i]][2]) >= rhoDF$X4[i]) {
        posL <- "up"
        disL <- rhoDF$X4[i]-prom[[i]][1]
        if (disL < -20000) {
          disAltL <- disClass[1]+2
        } else if (disL > 20000) {
          disAltL <- disClass[length(disClass)]-2
        } else {
          disAltL <- disL
        }
        kbL <- which(unlist(lapply(seqs, function(x) disAltL %in% seq(x[1], x[2]-1))))
      } else {#if (as.integer(prom[[i]][1]) <= rhoDF$X3[i] & as.integer(prom[[i]][2]) <= rhoDF$X3[i]) {
        posL <- "down"
        disL <- rhoDF$X3[i]-prom[[i]][2]
        if (disL < -20000) {
          disAltL <- disClass[1]+2
        } else if (disL > 20000) {
          disAltL <- disClass[length(disClass)]-2
        } else {
          disAltL <- disL
        }

        kbL <- which(unlist(lapply(seqs, function(x) disAltL %in% seq(x[1], x[2]-1))))
      }
    }
    c(posL, disL, disAltL, kbL)
  }
  rhoDF$pos <- unlist(statDF[,1])
  #print(1)
  rhoDF$dis <- unlist(statDF[,2])
  #print(2)
  rhoDF$disAlt <- unlist(statDF[,3])
  #print(3)
  rhoDF$kb <- unlist(statDF[,4])
  return(rhoDF)
}

adjustAFEtoProm.par <- function(promDF) {
  cl <- makePSOCKcluster(detectCores() - 4)
  registerDoParallel(cl)
  p.chromStart <- as.integer(promDF$chromStart)
  p.chromEnd <- as.integer(promDF$chromEnd)


  df.t <- foreach (i = 1:length(promDF$gene), .combine = rbind) %dopar% {
    if (promDF$strand[i] == "+") {

      chromStart <- (p.chromStart[i]-300)
      chromEnd <- (p.chromStart[i]+100)

    } else {

      chromStart <- (p.chromEnd[i]-100)
      chromEnd <- (p.chromEnd[i]+300)
    }

    c(chromStart, chromEnd)
  }


  promDF$chromStart <- as.character(df.t[,1])
  promDF$chromEnd <- as.character(df.t[,2])
  return(promDF)
}

quantExt_toBED <- function(rho.df, HITfor08) {
  rho.df <- rho.df[!is.na(rho.df$rho),]
  quantSep <- list()
  qs <- seq(0, 1, .166666666)
  quantSep[[1]] <- rho.df[rho.df$rho >= quantile(rho.df$rho, qs[1], na.rm = TRUE) & rho.df$rho < quantile(rho.df$rho, qs[2], na.rm = TRUE),]
  quantSep[[2]] <- rho.df[rho.df$rho > quantile(rho.df$rho, qs[2], na.rm = TRUE) & rho.df$rho < quantile(rho.df$rho, qs[3], na.rm = TRUE),]
  quantSep[[3]] <- rho.df[rho.df$rho >= quantile(rho.df$rho, qs[3], na.rm = TRUE) & rho.df$rho < quantile(rho.df$rho, qs[4], na.rm = TRUE),]
  quantSep[[4]] <- rho.df[rho.df$rho > quantile(rho.df$rho, qs[4], na.rm = TRUE) & rho.df$rho <= quantile(rho.df$rho, qs[5], na.rm = TRUE),]
  quantSep[[5]] <- rho.df[rho.df$rho >= quantile(rho.df$rho, qs[5], na.rm = TRUE) & rho.df$rho < quantile(rho.df$rho, qs[6], na.rm = TRUE),]
  quantSep[[6]] <- rho.df[rho.df$rho > quantile(rho.df$rho, qs[6], na.rm = TRUE) & rho.df$rho <= quantile(rho.df$rho, qs[7], na.rm = TRUE),]

  chrHIT <- HITfor08[,c(1, 7)]
  chrHIT <- chrHIT[!duplicated(chrHIT), ]
  chrNum <- list()
  for (i in 1:6) {
    chrNum[[i]] <- merge(quantSep[[i]], chrHIT, by="gene")
  }

  boxplotify <- data.frame(rho = c(quantSep[[1]]$rho, quantSep[[2]]$rho, quantSep[[3]]$rho, quantSep[[4]]$rho, quantSep[[5]]$rho, quantSep[[6]]$rho),
                           quantile = c(rep("Sextile 1", length(quantSep[[1]]$rho)), rep("Sextile 2", length(quantSep[[2]]$rho)),
                                        rep("Sextile 3", length(quantSep[[3]]$rho)), rep("Sextile 3", length(quantSep[[4]]$rho)),
                                        rep("Sextile 4", length(quantSep[[5]]$rho)), rep("Sextile 5", length(quantSep[[6]]$rho))))


  (p3 <- ggplot(boxplotify, aes(x=quantile, y=rho, fill =  quantile)) +
      geom_boxplot(fill = c("deeppink4", rep("grey", 3), "deeppink4")) + theme_bw() +
      theme(axis.text.x = element_text(angle = 90)) + xlab("Sextile (Rho)") + ylab("Rho"))
  print(p3)
  saveRDS(p3, paste(out_dir, "sextileRho.RDS", sep = ""))
  listBED <- list()

  for (i in 1:6) {
    listBED[[i]] <- data.frame(chrom = chrNum[[i]]$chr,
                               chromStart = chrNum[[i]]$exon_start,
                               chromEnd = chrNum[[i]]$exon_end,
                               gene = chrNum[[i]]$gene,
                               name = seq(length(chrNum[[i]]$gene)),
                               strand = chrNum[[i]]$strand
    )
    listBED[[i]] <- listBED[[i]][!(duplicated(listBED[[i]][, c(1:3, 6)])),]
  }

  return(listBED)

}

SEquantExt_toBED <- function(rho.df, HITfor08) {
  rho.df <- rho.df[!is.na(rho.df$rho),]
  quantSep <- list()
  qs <- seq(0, 1, .166666666)
  quantSep[[1]] <- rho.df[rho.df$rho >= quantile(rho.df$rho, qs[1], na.rm = TRUE) & rho.df$rho < quantile(rho.df$rho, qs[2], na.rm = TRUE),]
  quantSep[[2]] <- rho.df[rho.df$rho > quantile(rho.df$rho, qs[2], na.rm = TRUE) & rho.df$rho < quantile(rho.df$rho, qs[3], na.rm = TRUE),]
  quantSep[[3]] <- rho.df[rho.df$rho >= quantile(rho.df$rho, qs[3], na.rm = TRUE) & rho.df$rho < quantile(rho.df$rho, qs[4], na.rm = TRUE),]
  quantSep[[4]] <- rho.df[rho.df$rho > quantile(rho.df$rho, qs[4], na.rm = TRUE) & rho.df$rho <= quantile(rho.df$rho, qs[5], na.rm = TRUE),]
  quantSep[[5]] <- rho.df[rho.df$rho >= quantile(rho.df$rho, qs[5], na.rm = TRUE) & rho.df$rho < quantile(rho.df$rho, qs[6], na.rm = TRUE),]
  quantSep[[6]] <- rho.df[rho.df$rho > quantile(rho.df$rho, qs[6], na.rm = TRUE) & rho.df$rho <= quantile(rho.df$rho, qs[7], na.rm = TRUE),]
  # quantSep <- lapply(list(c(0.00, 0.25), c(0.25, 0.50), c(0.50, 0.75), c(0.75, 1.00)), function(x) rho.df[rho.df$rho > quantile(rho.df$rho, x[1], na.rm = TRUE) & rho.df$rho <= quantile(rho.df$rho, x[2], na.rm = TRUE),])






  chrHIT <- HITfor08[,c(1, 7)]
  chrHIT <- chrHIT[!duplicated(chrHIT), ]
  chrNum <- list()
  for (i in 1:6) {
    chrNum[[i]] <- merge(quantSep[[i]], chrHIT, by="gene")
  }

  listBED <- list()

  for (i in 1:6) {
    listBED[[i]] <- data.frame(chrom = chrNum[[i]]$chr,
                               chromStart = chrNum[[i]]$exon_start,
                               chromEnd = chrNum[[i]]$exon_end,
                               gene = chrNum[[i]]$gene,
                               name = seq(length(chrNum[[i]]$gene)),
                               strand = chrNum[[i]]$strand,
                               dis = chrNum[[i]]$dis
    )
    listBED[[i]] <- listBED[[i]][!(duplicated(listBED[[i]][, c(1:3, 6)])),]
  }
  t.BED <- listBED
  TSSdist.df <- data.frame(DistanceToTSSfromSE =
                             c(t.BED[[1]]$dis/1000,
                               t.BED[[2]]$dis/1000,
                               t.BED[[3]]$dis/1000,
                               t.BED[[4]]$dis/1000,
                               t.BED[[5]]$dis/1000,
                               t.BED[[6]]$dis/1000),
                           #unlist(lapply(t.BED, function(x) x$dis)),
                           Sextile =
                             c(rep("S1", length(t.BED[[1]]$dis)),
                               rep("S2", length(t.BED[[2]]$dis)),
                               rep("S3", length(t.BED[[3]]$dis)),
                               rep("S3", length(t.BED[[4]]$dis)),
                               rep("S5", length(t.BED[[5]]$dis)),
                               rep("S6", length(t.BED[[6]]$dis))))

  (p<-ggplot(TSSdist.df, aes(x=Sextile, y=DistanceToTSSfromSE/1000, fill=Sextile)) +
      geom_boxplot(fill = c("deeppink4", rep("grey", 3), "deeppink4")) + ylab("KB to AFE TSS from SE") + theme_bw())

  print(p)
  saveRDS(p, paste(out_dir, "TSSdis.RDS", sep = ""))
  return(listBED)

}

concatHITfor08 <- function(sampleNames, samp_dir) {
  totAFE <- read.delim(paste(samp_dir, sampleNames[1], "/hit/", sampleNames[1],".AFEPSI", sep = ""))
  totAFE$sampleNum <- 1
  for (i in 2:length(sampleNames)) {
    tempAFE <- read.delim(paste(samp_dir, sampleNames[i], "/hit/", sampleNames[i],".AFEPSI", sep = ""))
    tempAFE$sampleNum <- i
    totAFE <- rbind(totAFE, tempAFE)
  }
  totAFE <- totAFE %>% dplyr::select(gene, sampleNum, exon, strand, AFEPSI) %>% arrange(gene, exon)
  totAFE$bedr.loc <- unlist(lapply(totAFE$exon, function(x) paste("chr", x,sep="")))
  totAFE <- totAFE %>% separate(bedr.loc, c("chr", "loc"), ":", remove = FALSE)

  return(totAFE = totAFE)
}





SE_AFE_Corr <- function(samp_dir, rho_REM = c("none", "0", "1", "0,1")[3], out_dir) {
  samples <- list.files(samp_dir)[grep("SRR",list.files(samp_dir),fixed=TRUE)]

  concatHIT <- makeConcatHIT(samples, samp_dir = samp_dir)
  saveRDS(concatHIT, paste(out_dir, "concatAFE.RDS", sep = ""))

  concatRmats <- makeConcatrMATS(samples, samp_dir = samp_dir)
  saveRDS(concatRmats, paste(out_dir, "concatSE.RDS", sep = ""))

  print("Concat Done")
  finSE <- mergeSE.par(concatRmats)
  saveRDS(finSE, paste(out_dir, "mergeSE.RDS", sep = ""))

  finAFE <- mergeProm.par(concatHIT)
  saveRDS(finAFE, paste(out_dir, "mergeAFE.RDS", sep = ""))
  print("merges Done")
  dfPerSE <- makeDFperSE.par(concatRmats, concatHIT, finSE, finAFE)
  saveRDS(dfPerSE, paste(out_dir, "dfPerSE.RDS", sep = ""))

  print("DF-ify Done")
  if (rho_REM == "none") {
    r <- calcRho(dfPerSE)
  } else if (rho_REM == "0") {
    r <- calcRho_noZ(dfPerSE)
  } else if (rho_REM == "1") {
    r <- calcRho_noO.par(dfPerSE)
  } else if (rho_REM == "0,1") {
    r <- calcRho_noZO(dfPerSE)
  }
  print("Rho Done")
  saveRDS(r, paste(out_dir, "rho.RDS", sep = ""))
  # r <- readRDS("/projectnb/evolution/zwakefield/GTEx_dl/upsampleOut/rho.RDS")
  newRHO <- r
  newRHO <- newRHO[2:length(newRHO$X2),]
  newRHO$X3 <- as.integer(newRHO$X3)
  newRHO$X4 <- as.integer(newRHO$X4)
  newRHO$X6 <- as.double(newRHO$X6)

  intL <- list()
  for (i in 1:length(newRHO$X1)) {
    intL[[i]] <- twoInt(newRHO$X3[i], newRHO$X4[i], newRHO$X5[i])
  }
  nnRHO <- newRHO[!unlist(intL),]

  classifyRho <- se.classifyDF.par(nnRHO)
  saveRDS(classifyRho, paste(out_dir, "classifyRho.RDS", sep = ""))


  classifyRho$kb <- as.integer(classifyRho$kb)
  classifyRho$kb <- as.integer(classifyRho$kb)
  classifyRho$dis <- as.integer(classifyRho$dis)
  classifyRho$disAlt <- as.double(classifyRho$disAlt)
  classifyRho <- classifyRho %>% arrange(as.integer(kb))
  classifyRho$groups <- factor(classifyRho$kb, levels=order(unique(classifyRho$kb)))
  classifyRho$posgroups <- factor(classifyRho$pos, levels=sort(unique(classifyRho$pos), decreasing = TRUE))
  classifyRho$color <- c(rep("deeppink4", as.integer(table(classifyRho$kb < 22)[2])), rep("gray", as.integer(table(classifyRho$kb < 22)[1])))
  # classifyRho <- readRDS("/Users/zacharywakefield/Desktop/samplesDown/ResultsDown/classifyRho.RDS")
  bpL <- c(
    "(-X, -20)", "(-20, -19)", "(-19, -18)","(-18, -17)","(-17, -16)","(-16, -15)","(-15, -14)",
    "(-14, -13)","(-13, -12)","(-12, -11)","(-11, -10)",
    "(-10, -9)", "(-9, -8)","(-8, -7)","(-7, -6)","(-6, -5)","(-5, -4)",
    "(-4, -3)","(-3, -2)","(-2, -1)","(-1, 0)","(0, 1)","(1, 2)","(2, 3)",
    "(3, 4)","(4, 5)","(5, 6)","(6, 7)","(7, 8)",
    "(8, 9)","(9, 10)",

    "(10, 11)","(11, 12)","(12, 13)",
    "(13, 14)","(14, 15)","(15, 16)","(16, 17)","(17, 18)",
    "(18, 19)","(19, 20)", "(20, +X)")

  p.20kb <- ggplot(classifyRho, aes(x=groups, y=X6, fill = color)) +
    geom_boxplot(alpha = c(seq(.2, 1, length = 21), seq(1, .2, length = 21)), fill = c(rep("deeppink4", 21), rep("gray", 21)))+
    theme(axis.text.x = element_text(angle = 90)) +                      # Add median line to boxplot
    stat_summary(fun = mean,
                 geom = "line",
                 aes(group = 1),
                 col = "black") + theme_bw() +
    theme(axis.text.x = element_text(angle = 90)) + xlab("AFE Kb relative to SE") +
    scale_x_discrete(labels=bpL) + ylab("Rho")
  print(p.20kb)

  saveRDS(classifyRho, paste(out_dir, "classifyRho.RDS", sep = ""))
  saveRDS(p.20kb, paste(out_dir, "p20kb.RDS", sep = ""))

  classifyRho$kb <- as.integer(classifyRho$kb)
  prom2Rho <- classifyRho[,c(1,2,5,6)]
  prom2Rho <- prom2Rho %>% separate(X5, into = c("start", "end"))
  colnames(prom2Rho) <- c("gene", "strand", "chromStart", "chromEnd", "rho")
  adjPromRho <- adjustAFEtoProm.par(prom2Rho)
  adjPromRho$chromStart <- as.character(adjPromRho$chromStart)
  adjPromRho$chromEnd <- as.character(adjPromRho$chromEnd)
  colnames(adjPromRho) <- c("gene", "strand", "exon_start", "exon_end", "rho")
  saveRDS(adjPromRho, paste(out_dir, "adjPromRho.RDS", sep =""))
  # adjPromRho <- readRDS("/Users/zacharywakefield/Desktop/samplesDown/ResultsDown/adjPromRho.RDS")
  HITfor08 <- concatHITfor08(samples, samp_dir = samp_dir)
  t.BED <- quantExt_toBED(adjPromRho, HITfor08 = HITfor08)

  for (i in 1:6) {
    t.BED[[i]] <- t.BED[[i]][(abs(as.integer(t.BED[[i]]$chromEnd)-as.integer(t.BED[[i]]$chromStart)) > 2),]
  }
  write.table(t.BED[[1]], paste(out_dir, "afe.s1.bed", sep=""), row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
  write.table(t.BED[[2]], paste(out_dir, "afe.s2.bed", sep=""), row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
  write.table(t.BED[[3]], paste(out_dir, "afe.s3.bed", sep=""), row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
  write.table(t.BED[[4]], paste(out_dir, "afe.s4.bed" , sep=""), row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
  write.table(t.BED[[5]], paste(out_dir, "afe.s5.bed", sep=""), row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
  write.table(t.BED[[6]], paste(out_dir,"afe.s6.bed" , sep=""), row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")


  rhoPreQ <- classifyRho[,c(1:4,6,9)]
  colnames(rhoPreQ) <- c("gene", "strand", "exon_start", "exon_end", "rho", "dis")
  t.se.BED <- SEquantExt_toBED(rhoPreQ, HITfor08 = HITfor08)

  for (i in 1:6) {
    t.se.BED[[i]] <- t.se.BED[[i]][(abs(as.integer(t.se.BED[[i]]$chromEnd)-as.integer(t.se.BED[[i]]$chromStart)) > 7),]
  }
  write.table(t.se.BED[[1]], paste(out_dir, "se.s1.bed", sep=""), row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
  write.table(t.se.BED[[2]], paste(out_dir, "se.s2.bed", sep=""), row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
  write.table(t.se.BED[[3]], paste(out_dir, "se.s3.bed", sep=""), row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
  write.table(t.se.BED[[4]], paste(out_dir, "se.s4.bed" , sep=""), row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
  write.table(t.se.BED[[5]], paste(out_dir, "se.s5.bed", sep=""), row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
  write.table(t.se.BED[[6]], paste(out_dir,"se.s6.bed" , sep=""), row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

  return(list(concatAFE = concatHIT,
              concatSE = concatRmats,
              mergeAFE = finAFE,
              mergeSE = finSE,
              dfPerSE = dfPerSE,
              rho = r,
              classifyRho = classifyRho,
              p.20kb = p.20kb,
              adjAFE = adjPromRho,
              AFEquant = t.BED,
              SEquant = t.se.BED
  ))
}




fullRun <- SE_AFE_Corr(samp_dir = args[1],
                       rho_REM = "1",
                       out_dir = args[2])
saveRDS(fullRun, paste(args[2], "PipeOut.RDS", sep =""))
print("Corr Pipe Done!")
