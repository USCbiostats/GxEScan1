ReadGxEOut <- function(basefilename, extension) {
  res <- data.table::fread(paste(basefilename, "_", extension, ".gxeout", sep = ""), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  res.sub <- with(res, res[!is.na(P)])
  return (res.sub)
}

QQ <- function(p, t, max, PlotAll) {
  data.table::setnames(p, "p")
  data.table::setkeyv(p, cols="p")
  num.snp <- nrow(p)
  
  if (PlotAll == TRUE) {
    qqplot(-log10(seq(1, num.snp)/(num.snp + 1)), -log10(p[,p]), xlim=c(0, max), ylim = c(0, max), main = t, ylab = 'Observed', xlab = 'Expected')
    abline(0, 1)
  } else {
    exp <- -log10(seq(1, num.snp)/(num.snp + 1))
    grpNum <- floor(log(num.snp))
    pSub <- list()
    expSub <- list()
    pSub[[1]] <- p[1:5000]
    expSub[[1]] <- exp[1:5000]
    x <- 5000
    y <- 15000
    i = 2
    while (i < grpNum & y <= num.snp) {
      pSub[[1]] <- p[seq((x + 1), y, 2^(i - 1))]
      expSub[[1]] <- exp[seq((x + 1), y, 2^(i - 1))]
      x <- y
      y <- y + (1e4) * 2^(i - 1)
      i = i + 1
    }
    pSub[[i]] <- p[seq((x + 1), num.snp, 2^(i - 1))]
    expSub[[i]] <- exp[seq((x + 1), num.snp, 2^(i - 1))]
    p <- c(NA)
    exp <- c(NA)
    for (j in 1:i) {
      p <- c(p, pSub[[j]][,p])
      exp <- c(exp, expSub[[j]])
    }
    p <- p[-1]
    exp <- exp[-1]
    num.snp <- length(p)
    qqplot(exp, -log10(p), xlim = c(0, max), ylim = c(0, max), main = t, ylab = 'Observed', xlab = 'Expected')
    abline(0, 1)
  }
}

manhattan <- function(data, sig, min.p, scale, PlotAll, step2) {
  if (PlotAll == TRUE) {
    cutoff = 0
  } else {
    cutoff = 1
  }
  
  data.table::setnames(data, c(1:3), c('p', 'Chr', 'MapInfo'))
  data$y <- -log10(data$p)
  Chr <- list()
  for (i in 1:22) {
    t <- data[data$Chr == i,]
    Chr[[i]] <- t[t$y > cutoff,]
  }
  
  x <- list()
  for (i in 1:22) {
    x[[i]] <- (Chr[[i]]$MapInfo - min(Chr[[i]]$MapInfo)) / ((max(Chr[[i]]$MapInfo) + 0.0001) - min(Chr[[i]]$MapInfo))*scale + i - 1
  }
  
  colset <- c(rep(c("blue","olivedrab4", "red", "darkorange3", "purple", "orange"), 3), "blue", "olivedrab4", "darkorange3")
  plot(x[[1]], Chr[[1]]$y, col="blue", xlab = "Chromosome", ylab = "-log10(p)", xlim = c(0, 22), ylim = c(cutoff, min.p), axes = FALSE, pch = 20)
  for (i in 2:22) {
    points(x[[i]], Chr[[i]]$y, col = colset[i], pch = 20)
  }
  
  if (step2 == TRUE) {
    abline(-1*log10(sig), 0, lwd = 1)
  }
  
  axis(1, at = seq(scale/2 - 1, 21 + scale/2 , 2), labels = seq(0, 22, 2), cex.axis = 0.6)
  
  axis(2, at = c(cutoff:min.p), labels = c(cutoff:min.p), cex.axis = 0.65)
}

saveimage <- function(wb, fname, sheetname, picname, row, col) {
  sheet = sheetname
  XLConnect::createName(wb, name = fname, formula = paste(sheet, XLConnect::idx2cref(c(row, col)), sep = "!"), overwrite = TRUE)
  # Note: idx2cref converts indices (row, col) to Excel cell references
  # Put the image created above at the corresponding location
  XLConnect::addImage(wb, filename = picname, name = fname, originalSize = TRUE)
}

wt <- function(pv, alpha, num) {
  rk.pv <- c(1:nrow(pv))
  grp = ceiling(log(rk.pv/num + 1, base = 2))
  with(pv, pv[,Bin := grp])
#  with(pv, pv[,`:=`(Bin = grp)])
  data.table::setkeyv(pv, cols="Bin")
  for (i in 1:max(grp)) {
    with(pv, pv[list(i), Threshold := alpha*2^(-i)/nrow(pv[list(i)])])
  }
#  pv[0, 0] = pv[0, 0]
}

wtplot <- function(results, min.p, title, scale, last.sig, num, PlotAll, fileOutput = TRUE, basefilename) {

  if (fileOutput == TRUE)
    png(paste(basefilename, "Weighted.png", sep = "_"))
  
  if(PlotAll == TRUE){
    cutoff = 0
  } else {
    cutoff = 1
  }
#  setnames(results,c(1:4), c('p', 'grp','wt','MapInfo'))
  with(results, results[,y := -log10(results[,Step2_P])])
#  results[,data.table::`:=`(y = -log10(results[,Step2_P]))]
#  results[,data.table::`:=`(y = -log10(results$Step2_P))]
  
  glist <- list()
  for(i in 1:num){
    t <- results[list(i)]
    with(t, t[,ref := -1 * log10(min(t[,Threshold]))])
    with(t, t[,x := (t[,Rank] - min(t[,Rank])) / ((max(t[,Rank]) + 0.0001) - min(t[,Rank]))*scale + i - 1])
    glist[[i]] <- t
    rm(t)
  }
  ## Scale mapinfo for each Bin to range between 0-1 and add a unit increase for successive Bin

  color<-rep(c("blue","olivedrab4"),100)
  
  with(glist[[1]], plot(glist[[1]][,x], glist[[1]][,y], col="blue", xlab="Bin # for step1 p-values", ylab="-log10(step2 p-values)", xlim=c(0,num), ylim=c(0,min.p), axes=F,pch=19))
  with(glist[[1]], lines(glist[[1]][,x], glist[[1]][,ref], col="black",lwd=2))
  
  if(PlotAll == FALSE){
    for(i in 2:num-1){
      if(nrow(glist[[i]]) >= 5115){
        with(glist[[i]], points(glist[[i]][y > cutoff][,x], glist[[i]][y > cutoff][,y], col=color[i],pch=19))
        with(glist[[i]], rect(-.1 + i - 1, -0.07, 0.8 + i - 1, min(glist[[i]][y > cutoff][,y]) + .1, col = color[i], border = color[i]))
      } else {
        with(glist[[i]], points(glist[[i]][,x], glist[[i]][,y], col = color[i], pch = 19))
      }
      with(glist[[i]], lines(glist[[i]][,x], glist[[i]][,ref], col = "black",lwd = 2))
    }
    
    with(glist[[i]], points(glist[[num]][y > cutoff][,x], glist[[num]][y > cutoff][,y], col = color[num], pch = 19))
    with(glist[[i]], rect(-.1 + num - 1, -0.07, 0.8 + num - 1, min(glist[[num]][y > cutoff][,y]) + .1, col = color[num], border = color[num]))
    with(glist[[i]], lines(glist[[num]][,x], last.sig, col = "black", lwd = 2))
  } else {
    for(i in 2:num-1){
      with(glist[[i]], points(glist[[i]][,x], glist[[i]][,y], col = color[i], pch = 19))
      with(glist[[i]], lines(glist[[i]][,x], glist[[i]][,ref], col = "black", lwd = 2))
    }
    
    
    with(glist[[num]], points(glist[[num]][,x], glist[[num]][,y], col = color[num], pch = 19))
    with(glist[[num]], lines(glist[[num]][,x], last.sig,col = "black",lwd = 2))
    
  }	
  
  decix <- (num / 2) - floor(num / 2)
  if(decix > 0){
    axis(1, at = c(-1.5, seq(.5, num - 0.5, 2)), labels = c(0, seq(1, num, 2)), cex.axis = 0.6)
  } else {
    
    axis(1, at = seq(-.5, num - 0.5, 2), labels = seq(0, num, 2), cex.axis = 0.6)
  }
  
  deciy <- (min.p / 2) - floor(min.p / 2)
  
  axis(2, at = c(0:floor(min.p)), labels = c(0:min.p), cex.axis = 0.65)
  
  title (main = title, sub = "Bin Size=5")
  if (fileOutput == TRUE)
    dev.off()
}

savetitle <- function(wb,sheetname,ref,srow,scol, sideBorder = FALSE) {
  
  XLConnect::mergeCells(wb,sheet=sheetname,reference=ref)
  cs <- XLConnect::createCellStyle(wb)
  XLConnect::setFillForegroundColor(cs, color = XLConnect::XLC$"COLOR.GREY_25_PERCENT")
  XLConnect::setFillPattern(cs, fill = XLConnect::XLC$"FILL.SOLID_FOREGROUND")
  XLConnect::setBorder(cs, side = c("bottom"), type = XLConnect::XLC$"BORDER.MEDIUM",color = c(XLConnect::XLC$"COLOR.BLACK"))
  if (sideBorder)
    XLConnect::setBorder(cs, side = c("right"), type = XLConnect::XLC$"BORDER.MEDIUM",color = c(XLConnect::XLC$"COLOR.BLACK"))
  XLConnect::setCellStyle(wb, sheet = sheetname, row = srow, col = scol, cellstyle = cs)
  
}

# Function to change background of cell to gray
# 
# Function to change background of cell to gray
# 
# @param workbook
# Excel workbook to use
# @param sheet
# Sheet in workbook to use
# @param row
# Row number
# @param column
# Column number
GrayBackground <- function(workbook, sheet, row, column) {
  cs <- XLConnect::createCellStyle(workbook)
  XLConnect::setFillForegroundColor(cs, color = XLConnect::XLC$"COLOR.GREY_25_PERCENT")
  XLConnect::setFillPattern(cs, fill = XLConnect::XLC$"FILL.SOLID_FOREGROUND")
  XLConnect::setCellStyle(object = workbook, sheet = sheet, row = row, col = column, cellstyle = cs)
}

#--# Function to process a two step test from GxEScan
#--# 
#--# Function to take the results for two GxEScan tests that were read in
#--# using the ProcessGxE and produce results for a two step test
#--# 
#--# @param test1
#--# The gxeout output from ProcessGxEOut to be used as the screening test
#--# @param test2
#--# The gxeout output from ProcessGxEOut to be in the second step
#--# @param basefilename
#--# Base filename of output files from GxEScan or GxEMerge also used for plots
#--# @param alpha
#--# Type I error rate
#--# @param screen
#--# Significance level required to be included in second step
#--# @param bins
#--# Number of bins to use in weighted test
#--# @param topnum
#--# Number of top hits to report
#--# @param includeRanks
#--# Indicator to save ranks for Excel output
#--# @param rankName
#--# Name of the column containing the rank
#--# @param ptitle
#--# First part of title for plots
#--# @param generatePlot
#--# Indicator to generate qq and Manhattan plots
#--# @param fileOuput
#--# Indicator to save plots to file - Ignored if generatePlot is FALSE
#--# 
Process2Step <- function(test1, test2, basefilename, alpha, screen, bins,
                         topnum, includeRanks, rankName, ptitle,
                         generatePlot = TRUE, fileOutput = TRUE) {
  data.table::setnames(test1, c("STAT", "P"), c("Step1_STAT", "Step1_P"))
  wt(test1, alpha, bins[1])
  data.table::setkeyv(test1, cols="SNPID")
  data.table::setnames(test2, c("NMISS", "STAT", "P"), c("N", "Step2_STAT", "Step2_P"))
  combinedTests <- test2[test1, nomatch = 0]
  data.table::setkeyv(combinedTests, cols="Step1_P")

  M <- nrow(combinedTests)
  if (M < 1e5) {
    PlotAll = TRUE
  } else {
    PlotAll = FALSE
  }
  data.table::setkeyv(combinedTests, cols="Step1_P")
  s2 <- with(combinedTests, combinedTests[Step1_P < screen])
  with(s2, s2[,c("Bin", "Threshold") := NULL])
  with(combinedTests, combinedTests[,c("SNPID") := NULL])
  s2.size <- nrow(s2)
  
  with(combinedTests, combinedTests[,SigTemp := (combinedTests[,Step2_P] < combinedTests[,Threshold]) * 1])
  with(combinedTests, combinedTests[SigTemp == 1, Sig := "***"])
  with(combinedTests, combinedTests[SigTemp == 0, Sig := ""])
  with(combinedTests, combinedTests[,SigTemp := NULL])
  with(combinedTests, combinedTests[,Rank := 1:M])
  
  wt.hits <- with(combinedTests, combinedTests[Sig == "***"])
  wt.top <- combinedTests[1:topnum]
  wt.top <- rbind(wt.hits[!wt.hits$Rank %in% wt.top$Rank,], wt.top)
  
  rm(wt.hits)
  
  step2.min <- min(combinedTests[,c("Step1_P")])
  maxlimqq <- ceiling(1 - log10(step2.min))
  if (generatePlot == TRUE) {
    if (fileOutput == TRUE)
      png(paste(basefilename, "Step1_QQ.png", sep = "_"))
    with(combinedTests, QQ(combinedTests[,list(Step1_P)], paste(ptitle, ": Step 1 screen", sep = ""), maxlimqq, PlotAll))
    if (fileOutput == TRUE)
      dev.off()
    if (fileOutput == TRUE)
      png(paste(basefilename, "Step1_Manhattan.png", sep = "_"))
    with(combinedTests, manhattan(combinedTests[,list(Step1_P, CHR, BP)], 0, maxlimqq, 0.7, PlotAll, 0))
    title(paste(ptitle,": Step 1 screen", sep = ""))
    if (fileOutput == TRUE)
      dev.off()
  }
  
  with(combinedTests, combinedTests[,c("CHR","SNP","BP","A1","N","Step1_STAT","Step2_STAT","Step1_P","Sig") := NULL])
  
  step2.rank <- NULL
  if (s2.size > 0) {        ## If we have SNPs pass to step2
    step2.min<-with(s2,min(s2[,Step2_P]))
    sig_step2_st2 <- alpha / s2.size
    step2.min2 <- min(step2.min, sig_step2_st2)
    maxlimqq <- ceiling(1 - log10(step2.min))
    maxlimman <- ceiling(1 - log10(step2.min2))
    
    if(s2.size < 1e5){
      PlotAll = TRUE
    } else {
      PlotAll = FALSE
    }
    if (generatePlot == TRUE) {
      if (fileOutput == TRUE)
        png(paste(basefilename, "Step2_QQ.png", sep = "_"))
      with(s2, QQ(s2[,list(Step2_P)], paste(ptitle,": Step 2", sep = ""), maxlimqq, PlotAll))
      if (fileOutput == TRUE)
        dev.off()
      if (fileOutput == TRUE)
        png(paste(basefilename, "Step2_Manhattan.png", sep = "_"))
      with(s2, manhattan(s2[,list(Step2_P, CHR, BP)], sig_step2_st2, maxlimman, 0.7, PlotAll, 1))
      title(paste(ptitle,": Step 2", sep = ""))
      if (fileOutput == TRUE)
        dev.off()
    }

    data.table::setkeyv(s2,cols="Step2_P")
    with(s2, s2[,SigTemp := (s2[,Step2_P] < sig_step2_st2)*1])
    with(s2, s2[SigTemp == 1, Sig := "***"])
    with(s2, s2[SigTemp == 0, Sig := ""])
    with(s2, s2[,SigTemp := NULL])
    
    if(includeRanks == TRUE) {
      step2.rank <- with(s2, s2[,list(SNPID)])
      with(step2.rank, step2.rank[,c(rankName) := 1:s2.size])
    }
    s2 <- s2[1:topnum]
  } else {
    sig_step2_st2<-NA
  }

  data.table::setkeyv(combinedTests, cols="Bin")
  
  step2.last <- with(combinedTests, max(combinedTests[,Bin]))
  k <- bins[1] * (2^(step2.last - 1))
  wt.binsig <- alpha * ((1/2)^step2.last)
  wt.sig <- wt.binsig / k
  wt.lnum <- nrow(combinedTests[list(step2.last)])
  wt.ref <- rep(-1 * log10(wt.sig), wt.lnum)
  
  step2.ref <- max(ceiling(1 - log10(with(combinedTests, min(combinedTests[,Step2_P])))), wt.ref[1])
  wtplot(combinedTests, step2.ref + 1, paste(ptitle, "weighted", sep = " "), 0.7, wt.ref, step2.last, PlotAll,
         basefilename = basefilename)
  data.table::setnames(test2, c("N", "Step2_STAT", "Step2_P"), c("NMISS", "STAT", "P"))
  data.table::setnames(test1, c("Step1_STAT", "Step1_P"), c("STAT", "P"))
  
  return (list(M = M, m = s2.size, sig = sig_step2_st2, ranks = step2.rank, top = s2))
}

#--# Function to process a single gxeout file
#--# 
#--# Function to process results from GxEScan for a one step test. The
#--# outputs from two tests can be passed to the Process2Step function
#--# to generate results for two step tests.
#--# 
#--# @param snpinfo
#--# Data table with SNP information
#--# @param basefilename
#--# Base filename of output files from GxEScan or GxEMerge also used for plots
#--# @param extension
#--# Extenstion to base file name for plot
#--# @param df
#--# Degrees of freedom in test
#--# @param alpha
#--# Type I error rate
#--# @param topnum
#--# Number of top hits to report
#--# @param ptitle
#--# First part of title for plots
#--# @param includeRanks
#--# Indicator to save ranks for Excel output
#--# @param rankName
#--# Name of the column containing the rank
#--# @param generatePlot
#--# Indicator to generate qq and Manhattan plots
#--# @param fileOuput
#--# Indicator to save plots to file - Ignored if generatePlot is FALSE
#--# 
#--# @export
ProcessGxEOut <- function(snpinfo, basefilename, extension, df, alpha,
                          topnum, ptitle, includeRanks, rankName,
                          generatePlot = TRUE, fileOutput = TRUE) {
  gxeout <- ReadGxEOut(basefilename = basefilename, extension = extension)
  M.gxe <- nrow(gxeout)
  sig_gxe <- alpha / M.gxe;
  
  if (M.gxe < 1e5)
    PlotAll = 1
  else
    PlotAll = 0
  
  gxe.min <- min(gxeout[,"P"])
  maxlimqq <- ceiling(1 - log10(gxe.min))
  gxe.minman <- min(gxe.min, sig_gxe)
  maxlimman <- ceiling(1 - log10(gxe.minman))
  
  if (includeRanks == TRUE) {
    gxe.rank <- with(gxeout, gxeout[,list(SNPID)])
    with(gxe.rank, gxe.rank[,c(rankName) := 1:M.gxe])
  } else {
    gxe.rank <- NULL;
  }
  
  gxe.top <- gxeout[1:topnum]
  
  with(gxe.top, gxe.top[, SigTemp := (gxe.top[,P] < sig_gxe) * 1])
  with(gxe.top, gxe.top[SigTemp == 1, Sig := "***"])
  with(gxe.top, gxe.top[SigTemp == 0, Sig := ""])
  with(gxe.top, gxe.top[,SigTemp := NULL])
  
  data.table::setkeyv(gxeout,cols="SNPID")
  if (includeRanks == TRUE)
    data.table::setkeyv(gxe.top,cols="SNPID")
  
  gxeout <- snpinfo[gxeout]
  gxe.top <- snpinfo[gxe.top]
  
  if (generatePlot == TRUE) {
    if (fileOutput == TRUE)
      png(paste(basefilename, "_", extension, "_QQ.png", sep = ""))
    with(gxeout, QQ(gxeout[,list(P)], ptitle, maxlimqq, PlotAll))
    if (fileOutput == TRUE)
      dev.off()
    
    if (fileOutput == TRUE)
      png(paste(basefilename, "_", extension, "_Manhattan.png", sep = ""))
    with(gxeout, manhattan(gxeout[,list(P, CHR, BP)], sig_gxe, maxlimman, 0.7, PlotAll, TRUE))
    title(ptitle)
    if (fileOutput == TRUE)
      dev.off()
  }
  
  if (df == 1)
    data.table::setnames(gxe.top, c("NMISS", "BETA", "STAT", "P"), c("N", "Beta", "tTest", "Pvalue"))
  else
    data.table::setnames(gxe.top, c("NMISS", "STAT", "P"), c("N", paste("df",df,"_Chisq", sep = ""), "Pvalue"))
  
  return(list(top = gxe.top, gxeout = gxeout, M = M.gxe, sig = sig_gxe, ranks = gxe.rank))
}

#--# Function to write GxEScan summary sheet
#--# 
#--# Function to write GxEScan summary sheet
#--# listing tests performed, number of SNPs tests converged for,
#--# overall significance level, and required significance
#--#
#--# @param filename
#--# Excel workbook to write summary sheet to
#--# @param exhaustiveTests
#--# List of indicator of exhaustive tests performed
#--# @param exhaustiveNames
#--# Names of exhaustive test names
#--# @param numSNPs1
#--# Number of SNPs that converged for each test
#--# @param sig1
#--# Overall significance level for each test
#--# @param twoStepTests
#--# Indicator of two step tests performed
#--# @param twoStepNames
#--# Names of two step tests
#--# @param numSNPs2
#--# Number of SNPs that made it to step 2
#--# @param sig2
#--# Overall significance level for two step test
#--# @param bins
#--# Number of bins for weighted tests
#--# @param step1Test
#--# Tests performed for step 1 of two step test
#--# @param stacked
#--# Indicator if two step steps are underneath exhaustive test, otherwise side by side
#--# @param decrip2side
#--# Inidicator if method decriptions are written to side of table, otherwise written below
WriteSummaryWorksheet <- function(workbook, exhaustiveTests, exhaustiveNames, numSNPs1, sig1,
                                  twoStepTests, twoStepNames, numSNPs2, sig2, bins, step1Test,
                                  stacked = TRUE, descrip2side = TRUE) {
  x1 <- c("alpha: Overall type I error rate",
          "M:  total number of SNPs tested ( After removing unconverged SNPs )",
          "alpha/M:  Bonferroni corrected significance level for exhaustive testing")
  x2 <- c(x1[1:2],
          "alpha1:  Step 1 significance threshold for 2-step methods",
          "m: Number of SNPs that passed Step1 based on alpha1 level",
          "alpha/m:  Significance level for SNPs tested in Step 2; subset testing",
          "Bin Size:  Initial Bin size for weighted hypothesis testing in Step 2")
  x3 <- c("Marginal G:  Standard GWAS analysis of marginal G effect for all M SNPs",
          "GxE:  Standard case-control analysis of GxE interaction for all M SNPs",
          "G|E: ",
          "Case Only: Standard case-only analysis of GxE interaction for all M SNPs",
          "Control Only: Analysis of E-G association in controls only",
          "2df: Two degree of freedom joint test of G,GxE for all M SNPs (Kraft et al.2007)",
          "3df: Three degree of freedom joint test of G,GxE and E-G association for all M SNPs",
          "DGGE: Two degree of freedom joint test of D|G, and G|E")
  x4 <- c("EG|2df: 2-step with Step-1 screen based on E-G association, Step-2 is 2df joint test of G,GxE",
          "EG|GxE: 2-step with Step-1 screen based on E-G association, Step-2 test of GxE using CC analysis (Murcray et al., 2009)",
          "DG|GxE: 2-step with Step-1 screen based on D-G association, Step-2 test of GxE using CC analysis (Kooperberg and LeBlanc, 2009)",
          "EDGE: 2-step with Step-1 screen based on joint test of E-G and D-G association, Step-2 test of GxE using CC analysis (Gauderman et al., 2013)")
  x5 <- c("   Standard GWAS analysis of marginal G effect for all M SNPs",
          "   Standard case-control analysis of GxE interaction for all M SNPs",
          "   Test of association between G and E",
          "   Standard case-only analysis of GxE interaction for all M SNPs",
          "   Analysis of E-G association in controls only",
          "   Two degree of freedom joint test of G,GxE for all M SNPs (Kraft et al.2007)",
          "   Three degree of freedom joint test of G,GxE and E-G association for all M SNPs",
          "   Two degree of freedom joint test of D|G, and G|E")
  x6 <- c("   2-step with Step-1 screen based on E-G association, Step-2 is 2df joint test of G,GxE",
          "   2-step with Step-1 screen based on E-G association, Step-2 test of GxE using CC analysis (Murcray et al., 2009)",
          "   2-step with Step-1 screen based on D-G association, Step-2 test of GxE using CC analysis (Kooperberg and LeBlanc, 2009)",
          "   2-step with Step-1 screen based on joint test of E-G and D-G association, Step-2 test of GxE using CC analysis (Gauderman et al., 2013)")
  exColNames <- c("Method", "alpha", "M", "alpha/M")
  twoStepColNames <- c("Method", "alpha", "M", "alpha1", "m", "alpha/m", "Bin Size")
  XLConnect::createSheet(object = workbook, name = "Summary")
  XLConnect::setColumnWidth(object = workbook, sheet = "Summary", column = 1, width = 3600)
  
  if (descrip2side)
    stacked = TRUE
  startRow <- 1
  startColumn <- 1
  if (sum(exhaustiveTests) > 0) {
    wrap( wb = workbook, x = "Exhaustive tests", sheetname = "Summary",
          srow = startRow, scol = startColumn, h = FALSE, wrap = FALSE)
    startRow <- startRow + 1

    reqsig <- sig1 / numSNPs1
    summaryDF <- data.frame(exhaustiveNames, rep(sig1, length(exhaustiveTests)), numSNPs1, signif(reqsig, 3))
    colnames(summaryDF) <- exColNames
    if (descrip2side)
      summaryDF$Description <- x5
    XLConnect::writeWorksheet(object = workbook, data = summaryDF[exhaustiveTests,], sheet = "Summary",
                   startRow = startRow, startCol = startColumn, header = TRUE)
    if (descrip2side)
      XLConnect::mergeCells(object = workbook, sheet = "Summary", paste("E", startRow, ":L", startRow, sep = ""))
    if (sum(twoStepTests) == 0 || stacked)
      startRow <- startRow + sum(exhaustiveTests) + 2
    else
      startRow <- startRow + max(sum(exhaustiveTests), sum(twoStepTests)) + 2
      
    XLConnect::writeWorksheet(object = workbook, data = x1, sheet = "Summary",
                   startRow = startRow, startCol = startColumn, header = FALSE)
    if (sum(twoStepTests) == 0 || stacked)
      startRow <- startRow + length(x1) + 1
    else
      startRow <- startRow + max(length(x1), length(x2)) + 1

    if (descrip2side == FALSE) {    
      XLConnect::writeWorksheet(object = workbook, data = x3[exhaustiveTests], sheet = "Summary",
                     startRow = startRow, startCol = startColumn, header = FALSE)
      startRow <- startRow + sum(exhaustiveTests) + 2
    }
    if (stacked == FALSE) {
      startColumn = 9
      startRow = 1
    } else {
      startRow <- startRow + 1
    }
  }

  if (sum(twoStepTests) > 0) {
    wrap(workbook, "Two step tests", "Summary", startRow, startColumn, FALSE, FALSE)
    startRow <- startRow + 1

    reqsig <- sig2 / numSNPs2
    summaryDF <- data.frame(twoStepNames, rep(sig1, length(twoStepTests)), numSNPs1[step1Test],
                            sig2, numSNPs2, signif(reqsig, 3), bins)
    colnames(summaryDF) <- twoStepColNames
    if (descrip2side)
      summaryDF$Description <- x6
    XLConnect::writeWorksheet(object = workbook, data = summaryDF[twoStepTests,], sheet = "Summary",
                   startRow = startRow, startCol = startColumn, header = TRUE)
    if (descrip2side)
      XLConnect::mergeCells(object = workbook, sheet = "Summary", paste("H", startRow, ":T", startRow, sep = ""))
    if (sum(twoStepTests) == 0 || stacked)
      startRow <- startRow + sum(twoStepTests) + 2
    else
      startRow <- startRow + max(sum(exhaustiveTests), sum(twoStepTests)) + 2

    XLConnect::writeWorksheet(object = workbook, data = x2, sheet = "Summary",
                   startRow = startRow, startCol = startColumn, header = FALSE)
    if (stacked == FALSE) {
      if (sum(twoStepTests) == 0 || stacked)
        startRow <- startRow + length(x2) + 1
      else
        startRow <- startRow + max(length(x1), length(x2)) + 1
    }
    if (descrip2side == FALSE) {
      XLConnect::writeWorksheet(object = workbook, data = x4[twoStepTests], sheet = "Summary",
                     startRow = startRow, startCol = startColumn, header = FALSE)
    }
  }
}

wrap <- function(wb,x, sheetname,srow,scol,h,wrap) {
  
  XLConnect::writeWorksheet(wb, x, sheet = sheetname, startRow = srow, startCol = scol,header=h)
  cs <- XLConnect::createCellStyle(wb)
  XLConnect::setWrapText(cs, wrap =wrap)
  XLConnect::setCellStyle(wb, sheet = sheetname, row = srow, col = scol,cellstyle = cs)
  
}

InitializeExhaustiveSummary <- function(alpha) {
  method = c("MarG", "CC", "Case only", "Cntl only", "2df", "3df")
  M <- c(rep(0,6))
  alpha1 <- c(rep(alpha, 6))
  alpha2 <- c(rep(0, 6))
  x <- data.frame(method, M, alpha1, alpha2)
  colnames(x) <-c ("Method","M","alpha","alpha/M")

  return (x)
}

InitializeTwoStepSummary <- function (alphas, bins) {
  method = c("EDGxE", "DG|GxE", "EG|GxE", "EG|2df")
  M <- c(rep(0, 4))
  m <- M
  alpha2 <- M
  x <- data.frame(method, M, alphas, m, alpha2, bins)
  colnames(x) <-c ("Method","M","alpha1","m","alpha/m","Bin Size")
  
  return (x)
}

WriteSummaryDescriptions <- function (wb, startRow, numExhaustive, exhaustiveTests, numTwoStep, twoStepTests) {
  x1 <- c("M:  total number of SNPs tested ( After removing unconverged SNPs )",
          "alpha/M:  Bonferroni corrected significance level for exhaustive testing")
  x2 <- c("alpha1:  Step1 significance threshold for 2-step methods",
          "m: Number of SNPs that passed Step1",
          "alpha1/m:  Significance level for SNPs tested in Step 2; subset testing",
          "Bin Size:  Initial Bin size for weighted hypothesis testing in Step 2")
  x3 <- c("MarG:  Standard GWAS analysis of marginal G effect for all M SNPs",
          "CC:  Standard case-control analysis of GxE interaction for all M SNPs",
          "CO: Standard case-only analysis of GxE interaction for all M SNPs",
          "CntlOnly: Analysis of E-G association in controls only",
          "2df: Two degree of freedom joint test of G,GxE for all M SNPs( Kraft et al.2007)",
          "3df: Three degree of freedom joint test of G,GxE and E-G association for all M SNPs (Gauderman et al., 2013)")
  x4 <- c("EDGxE: 2-step with Step-1 screen based on joint test of E-G and D-G association, Step-2 test of GxE using CC analysis (Gauderman et al., 2013)",
          "DG|GxE: 2-step with Step-1 screen based on D-G association, Step-2 test of GxE using CC analysis (Kooperberg and LeBlanc, 2009)",
          "EG|GxE: 2-step with Step-1 screen based on E-G association, Step-2 test of GxE using CC analysis (Murcray et al., 2009)",
          "EG|2df: 2-step with Step-1 screen based on E-G association, Step-2 is 2df joint test of G,GxE (Gauderman et al., 2013)")
  x <- ""
  if (numExhaustive > 0)
    x <- c(x, "Exhaustive tests", x1, "")
  if (numTwoStep > 0)
    x <- c(x, "2-step tests", x2, "")
  x <- c(x, x3[exhaustiveTests], x4[twoStepTests])
  XLConnect::writeWorksheet(wb, data.frame(x), sheet = "Summary", startRow = startRow, startCol = 1, header = FALSE)
}

#--# Function to write exhaustive test results to Excel spreedsheet
#--# 
#--# Function to write exhaustive test results to Excel spreedsheet
#--# 
#--# @param workbook
#--# Workbook to write results to
#--# @param testName
#--# Name of column in testRanks for the test being output
#--# Also is the name of the sheet created for the results
#--# @param basefilename
#--# Base filename of plots
#--# @param extension
#--# Extension to base filename for plot
#--# @param  description
#--# Description of test
#--# @param M
#--# Number of SNPs that test converged for
#--# @param sig
#--# Required significance level
#--# @param topHits
#--# Data table of top hits for test
#--# @param includeRanks
#--# Indicator to include ranks on sheet
#--# @param numMethods
#--# Number of methods included in analysis
#--# @param testRanks
#--# Ranks for all tests performed
#--# @param rankOrder
#--# Order to display rank columns
WriteExhaustiveSheet <- function(workbook, testName, basefilename, extension, description, M, sig,
                                 topHits, includeRanks, numMethods, testRanks, rankOrder) {
  alphabet<- c("I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z")
  startcol <- NCOL(topHits)
  lastcol <- NCOL(topHits) + NCOL(testRanks) - 2
  firstLetter <- alphabet[startcol - 8]
  lastLetter <- alphabet[lastcol - 8]

  sheet = testName
  XLConnect::createSheet(workbook, name = sheet)
  XLConnect::setColumnWidth(workbook, sheet = sheet, column = 1, width = 3000)
  XLConnect::setColumnWidth(workbook, sheet = sheet, column = 6, width = 2600)
  
  wrap(workbook, description, sheet, 26, 1, FALSE, FALSE)
  wrap(workbook, paste("Number of SNPs tested = ", M, "; significance threshold = ", signif(sig, 3), sep=""),
       sheet, 27, 1, FALSE, FALSE)
  
  if (includeRanks == FALSE) {
    with(topHits, topHits[,SNPID := NULL])
    XLConnect::writeWorksheet(workbook, topHits, sheet = sheet, startRow = 29, startCol = 1, header = TRUE)
  } else {
    for (i in startcol:lastcol)
      XLConnect::setColumnWidth(object = workbook, sheet = sheet, column = i, width = 2800)
    topHits.names <- colnames(topHits[,-1])
    topHits <- testRanks[topHits]
    with(topHits, topHits[,SNPID:=NULL])
#    print(c(topHits.names, testName, rankOrder[-which(rankOrder %in% c(testName))]))
    data.table::setcolorder(topHits, c(topHits.names, testName, rankOrder[-which(rankOrder %in% c(testName))]))
    data.table::setkeyv(topHits, c(testName))
    wrap(workbook,
         "Ranks:  A blank means the SNP did not converge or did not pass the Step-1 screen for the 2-step approach",
         sheet, 28, startcol, FALSE, FALSE)
    
    XLConnect::writeWorksheet(workbook, "Ranks", sheet = sheet, startRow = 29, startCol = startcol, header = FALSE)

    savetitle(workbook, sheet, paste(firstLetter, "29:", lastLetter, "29", sep=""), 29, c(startcol:lastcol))
    XLConnect::writeWorksheet(workbook, topHits, sheet = sheet, startRow = 30, startCol = 1, header = TRUE)
  }
  # written at end to avoid issues with changing column widths
  saveimage(workbook, "graph", sheet, paste(basefilename, "_", extension, "_QQ.png", sep = ""), 1, 1)
  saveimage(workbook, "graph", sheet, paste(basefilename, "_", extension, "_Manhattan.png", sep = ""), 1, 11)
}
#--# Function to write 2 step results to Excel worksheet
#--# 
#--# Function to write 2 step results to Excel worksheet
#--# 
#--# @param workbook
#--# Workbook to write results to
#--# @param testName
#--# Name of column in testRanks for the test being output
#--# @param testType1
#--# Type of statistic for step 1, e.g., ttest, chisq
#--# @param testType2
#--# Type of statistic for step 2, e.g., ttest, chisq
#--# @param basefilename
#--# Base filename of plots
#--# @param extension
#--# Extension to base filename for plot
#--# Also used to create the name of the sheet created for the results
#--# @param  description
#--# Description of test
#--# @param m
#--# Number of SNPs that test passed step 1
#--# @param sig
#--# Required significance level
#--# @param topHits
#--# Data table of top hits for test
#--# @param includeRanks
#--# Indicator to include ranks on sheet
#--# @param testRanks
#--# Ranks for all tests performed
#--# @param rankOrder
#--# Order to display rank columns
Write2StepSheet <- function(workbook, testName, testType1, testType2, basefilename, extension,
                            description, m , sig, topHits, includeRanks, testRanks, rankOrder) {
  alphabet<- c("I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z")
#  firstLetter <- alphabet[NCOL(top) - 7]
#  startcol <- NCOL(top)
#  letter <- alphabet[numMethod + 1]
#  lastcol <- numMethod + NCOL(top)
  
  startcol <- NCOL(topHits)
  lastcol <- NCOL(topHits) + NCOL(testRanks) - 2
  firstLetter <- alphabet[startcol - 8]
  lastLetter <- alphabet[lastcol - 8]
  
  sheet = paste(extension, "Subset", sep = "_")
  XLConnect::createSheet(workbook, name = sheet)

  if (m > 0) {
    XLConnect::setColumnWidth(workbook, sheet = sheet, column = 6, width = 2600)
    XLConnect::setColumnWidth(workbook, sheet = sheet, column = 8, width = 2600)
    XLConnect::setColumnWidth(workbook, sheet = sheet, column = 9, width = 2600)
    XLConnect::setColumnWidth(workbook, sheet = sheet, column = 1, width = 3000)
    
    wrap(workbook, description, sheet, 51, 1, FALSE, FALSE)
    wrap(workbook, paste("Number of SNPs tested in Step 2= ", m,"; Step2 significance threshold = ",
                         signif(sig,3), sep=""), sheet, 52, 1, FALSE, FALSE)

    saveimage(workbook, "graph", sheet, paste(basefilename, extension, "Step2_QQ.png", sep = "_"), 26, 1)
    saveimage(workbook, "graph", sheet, paste(basefilename, extension, "Step2_Manhattan.png", sep = "_"), 26, 11)

    XLConnect::writeWorksheet(workbook, "Step 1", sheet = sheet, startRow = 53, startCol = 6, header = FALSE)
    savetitle(workbook, sheet, paste("F", "53:", "G", "53", sep=""), 53, c(6:7), sideBorder = TRUE)
    XLConnect::writeWorksheet(workbook, "Step 2", sheet = sheet, startRow = 53, startCol = 8, header = FALSE)
    if (startcol == 11) {
      savetitle(workbook, sheet, paste("H", "53:", "I", "53", sep=""), 53, c(8:9))
      clst <- c(1:6, 9, 10, 7, 8, 11)
    } else {
      savetitle(workbook, sheet, paste("H", "53:", "J", "53", sep=""), 53, c(8:10))
      clst <- c(1:6, 10, 11, 7:9, 12)
    }
    topHits <- topHits[, clst, with = FALSE]
    if (includeRanks == FALSE) {
      XLConnect::writeWorksheet(workbook, topHits[,-1], sheet = sheet, startRow = 54, startCol = 1, header = TRUE)
    } else {
      data.table::setkeyv(topHits, cols="SNPID")
      topHits.names <- colnames(topHits[,-1])
      topHits <- testRanks[topHits]
      with(topHits, topHits[,SNPID:=NULL])
      data.table::setcolorder(topHits, c(topHits.names, testName, rankOrder[-which(rankOrder %in% c(testName))]))

      wrap(workbook,
           "Ranks:  A blank means the SNP did not converge or did not pass the Step-1 screen for the 2-step approach",
           sheet, 52, startcol, FALSE, FALSE)
      XLConnect::writeWorksheet(workbook, "Ranks", sheet = sheet, startRow = 53, startCol = startcol, header = FALSE)
      savetitle(workbook, sheet, paste(firstLetter, "53:", lastLetter, "53", sep=""), 53, c(startcol:lastcol))
      
      data.table::setkeyv(topHits, c(testName))
      XLConnect::writeWorksheet(workbook, topHits, sheet = sheet, startRow = 54, startCol = 1, header = TRUE)
    }
    sc <- 6
    XLConnect::writeWorksheet(workbook, testType1, sheet = sheet, startRow = 54, startCol = sc, header = FALSE)
    sc <- sc + 1
    XLConnect::writeWorksheet(workbook, "p-value", sheet = sheet, startRow = 54, startCol = sc, header = FALSE)
    sc <- sc + 1
    if (startcol == 12) {
      XLConnect::writeWorksheet(workbook, "beta", sheet = sheet, startRow = 54, startCol = sc, header = FALSE)
      sc <- sc + 1
    }
    XLConnect::writeWorksheet(workbook, testType2, sheet = sheet, startRow = 54, startCol = sc, header = FALSE)
    sc <- sc + 1
    XLConnect::writeWorksheet(workbook, c("p-value"), sheet = sheet, startRow = 54, startCol = sc, header = FALSE)
    GrayBackground(workbook = workbook, sheet = sheet, row = 54, column = c(6:sc))
  }
  saveimage(workbook, "graph", sheet, paste(basefilename, extension, "Step1_QQ.png", sep = "_"), 1, 1)
  saveimage(workbook, "graph", sheet, paste(basefilename, extension, "Step1_Manhattan.png", sep = "_"), 1, 11)
}
#--# Function to write weighted 2 step results to Excel worksheet
#--# 
#--# Function to write weighted 2 step results to Excel worksheet
#--# 
#--# @param workbook
#--# Workbook to write results to
#--# @param basefilename
#--# Base filename of plots
#--# @param extension
#--# Extension to base filename for plot
#--# Also used to create the name of the sheet created for the results
#--# @param  description
#--# Description of test
#--# @param topHits
#--# Data table of top hits for test
#--# @param testType1
#--# Type of statistic for step 1, e.g., ttest, chisq
#--# @param testType2
#--# Type of statistic for step 2, e.g., ttest, chisq
WriteWeightedSheet <- function(workbook, basefilename, extension, description, topHits,
                               testType1, testType2) {
  sheet = paste(extension, "Weighted", sep = "_")
  XLConnect::createSheet(workbook, name = sheet)

  XLConnect::setColumnWidth(workbook, sheet = sheet, column = 7, width = 2600)
  XLConnect::setColumnWidth(workbook, sheet = sheet, column = 9, width = 2600)
  XLConnect::setColumnWidth(workbook, sheet = sheet, column = 10, width = 2600)
  XLConnect::setColumnWidth(workbook, sheet = sheet, column = 12, width = 2600)
  XLConnect::setColumnWidth(workbook, sheet = sheet, column = 3, width = 3000)

  wrap(workbook, description, sheet, 26, 1, FALSE, FALSE)
  
  XLConnect::writeWorksheet(workbook, "Step 1", sheet = sheet, startRow = 27, startCol = 6, header = FALSE)
  savetitle(workbook, sheet, paste("F", "27:", "G", "27", sep=""), 27, c(6:7), sideBorder = TRUE)
  XLConnect::writeWorksheet(workbook, "Step 2", sheet = sheet, startRow = 27, startCol = 8, header = FALSE)
  startcol <- NCOL(topHits)
  if (startcol == 11) {
    savetitle(workbook, sheet, paste("H", "27:", "I", "27", sep=""), 27, c(8:9))
    clst <- c(1:6, 9, 10, 7, 8, 11)
  } else {
    savetitle(workbook, sheet, paste("H", "27:", "J", "27", sep=""), 27, c(8:10))
    clst <- c(1:6, 10, 11, 7:9, 12)
  }
  topHits <- topHits[, clst, with = FALSE]
  XLConnect::writeWorksheet(workbook, topHits[,-1], sheet = sheet, startRow = 28, startCol = 1, header = TRUE)
  
  sc <- 6
  XLConnect::writeWorksheet(workbook, testType1, sheet = sheet, startRow = 28, startCol = sc, header = FALSE)
  sc <- sc + 1
  XLConnect::writeWorksheet(workbook, "p-value", sheet = sheet, startRow = 28, startCol = sc, header = FALSE)
  sc <- sc + 1
  if (startcol == 12) {
    XLConnect::writeWorksheet(workbook, "beta", sheet = sheet, startRow = 28, startCol = sc, header = FALSE)
    sc <- sc + 1
  }
  XLConnect::writeWorksheet(workbook, testType2, sheet = sheet, startRow = 28, startCol = sc, header = FALSE)
  sc <- sc + 1
  XLConnect::writeWorksheet(workbook, c("p-value"), sheet = sheet, startRow = 28, startCol = sc, header = FALSE)
  GrayBackground(workbook = workbook, sheet = sheet, row = 28, column = c(6:sc))
  
  saveimage(workbook, "graph", sheet, paste(basefilename, extension, "Weighted.png", sep = "_"), 1, 1)
}