#' Function to process results from GxEScan
#' 
#' Function to process results from GxEScan to perform the two-step steps
#' and produce summaries and plots in an Excel spreadsheet
#' 
#' @param basefilename
#' Base filename of output files from GxEScan or GxEMerge
#' @param alpha
#' Type I error rate
#' @param topnum
#' Number of top hits to report
#' @param includeRanks
#' Include ranks
#' @export
GxEResults <- function(homedir, basefilename,
                       alpha = 0.05, exhaustiveTests,
                       topnum = 10, includeRanks = FALSE,
                       bins, screens, twoStepTests) {
  # Values need for each test
  # Exhaustive tests
  extensions1 <- c("CC_DG", "CC_GxE", "CC_GE", "Case_GE", "Cntl_GE", "CC_2df", "CC_3df", "CC_DGGE")
  df <- c(1, 1, 1, 1, 1, 2, 3, 2)
  ptitle1 <- c("Marginal G", "GxE", "G|E", "Case Only", "Control Only", "2df (G,GxE)", "3df (G,GxE,G|E)", "2df (D|G, G|E)")
  rankName1 <- c("MarginalG", "GxE", "G|E", "CaseOnly", "CntlOnly", "TwoDF", "ThreeDF", "DGGE")
  exhaustiveDescriptions <- c("Marginal G:Standard GWAS analysis of marginal G effect for all M SNPs",
                              "GxE: Standard case-control analysis of GxE interaction for all M SNPs",
                              "G|E: Gene-environment association test",
                              "Case Only: Standard case-only analysis of GxE interaction for all M SNPs",
                              "Control Only: Analysis of E-G association in controls only",
                              "2df: Two degree of freedom joint test of G,GxE for all M SNPs",
                              "3df: Three degree of freedom joint test of G,GxE and E-G association for all M SNPs",
                              "DGGE: Two degree of freedom joint test of D|G, and G|E")
  # 2 step tests
  extensions2 <- c("EG2DF", "EGGxE", "DGGxE", "EDGE")
  twoStepSteps <- data.frame(c(3, 3, 1, 8), c(6, 2, 2, 2))
  ptitle2 <- c("EG|2df", "EG|GxE", "DG|GxE", "EDGE")
  rankName2 <- c("EG|2df", "EG|GxE", "DG|GxE", "EDGE")
  testType <- c("tTest", "2df ChiSq")
  testType1 <- c(1, 1, 1, 2)
  testType2 <- c(2, 1, 1, 1)
  twoStepDescriptions1 <- c(
    "EG|2df: 2-step with screening based on E-G association; subset testing of G, GxE in Step 2",
    "EG|GxE: 2-step with screening based on E-G association; subset testing of GxE in Step 2",
    "DG|GxE: 2-step with screening based on D-G association; subset testing of GxE in Step 2",
    "EDGxE: 2-step with screening based on joint test of E-G and D-G association; subset testing of GxE in Step 2"
  )
  twoStepDescriptions2 <- c(
    "EG|2df: 2-step with screening based on E-G association; weighted testing of G, GxE in Step 2",
    "EG|GxE: 2-step with screening based on E-G association; weighted testing of GxE in Step 2",
    "DG|GxE: 2-step with screening based on D-G association; weighted testing of GxE in Step 2",
    "EDGxE: 2-step with screening based on joint test of E-G and D-G association; weighted testing of GxE in Step 2"
  )
  
  rankOrder <- c(rankName1[exhaustiveTests], rankName2[twoStepTests])
  # Switch to directory with results - Probably not a good idea
  olddir <- getwd();
  setwd(homedir)
  
#  exSummary <- InitializeExhaustiveSummary(alpha)
#  twoStepSummary <- InitializeTwoStepSummary(screens, bins)
# Read in the SNP information  
  info <- data.table::fread(paste(basefilename, ".snpinfo", sep=""), header = TRUE, sep = '\t', stringsAsFactors = FALSE)
  M.info <- nrow(info)
  setkey(info, SNPID)
# Determine if all values will be displayed in QQ and Manhattan plots
  # Not sure this is a good idea - Not every SNP will have converged
  if (M.info < 1e5)
    PlotAll = 1
  else
    PlotAll = 0
 
  numSNPs1 <- rep(NA, 8)
  reqSig1 <- rep(NA, 8)
  numSNPs2 <- rep(NA, 4)
  reqSig2 <- rep(NA, 4)
  testRanks <- data.table(info[,SNPID])
  setnames(testRanks, 1, "SNPID")
  setkey(testRanks, SNPID)
  top <- list()
  top2 <- list()
  gxeData <- list()
  # Two step tests - loads required exhaustive test as needed
  x <- list()
  for (i in 1:4) {
    if (twoStepTests[i]) {
      for (j in 1:2) {
        k <- twoStepSteps[i, j]
        if (is.na(numSNPs1[k])) {
          x1 <- ProcessGxEOut(snpinfo = info, basefilename = basefilename, extension = extensions1[k],
                              df = df[k], alpha = alpha, topnum = topnum, includeRanks = includeRanks,
                              ptitle = ptitle1[k], rankName = rankName1[k], generatePlot = exhaustiveTests[k])
          numSNPs1[k] <- x1$M
          if (exhaustiveTests[k]) {
            reqSig1[k] <- x1$sig
            top[[k]] <- x1$top
            if (includeRanks) {
              setkey(x1$ranks, SNPID)
              testRanks <- x1$ranks[testRanks]
            }
          }
          if (j == 1)
            x[[j]] <- x1$gxeout[,c("SNPID", "STAT", "P")]
          else
            x[[j]] <- x1$gxeout
          rm(x1)
        }
      }
      y <- Process2Step(test1 = x[[1]], test2 = x[[2]],
                        basefilename = paste(basefilename, extensions2[i], sep = "_"),
                        alpha = alpha, screen = screens[i], bins = bins[i],
                        topnum = topnum, includeRanks = includeRanks, rankName = rankName2[i],
                        ptitle = ptitle2[i])
      numSNPs2[i] <- y$m
      reqSig2[i] <- y$sig
      top2[[i]] = y$top
      if (includeRanks) {
        setkey(y$rank, SNPID)
        testRanks <- y$rank[testRanks]
      }
    }
  }
  if (sum(twoStepTests) > 0) {
    rm(x)
    rm(y)
  }
  # Exhaustive tests - skips tests that were loaded for two step tests
  for (i in 1:8) {
    if (exhaustiveTests[i] & is.na(numSNPs1[i])) {
      x <- ProcessGxEOut(snpinfo = info, basefilename = basefilename, extension = extensions1[i],
                         df = df[i], alpha = alpha, topnum = topnum, includeRanks = includeRanks,
                         ptitle = ptitle1[i], rankName = rankName1[i])
      if (exhaustiveTests[i]) {
        numSNPs1[i] <- x$M
        reqSig1[i] <- x$sig
        top[[i]] <- x$top
        if (includeRanks) {
          setkey(x$ranks, SNPID)
          testRanks <- x$ranks[testRanks]
        }
      }
      rm(x)
    }
  }
  # Delete existing Excel output file if it exists
  id <- grep(paste(basefilename, ".xlsx", sep = ""), dir(homedir))
  todelete <- dir(homedir, full.names = TRUE)[id]
  unlink(todelete)
  # Create Excel workbook
  numMethods = sum(exhaustiveTests) + sum(twoStepTests)
  wb <- loadWorkbook(paste(basefilename, ".xlsx", sep = ""), create = TRUE)
  WriteSummaryWorksheet(workbook = wb, exhaustiveTests = exhaustiveTests, exhaustiveNames = ptitle1,
                        numSNPs1 = numSNPs1, sig1 = alpha, twoStepTests = twoStepTests, twoStepNames = ptitle2,
                        numSNPs2 = numSNPs2, sig2 = screens, bins = bins, step1Test = twoStepSteps[,1])
  for (i in 1:8) {
    if (exhaustiveTests[i])
      WriteExhaustiveSheet(workbook = wb, testName = rankName1[i], basefilename = basefilename, extension = extensions1[i],
                           description = exhaustiveDescriptions[i],
                           M = numSNPs1[i], sig = reqSig1[i], topHits = top[[i]], includeRanks = includeRanks,
                           numMethods = numMethods, testRanks = testRanks, rankOrder = rankOrder)
  }
  for (i in 1:4) {
    Write2StepSheet(workbook = wb, testName = rankName2[i], testType1 = testType[testType1[i]],
                    testType2 = testType[testType2[i]], basefilename = basefilename, extension = extensions2[i],
                    description = twoStepDescriptions1[i], m = numSNPs2[i], sig = reqSig2[i], topHits = top2[[i]],
                    includeRanks = includeRanks, testRanks = testRanks, rankOrder = rankOrder)
    WriteWeightedSheet(workbook = wb, basefilename = basefilename, extension = extensions2[i],
                       description = twoStepDescriptions2[i],
                       topHits = top2[[i]], testType1 = testType[testType1[i]], testType2 = testType[testType2[i]])
  }
  #  Write2StepSheet(wb, "EG2df_Subset", paste(basefilename, "EG2DF", sep = "_"),
#                  "EG|2df: 2-step with screening based on E-G association; subset testing of G, GxE in Step 2",
#                 top2[[1]], numSNPs2[i], reqSig2[i], includeRanks, testRanks, "EG2df", numMethods)
  saveWorkbook(wb)
  setwd(olddir);
  return (0)
}