merge_baypass <- function(env, maxiter, loci) {

    loci <- read.table(loci,h=F,col.names=c("chrom","pos"))
    df <- data.frame(n = 1:nrow(loci))
    xtx <- list()
    betai <- list()
    J <- 1


for (i in c(1:as.integer(maxiter))) {
      
    betai[[J]] <- read.table(paste0("baypass_core_",env,"_",i,"_summary_betai_reg.out"),h=T)
    xtx[[J]] <- read.table(paste0("baypass_core_",env,"_",i,"_summary_pi_xtx.out"),h=T)
      
      betai.df <- betai[[J]]
      xtx.df <- xtx[[J]]
      #head(xtx.df, n = 1L)
      df <- cbind(df, betai.df[c("BF.dB.","M_Pearson")])        
      # names BF
      colnames(df)[ncol(df) - 1] <- paste0("BF",i)
      # names Beta
      colnames(df)[ncol(df)] <- paste0("Pearson",i)
      
      df <- cbind(df, xtx.df["M_XtX"])
      colnames(df)[ncol(df)] <- paste0("M_XtX",i)
      
      J <- J + 1 
}

df %>% rowwise() %>% 
    mutate(
        "{env}_BF" := median(c_across(starts_with("BF"))), 
        "{env}_Pearson" := median(c_across(starts_with("Pearson"))),
        "{env}_XtX" := median(c_across(starts_with("M_XtX")))) %>% 
    select(starts_with(env)) %>%
    write.table(x=.,file=paste0(env,"_med_baypass.txt"),row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
}