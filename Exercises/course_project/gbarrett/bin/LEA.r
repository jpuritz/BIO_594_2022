
lea <- function(env, ped,env_file) {
    if (is.character(env))
        if (env != "max.dist" | env != "min.dist"){
            latentFactors = 5
        } else {
            latentFactors = 3
        }
        env <- toString(env) %>% gsub("[[:space:]]","",.)
        geno <- ped2geno(ped, output.file=paste0(env,".geno"), force=TRUE)
        lfmm <- ped2lfmm(ped, output.file=paste0(env,".lfmm"), force=TRUE)
        #geno <- vcf2geno(vcf, output.file=paste0(env,".geno"), force=TRUE)
        #lfmm <- vcf2lfmm(vcf, output.file=paste0(env,".lfmm"), force=TRUE)

        ################
        # Differentiation
        proj.snmf <- snmf(geno,K=latentFactors,entropy=T,ploidy=2,project="new",alpha=10,tolerance=0.0001,repetitions=2,iterations=100,CPU=25,percentage=.75)
        # fst values
        best <- which.min(cross.entropy(proj.snmf, K = latentFactors))
        fst.values <- fst(proj.snmf, K = latentFactors, run = best)
        # z-scores
        n <- dim(Q(proj.snmf, K = latentFactors, run = best))[1]
        fst.values[fst.values<0] <- 0.000001
        GD_z_scores <- sqrt(fst.values*(n - latentFactors)/(1 - fst.values))
        GD_z_scores <- as.data.frame(GD_z_scores)
        colnames(GD_z_scores) <- paste0(env,"_GD_zscores")
        #write.table(x=z.scores,file=paste0(env,"_GD_zscores.txt"),quote=F,row.names=F,col.names=paste0(env,"_z"))
        ################
        ################
        # Association
        proj.lfmm <- lfmm(lfmm, env_file, K = latentFactors, repetitions = 2, project = "new", iterations = 100, burnin = 50, CPU = 25, missing.data = TRUE, random.init = TRUE)
        # z-scores from all repititions
        zv <- data.frame(z.scores(proj.lfmm, K = latentFactors))
        zv %>% rowwise() %>% mutate("{env}_EA_zscores" := median(c_across(everything()))) %>% select(paste0(env,"_EA_zscores")) %>% 
        cbind(.,GD_z_scores) %>% 
        write.table(x=., file = paste0(env,"_zscores.txt"),quote=F,row.names=F,col.names=T,sep="\t")
        #################
}