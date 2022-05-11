library(dplyr)
library(stringr)
library(tidyr)

meta2env <- function (meta, indv, grps, grp_order){
    
    # environmental data from channel 
    indvs <- read.table(indv, header = F, col.names = "id")
    env <- sub('.indv', '',basename(indv))
    
    # environmental data
    meta <- read.table(meta, header = T , sep = ',') %>% select(id='bioinformatics_id',env)
    
    # Baypass groupings
    grps <- read.table(grps, header = FALSE, col.names = c("id","pop"))
    grp_order <- read.table(grp_order, header = FALSE, sep="",col.names="order")

    # recode discrete variables to numeric for standardization
    if (env == "H1_coast_noncoast"){
        meta <- mutate(meta, env = as.numeric(str_replace_all(meta[,2], c("noncoast" = "0", "coast" = "1"))))
    } else if (env == "H2_Aufeis_above_below_coast"){
        meta <- mutate(meta, env = as.numeric(str_replace_all(meta[,2], c("above" = "-1", "below" = "0", "coast" = "1"))))
    } else if (env == "H3_WShd_up_down_lower"){
        meta <- mutate(meta, env = as.numeric(str_replace_all(meta[,2], c("up" = "-1", "down" = "0", "lower" = "1"))))
    } else {
        meta <- mutate(meta, env = meta[,2])
    }
    
    # Standardize all records
    meta <- mutate(meta, 
        standard_env = (as.numeric(env) - mean(as.numeric(env),na.rm=T)) / sd(as.numeric(env),na.rm=T))[,c(1,4)]
    
    # Remove problematic records to match indv in vcf
    correct <- left_join(grps, meta, by = "id")

    # Write lfmm format
    write.table(file = paste0(env,".env"), x = correct[,c(3)], row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

    correct %>% group_by(pop) %>% summarise(avg = mean(as.numeric(standard_env))) %>% arrange(factor(pop,levels=grp_order$order)) %>% tidyr::pivot_wider(.,names_from = "pop", values_from = avg) %>%
        write.table(x=.,file=paste0(env,".bayenv"), row.names = FALSE, col.names = FALSE, quote = FALSE, sep = " ")
}
