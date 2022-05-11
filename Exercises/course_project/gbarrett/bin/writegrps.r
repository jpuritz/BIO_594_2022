library(dplyr)

writegrps <- function(meta, indv) {

meta <- read.table(file = meta,sep = ",", header = T) %>% select(id='bioinformatics_id',pop=popmap_adjusted,latitude:RASTERVALU)
vcfindv <- read.table(file = indv, sep = ',', header = F, col.names = 'id')
# write env proceeding 2nd col. 
for (i in 3: ncol(meta)) {
  tmp <- meta[,c(1,2,i)] # subset to 3 col.
  
  title <- toString(names(tmp)[3])
  filename <- paste("./",trimws(title),".grp") %>% gsub("[[:space:]]","",.)
  
  #clean <- tmp[!is.na(tmp[,3]),]
  clean <- filter(tmp, !(tmp[,3] == "Na" | tmp[,3] == "" | tmp[,3] == "GCL"))
  #clean <- tmp[!(is.na(as.numeric(tmp[,3]))),]

  l=unique(c(as.character(clean$pop)))
  
  grp <- data.frame(id=clean$id,grp=as.numeric(factor(clean$pop, levels=l)))
  correct <- merge(grp, vcfindv,by = "id",all.x = FALSE, all.y = FALSE)
  #correct <- grp %>% filter(ind %in% vcfindv$)
  
  write.table(file = filename,x = correct, quote = F,sep = "\t",col.names = F,row.names = F)
  
  i <- i + 1 
  }
}