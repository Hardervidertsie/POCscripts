# script for plotting bargraphs of PI data
require(data.table)
require(plyr)
require(ggplot2)
require(reshape2)


rm(list=ls())
options(stringsAsFactors = FALSE)
dir()

# there are 20 niii plates and 1 nii plates and 1 new plate for btg2 and p21 . check this 
setwd("E:/POC screen/PI POC data/PI/fourthrun new pipeline/output niii ook met no masked all plates")
my.data <- read.table(file = "Summary Data with the 1016.txt", sep = "\t", header = T)

getwd()

head(my.data)


my.data <- as.data.table(my.data)
my.data[,plateID:=gsub("_CDKN1A", "", plateID)]
my.data[,plateID:=gsub("_BTG2", "", plateID)]
(unique(my.data[,plateID]))
(unique(my.data[,cell_line]))


#my.data[ treatment %in% "Staurosporin" & cell_line %in% "Srxn1",]

#my.data[ treatment %in% "Antimycin A deox" & cell_line %in% "P53",]


#summarize replicates


unique(my.data$plateID)
unique(my.data$cell_line)
# select plates from selected plates based on responses
#selPlateIDs <- c("2013-01-16", "2013-01-23", "2013-01-31", "2013-02-06",  "2013-03-13", "2013-03-20", "2012_12_19", "2014-06-30", "2014-06-27")

#my.data <- my.data[ plateID %in% selPlateIDs, ]
unique(my.data$variable)

vars.keep <- c("imageCountParentObj", 
               "count_pi_obj_masked_primaryId_AreaShape_Area.DIV.Nuclei_AreaShape_Area_larger_0.05")
my.data.t <- my.data[variable %in% vars.keep[1],]
my.data.t$variable <- NULL
setnames(my.data.t, "value", "imageCountParentObj")

my.data.t$value <- my.data[variable %in% vars.keep[2], value]

setnames(my.data.t, "value", vars.keep[2])



# te weinig cellen --> verwijderen
my.data <- my.data.t

my.data<-my.data[ !(imageCountParentObj < 50 & cell_line != "TP53BP1") ,   ]
my.data<-my.data[!(imageCountParentObj > 2000), ]
my.data<-my.data[ !(imageCountParentObj < 10 & cell_line == "TP53BP1") ,   ]


my.data$treatment <- gsub("Iodoacetamide", "IAA",my.data$treatment )

my.data[, treatment:= gsub("DMSO 25%", "DMSO", treatment)]
my.data[, treatment:= gsub("DMSO 50%", "DMSO", treatment)]
my.data[, treatment:= gsub("DMSO 75%", "DMSO", treatment)]
my.data[, treatment:= gsub("DMSO 100%", "DMSO", treatment)]
my.data[, treatment:= gsub("DMSO max", "DMSO", treatment)]

my.data[ treatment == "Menadione" & dose_uM == 5, "dose_uM"] <- 5
my.data[ treatment == "CDDO" & dose_uM == 0.05, "dose_uM"] <- 0.03
my.data[ treatment == "CDDO" & dose_uM == 0.025, "dose_uM"] <- 0.015
my.data[ treatment == "CDDO" & dose_uM == 0.0125, "dose_uM"] <- 0.0075

unique(my.data$treatment)

comps.keep<- c("BFA", "CDDO","Cisplatin","DEM","DMSO","Etoposide","H2O2",
               "Thapsigargin","Tunicamycin","IAA" ,"Staurosporin")
my.data <- my.data[ treatment %in% comps.keep,]

unique(my.data$cell_line)

my.data[ , cell_line:=gsub( "NRF2" , "NFE2L2", cell_line)]
my.data[ , cell_line:=gsub( "Keap1" , "KEAP1", cell_line)]


cells.keep <- c("Srxn1", "HSPA5", "NFE2L2", "DDIT3", "XBP1","P53", "ATF4", "TP53BP1" ,"KEAP1", "CDKN1A", "BTG2")
my.data <- my.data[ cell_line %in% cells.keep,]

my.data.m<-melt(my.data, measure.vars = c("imageCountParentObj", "count_pi_obj_masked_primaryId_AreaShape_Area.DIV.Nuclei_AreaShape_Area_larger_0.05"))

# remove failed plates (based on work presented in "PI viability selection of plates.pptx")
my.data.m[, plate_cell:= paste(cell_line, plateID)]
sort(unique(my.data.m[, plate_cell]))
plateCellKeep <- c("ATF4 2012_12_19", 
                   "BTG2 2015_02_13", "BTG2 2013-10-16",
                   "CDKN1A 2015_02_13", "CDKN1A 2013-10-16",
                   "DDIT3 2013-01-16", "DDIT3 2013-01-23", "DDIT3 2013-01-31", "DDIT3 2013-02-06", "DDIT3 2013-03-13", "DDIT3 2013-03-20",
                   "HSPA5 2012-11-08", "HSPA5 2013-01-09", "HSPA5 2013-01-16", "HSPA5 2013-01-23", "HSPA5 2013-01-31", "HSPA5 2013-02-06",
                   "KEAP1 2012_12_19", "KEAP1 2013-10-16",
                   "NFE2L2 2012-11-08", "NFE2L2 2013-01-09", "NFE2L2 2013-01-16", "NFE2L2 2013-01-23", "NFE2L2 2013-02-06", "NFE2L2 2013-02-27",
                   "NFE2L2 2013-03-13", "NFE2L2 2013-03-20", "NFE2L2 2014-06-27", "NFE2L2 2014-06-30",
                   "P53 2012-11-02", "P53 2012_12_19",
                   "Srxn1 2012-12-06", "Srxn1 2012-12-13", "Srxn1 2013-02-27", "Srxn1 2013-03-13", "Srxn1 2013-03-20",
                   "TP53BP1 2013-02-13", "TP53BP1 2013-02-20", "TP53BP1 2013-02-27", "TP53BP1 2013-03-13", "TP53BP1 2013-03-20", "TP53BP1 2014-06-30",
                   "XBP1 2013-01-09", "XBP1 2013-01-16", "XBP1 2013-01-23", "XBP1 2013-02-06")
my.data.m <- my.data.m[plate_cell %in% plateCellKeep, ]
head(my.data.m)

?t.test
my.data.m[ treatment == "DMSO"]
my.data.m[,  cell.treat.dose:= paste(treatment, cell_line, dose_uM)]


my.data_s[is.na(SD), ]
my.data.m <- my.data.m[ treatment != "Staurosporin"]
all.vars <-  unique(my.data.m$cell.treat.dose)

all.vars <- all.vars[ !grepl("DMSO", all.vars)]

all.vars <- all.vars[ !grepl("ATF4", all.vars)]


stats.out = alist()
stats.out.data = alist()
for( i in seq_along( all.vars )) {
  # 0.4 %  DMSO was used for propylthiouracil en captopril. The rest in 0.2% DMSO
  
  cur.cell_line <-unique( my.data.m[  cell.treat.dose %in% all.vars[ i ], cell_line ] )
  cur.treat <-unique( my.data.m[  cell.treat.dose %in% all.vars[ i ], treatment ] )
  subSet1 <- my.data.m[ cell.treat.dose %in% all.vars[ i ] & variable != "imageCountParentObj"]
  subSetDMSO <- my.data.m[ treatment %in% "DMSO" & cell_line %in% cur.cell_line & variable != "imageCountParentObj"]

  
  if(cur.treat == "DMSO") {
    subSet1[, treatment :=  paste(treatment, dose_uM) ]
    cur.treat <- unique(subSet1$treatment)
  
}
  i
  subSet <- rbind( subSetDMSO, subSet1)
  
  
  # one sided significance test
  subSet[ , treatment := factor(treatment, levels = c( "DMSO", cur.treat ))]
  #setnames(subSet, "V1", "meanDiff")
  
  if( nrow(subSet1) > 1) {
  stats.out[[i]] <- t.test(value ~ treatment, alternative = "two.sided", var.equal = FALSE, data =   subSet)
  
  m.out<- stats.out[[i]] 
  stats.out.data[[i]] <- data.frame(cell_line = unique(subSet$cell_line),
                                    treatment = unique(subSet1$treatment),
                                    dose_uM = unique(subSet1$dose_uM),
                                    
                                    diffMeans = m.out$estimate[2] - m.out$estimate[1], 
                                    pVal = m.out$p.value)
  } else 
  {
    stats.out[[i]] <- NULL
    stats.out.data <- NULL
  }
 
}

stats.out.data.df <- do.call('rbind', stats.out.data)


save.path <-  "C:/Users/steve_000/Documents/work/POC paper"
save.image(paste(save.path, "/PIPOC04-062016.Rdata", sep = ""))

write.table(stats.out.data.df, file = paste(save.path, "/statsPIdataPOCpaper.txt", sep = ""), sep = "\t")

my.data_s <- ddply( my.data.m, .(treatment, dose_uM, cell_line, variable), summarize,
                                                    m.value= mean(value, na.rm=TRUE),
                                                    SD = sd(value, na.rm = TRUE))


my.data_s<- as.data.table(my.data_s ) 




unique(my.data_s$variable)
unique(my.data_s$cell_line)
# create dose factors

my.data_s<-my.data_s[order(my.data_s$cell_line, my.data_s$variable, my.data_s$treatment, my.data_s$dose_uM ),]
lapply(my.data_s, class)
my.data_s$treatment<- as.character(my.data_s$treatment)



# create levels within each treatment for dose
counts.d <- ddply(my.data_s, .(treatment, cell_line, variable ), summarize,count.d.l = length(dose_uM))
head(counts.d)
my.data_s$dose.f <- NULL
my.data_s$dose.f <- NA
for (i in 1 : nrow(counts.d))
{
  
  my.data_s$dose.f[ my.data_s$treatment == counts.d$treatment[i] & 
                      my.data_s$variable == counts.d$variable[i] &
#  my.data_s$plateID == counts.d$plateID[i] &                    
  my.data_s$cell_line == counts.d$cell_line[i]] <- gl(counts.d$count.d.l[i], 1)
  
}

unique(my.data_s$variable)

my.data_s$dose.f<-factor(my.data_s$dose.f)

my.data_s[ ,variable:= gsub("count_pi_obj_masked_primaryId_AreaShape_Area.DIV.Nuclei_AreaShape_Area_larger_0.05", "count_pi_obj_0.05", variable)]

my.data_s$cell_line_plateID <- paste(my.data_s$cell_line, my.data_s$plateID)
unique(my.data_s$treatment)

#for( i in seq_along(cells.keep)){
  
#  my.data_s.cell <- my.data_s[ cell_line %in% cells.keep[i],]

  my.data_s <- my.data_s[ variable != "imageCountParentObj"]
  my.data_s <- my.data_s[ !treatment %in% c("H2O2", "Staurosporin")]
  
# dmso data ernaast plakken. verschil in means / pooled sds  
  dmso.d <- my.data_s[ treatment =="DMSO"]
  
#limits <- aes(ymax = m.value + 0.5*SD, ymin = m.value - 0.5*SD)
#my.data_s.cell
p <- ggplot(data = my.data_s, aes(x = treatment , y = m.value, fill = dose.f))  + geom_bar(stat = "identity", position = "dodge") +
    theme_set(theme_gray(base_size = 8)) 
p <- p +  theme( axis.text.x = element_text(angle = 90, hjust = 1, vjust =0.4, size = 8, 
                                            colour = "grey50") ) + theme( strip.text.x = element_text( )) +
  ggtitle( "Cell Viability" ) + 
  theme(plot.title = element_text(lineheight=.8, size = 14 ))

dodge <- position_dodge(width=0.9)

p <- p + facet_grid( ~cell_line_plateID , scales = "free_y" ) #+ #ylim(c(0,1))
  #geom_errorbar(limits, width = 0.05, position = dodge) 
p+ ylim(c(0,1))
pdf( file = paste(paste("figures2/", cells.keep[i], sep =""), ".pdf"), height = 8, width = 14)
print(p) #+ylim(c(0,1))
dev.off()
}
getwd()
