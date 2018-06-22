#R
#GJR 5/21/2018
#Eddy re-aligned around indels and added indels...analyzing the new data set as exported from mysql as:
# mysql -u root -p  PomUrbanaGrant -e "select * from snpFisherurbana" -B > poolseqSnpsUrbana.txt

################################################################################### 
############################ INITIAL DATA MANIPULATION ############################ 
################################################################################### 

#import, deterimine which rows are indels, and creates 2 new data sets, snps and indels
setwd('/media/raglandlab/ExtraDrive4/Urbana_PoolSeq')
data<-read.table('poolseqSnpsUrbana.txt',stringsAsFactors=F,row.names=NULL,header=T)
gt1<-function(x) {
    out<-F
    if (nchar(x) > 1) {out<-T}
    return(out)
}
isIndel<-sapply(data$ref,gt1) | sapply(data$alt,gt1)
snps<-data[!isIndel,]
indels<-data[isIndel,]
rm(data)
gc()

#fix problem with character columns
a<-apply(snps[,6:23],2,as.numeric)
snps[,6:23]<-a
rm(a)
gc()

#filter to min 10x per pop coverage
ind<-rowSums(snps[,6:7]) >= 10 & rowSums(snps[,8:9]) >= 10 
snps<-snps[ind,]

#estimate fdr for snps
snps$fdrHost<-p.adjust(snps$urbana_appleave_hawave_fisher_pvalue,method='BH')
snps$fdrEclApple<-p.adjust(snps$urbana_appleearly_applelate_fisher_pvalue,method='BH')
snps$fdrEclHaw<-p.adjust(snps$urbana_hawearly_hawlate_fisher_pvalue,method='BH')


################################################################################### 
################################################################################### 
################################################################################### 


############################################################################################################## 
######### Snp overrepresentation stats and allele frequencies, focusin on eclosion associations ############## 
############################################################################################################## 


#SNPs
############# 1.68X as many nominally significant differences b/t host races than expected at alpha = 0.05 #####################
sum(snps$urbana_appleave_hawave_fisher_pvalue < 0.05)/nrow(snps)
                                        #[1] 0.08418959 -- 8.4% sig at 5% Type I error

########### 2.5X times as many nominally significant differences b/t eclosion bulks in haw than expected at alpha = 0.05 ######
sum(snps$urbana_hawearly_hawlate_fisher_pvalue < 0.05)/nrow(snps)
#[1] 0.125069 -- 12.5% sig at 5% Type I error

########### 2.7X times as many nominally significant differences b/t eclosion bulks in haw than expected at alpha = 0.05 ######
sum(snps$urbana_appleearly_applelate_fisher_pvalue < 0.05)/nrow(snps)
#[1] 0.1361606 -- 13.6% sig at 5% Type I error

#Indels
############# 2X as many nominally significant differences b/t host races than expected at alpha = 0.05 #####################
sum(indels$urbana_appleave_hawave_fisher_pvalue < 0.05)/nrow(indels)
#[1] 0.1077189 -- 10.7% sig at 5% Type I error

########### 3.2X times as many nominally significant differences b/t eclosion bulks in haw than expected at alpha = 0.05 ######
sum(indels$urbana_hawearly_hawlate_fisher_pvalue < 0.05)/nrow(indels)
#[1] 0.1604757 -- 16% sig at 5% Type I error

########### 3.5X times as many nominally significant differences b/t eclosion bulks in haw than expected at alpha = 0.05 ######
sum(indels$urbana_appleearly_applelate_fisher_pvalue < 0.05)/nrow(indels)
#[1] 0.174041 -- 17.4% sig at 5% Type I error




#Calculate allele frequency diffs for snps, apple - haw and early - late
snps$pdifHost<-(snps$urbana_appleave_maj/(snps$urbana_appleave_maj+snps$urbana_appleave_min)) -
    (snps$urbana_hawave_maj/(snps$urbana_hawave_maj+snps$urbana_hawave_min))
snps$pdifEclHaw<-(snps$urbana_hawearly_maj/(snps$urbana_hawearly_maj+snps$urbana_hawearly_min)) -
    (snps$urbana_hawlate_maj/(snps$urbana_hawlate_maj+snps$urbana_hawlate_min))
snps$pdifEclApple<-(snps$urbana_appleearly_maj/(snps$urbana_appleearly_maj+snps$urbana_appleearly_min)) -
    (snps$urbana_applelate_maj/(snps$urbana_applelate_maj+snps$urbana_applelate_min))


# r = 0.96 correlation b/t allele frequency difs for haw ecl vs apple ecl
# but, r goes down to ~.65 - .69 when consiering only loci sig in apple (or haw)
ind<-snps$fdrEclApple < 0.05 & snps$fdrEclHaw < 0.05
cor.test(snps$pdifEclApple[ind],snps$pdifEclHaw[ind])

#but, non-sig correlation b/t ecl and host race for sig loci in both apple and haw
snps$meanEclFreqDif<-rowMeans(snps[,29:30])
ind<-snps$fdrEclHaw < 0.05 & snps$fdrHost < 0.05 & snps$fdrEclApple < 0.05
cor.test(snps$pdifHost[ind],snps$meanEclFreqDif[ind])

#slightly positive correlation for apple ecl and host race
ind<-snps$fdrEclApple < 0.05 & snps$fdrHost < 0.05
cor.test(snps$pdifHost[ind],snps$pdifEclApple[ind])

#but, stronger negative correlation b/t host race and haw ecl
ind<-snps$fdrEclHaw < 0.05 & snps$fdrHost < 0.05 
cor.test(snps$pdifHost[ind],snps$pdifEclHaw[ind])


sameSign<-function(x,y) {
   out<-(x*y) > 0
   return(out)
}

ind<-snps$fdrEclHaw < 0.05 & snps$fdrHost < 0.05 
isSame<-sameSign(snps$pdifHost[ind],snps$pdifEclHaw[ind])
#Only 34% of sig loci (Ecl haw, host race) go in expected direction
sum(isSame)/sum(ind)


ind<-snps$fdrEclApple < 0.05 & snps$fdrHost < 0.05 
isSame<-sameSign(snps$pdifHost[ind],snps$pdifEclApple[ind])
#Here, 60% of sig loci (Ecl apple, host race) go in expected direction
sum(isSame)/sum(ind)


a<-sum( snps$fdrEclHaw < 0.05 & snps$fdrHost < 0.05 )
b<-sum( !(snps$fdrEclHaw) < 0.05 & snps$fdrHost < 0.05 )
c<-sum( snps$fdrEclHaw < 0.05 & !(snps$fdrHost < 0.05) )
d<-sum( !(snps$fdrEclHaw) < 0.05 & !(snps$fdrHost < 0.05) )

#                   sigEcl nSigEcl
#       sigHost       a       b
#       nSigHost      c       d

fisher.test(rbind(c(a,b),c(c,d)))
# 2.8X more shared loci than expected by chance between ecl association and host association ##NOTE: OR of 31 if setting thres at 0.005
#	Fisher's Exact Test for Count Data

#data:  rbind(c(a, b), c(c, d))
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
# 2.556448 3.032937
#sample estimates:
#odds ratio 
#  2.786932 



#source generic function to query mysql databases
source('queryDb.r')

ind<-snps$fdrEclHaw < 0.005 &  snps$fdrEclApple < 0.005
ids<-snps$snpId[ind]
#query<-sprintf("select * from annotation where snpId = %s", id)
queryString<-"select * from annotation where snpId = %s"
geneList<-queryDb(ids,queryString)
locs<-geneList$loc[geneList$effect!='intergenic_region']
queryString<-"select * from feature_alias where loc ='%s'"
annoList<-queryDb(locs,queryString)
#candidates affecting eclosion timing
#write.table(annoList,'snpsInGenesFdrEclHawLt0.005FdrEclAppleLt0.005.txt',sep="\t",quote=F,row.names=F)




################################################################################### 
################################################################################### 
################################################################################### 


######################################################################################################## 
######### various analyses including indels, still focusing on ecl associations  ####################### 
######################################################################################################## 


#fix problem with character columns
a<-apply(indels[,18:23],2,as.numeric)
indels[,18:23]<-a
rm(a)
gc()

#filter to min 10x per pop coverage
getCounts<-function(x) {sum(as.numeric(unlist(strsplit(x,','))))}
ind<-( sapply(indels$urbana_appleave_maj,getCounts) + sapply(indels$urbana_appleave_min,getCounts) ) >= 10 &
    ( sapply(indels$urbana_hawave_maj,getCounts) + sapply(indels$urbana_hawave_min,getCounts) ) >= 10
a<-ind
a[is.na(a)]<-F
indels<-indels[a,]


#estimate fdr
indels$fdrHost<-p.adjust(indels$urbana_appleave_hawave_fisher_pvalue,method='BH')
indels$fdrEclApple<-p.adjust(indels$urbana_appleearly_applelate_fisher_pvalue,method='BH')
indels$fdrEclHaw<-p.adjust(indels$urbana_hawearly_hawlate_fisher_pvalue,method='BH')

indels$fdrHost[is.na(indels$fdrHost)]<-1
indels$fdrEclApple[is.na(indels$fdrEclApple)]<-1
indels$fdrEclHaw[is.na(indels$fdrEclHaw)]<-1


a<-sum( indels$fdrEclHaw < 0.05 & indels$fdrHost < 0.05 )
b<-sum( !(indels$fdrEclHaw) < 0.05 & indels$fdrHost < 0.05 )
c<-sum( indels$fdrEclHaw < 0.05 & !(indels$fdrHost < 0.05) )
d<-sum( !(indels$fdrEclHaw) < 0.05 & !(indels$fdrHost < 0.05) )

#                   sigEcl nSigEcl
#       sigHost       a       b
#       nSigHost      c       d

fisher.test(rbind(c(a,b),c(c,d)))
# 1.2X more shared loci than expected by chance between ecl association and host association
#	Fisher's Exact Test for Count Data

#data:  rbind(c(a, b), c(c, d))
#p-value = 0.001426
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
# 1.085663 1.406700
#sample estimates:
#odds ratio 
#  1.238157 

a<-sum( indels$fdrEclHaw < 0.01 & indels$fdrEclApple < 0.01 )
b<-sum( !(indels$fdrEclHaw) < 0.01 & indels$fdrEclApple < 0.01 )
c<-sum( indels$fdrEclHaw < 0.01 & !(indels$fdrEclApple < 0.01) )
d<-sum( !(indels$fdrEclHaw) < 0.01 & !(indels$fdrEclApple < 0.01) )

#                   sigEclHaw nSigEclHaw
#       sigEclAp       a       b
#       nSigEclAp      c       d

fisher.test(rbind(c(a,b),c(c,d)))

#	Fisher's Exact Test for Count Data
# 5.2X more shared loci than expected by chance between ecl association in haw vs. apple
#data:  rbind(c(a, b), c(c, d))
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
# 4.985778 5.558897
#sample estimates:
#odds ratio 
#  5.266769  


#source generic function to query mysql databases
source('queryDb.r')

ind<-indels$fdrEclHaw < 0.01 &  indels$fdrEclApple < 0.01
ids<-indels$snpId[ind]
#query<-sprintf("select * from annotation where snpId = %s", id)
queryString<-"select * from annotation where snpId = %s"
geneList<-queryDb(ids,queryString)
locs<-geneList$loc[geneList$effect!='intergenic_region']
queryString<-"select * from feature_alias where loc ='%s'"
annoList<-queryDb(locs,queryString)
#candidates affecting eclosion timing
#write.table(annoList,'indelsInGenesFdrEclHawLt0.01FdrEclAppleLt0.01.txt',sep="\t",quote=F,row.names=F)


#indels AND snps
fdrHost<-p.adjust(c(indels$urbana_appleave_hawave_fisher_pvalue,snps$urbana_appleave_hawave_fisher_pvalue),method='BH')
fdrEclApple<-p.adjust(c(indels$urbana_appleearly_applelate_fisher_pvalue,snps$urbana_appleearly_applelate_fisher_pvalue),method='BH')
fdrEclHaw<-p.adjust(c(indels$urbana_hawearly_hawlate_fisher_pvalue,snps$urbana_hawearly_hawlate_fisher_pvalue),method='BH')

fdrHost[is.na(fdrHost)]<-1
fdrEclApple[is.na(fdrEclApple)]<-1
fdrEclHaw[is.na(fdrEclHaw)]<-1


#get max snp positions for each scaffold
getMax<-function(scaf) {
    return( max( dat$position[ dat$scaffold==scaf ] ) )
}

uniScafs<-unique(snps$scaffold)
dat<-snps[,2:3]
library(parallel)
nCores=10
#uniScafs<-uniScafs[1:1000]
cl <- makeCluster(nCores)
clusterExport(cl=cl, varlist=c("getMax","dat","uniScafs"),envir=environment())
out<-unlist(parLapply(cl, uniScafs, getMax ))
stopCluster(cl)

scafSizes<-data.frame(scaffold=uniScafs,maxPos=out)



#save.image('UrbanPoolseq.Rdat')
setwd('/media/raglandlab/ExtraDrive4/Urbana_PoolSeq')
load('UrbanPoolseq.Rdat')


a<-sum( fdrEclHaw < 0.005 & fdrEclApple < 0.005 )
b<-sum( !(fdrEclHaw) < 0.005 & fdrEclApple < 0.005 )
c<-sum( fdrEclHaw < 0.005 & !(fdrEclApple < 0.005) )
d<-sum( !(fdrEclHaw) < 0.005 & !(fdrEclApple < 0.005) )

fisher.test(rbind(c(a,b),c(c,d)))
#12.5X more shared loci than expected by chance between haw and apple ecl association
#	Fisher's Exact Test for Count Data
#
#data:  rbind(c(a, b), c(c, d))
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
# 11.89312 12.60329
#sample estimates:
#odds ratio 
#  12.24394 

a<-sum( fdrEclHaw < 0.01 & fdrEclApple < 0.01 )
b<-sum( !(fdrEclHaw) < 0.01 & fdrEclApple < 0.01 )
c<-sum( fdrEclHaw < 0.01 & !(fdrEclApple < 0.01) )
d<-sum( !(fdrEclHaw) < 0.01 & !(fdrEclApple < 0.01) )

fisher.test(rbind(c(a,b),c(c,d)))
#10.3X more shared loci than expected by chance between haw and apple ecl association
#	Fisher's Exact Test for Count Data

#data:  rbind(c(a, b), c(c, d))
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
# 10.11910 10.53964
#sample estimates:
#odds ratio 
#  10.32781 

#source generic function to query mysql databases
source('queryDb.r')

ind<-fdrEclHaw < 0.01 & fdrEclApple < 0.01
snpIds<-c(indels$snpId,snps$snpId)
ids<-snpIds[ind]
#query<-sprintf("select * from annotation where snpId = %s", id)
queryString<-"select * from annotation where snpId = %s"
geneList<-queryDb(ids,queryString)
locs<-geneList$loc[geneList$effect!='intergenic_region']
queryString<-"select * from feature_alias where loc ='%s'"
annoList<-queryDb(locs,queryString)
#candidates affecting eclosion timing
#write.table(annoList,'SnpsAndIndelsInGenesFdrEclHawLt0.01FdrEclAppleLt0.01.txt',sep="\t",quote=F,row.names=F)




############### test enrichment for custom lists ############################
#not much going on here

#get all flybase ids from genome annotation associated with snps
## NOTE: may want to get all annotations, not just genes with SNPs, though feature_alias may provide that
select distinct(Flybase_FBgn) from feature_alias limit 10;
library(RMySQL)
db<-'PomUrbanaGrant'
mydb = dbConnect(MySQL(), user='raglandlab', password='pomonella', dbname=db)
query<-"select distinct(Flybase_FBgn) from feature_alias"
queryOb<-dbGetQuery(mydb, query)
lapply(dbListConnections( dbDriver( drv = "MySQL")), dbDisconnect)
allFlyIds<-queryOb[,1]
allFlyIds<-allFlyIds[!is.na(allFlyIds)]

source('/media/raglandlab/ExtraDrive1/RpomDiapauseRNAseqTraj_GJR/pomResultsOutliersRemoved/applyEnrichTests.r')
source('/media/raglandlab/ExtraDrive1/RpomDiapauseRNAseqTraj_GJR/pomResultsOutliersRemoved/applyEnrichTestsCustomOnly.r')

flyIds<-annoList$Flybase_FBgn
flyIds[is.na(flyIds)]<-'noHit'
a<-applyEnrichTestsCustomOnly(flyIds,allFlyIds)
a<-applyEnrichTests(flyIds,allFlyIds)

###############################################################


############### Test overrepresentation of Snps in DE RNAseq data #####################
############## also, print set of genes that are DE across time and significant for
############## GWAS with eclosion time in both apple and haw ####################

load('/media/raglandlab/ExtraDrive1/RpomDiapauseRNAseqTraj_GJR/pomResultsOutliersRemoved/clus.Rdat')
names(appledat)[11]<-'Lowest_FDR_from2M.apple.TimeSeries'
pomAll<-merge(appledat[,c(1:6,11)],hawdat[,c(1:5,10)])
names(pomAll)[2]<-'flyid'
pomAll$flyid[is.na(pomAll$flyid)]<-'noHit'
deFlyIds<-unique(pomAll$flyid[pomAll$Lowest_FDR_from2M.apple.TimeSeries < 0.05 & pomAll$Lowest_FDR_from.TimeSeries.haw < 0.05])
snpFlyIds<-unique(annoList$Flybase_FBgn)
snpFlyIds<-snpFlyIds[!is.na(snpFlyIds)]
commonFlyIds<-allFlyIds[allFlyIds %in% pomAll$flyid]

#              DE     notDE
#sigSNP         a        b
#notsigSNP      c        d
a<-sum( commonFlyIds %in% deFlyIds & commonFlyIds %in% snpFlyIds  )
b<-sum( !(commonFlyIds %in% deFlyIds) & commonFlyIds %in% snpFlyIds  )
c<-sum( commonFlyIds %in% deFlyIds & !(commonFlyIds %in% snpFlyIds)  )
d<-sum( !(commonFlyIds %in% deFlyIds) & !(commonFlyIds %in% snpFlyIds)  )
fisher.test(rbind(c(a,b),c(c,d)))
#1.3x more loci than expected by chance (or 20% more than expected) 
#  Fisher's Exact Test for Count Data
#data:  rbind(c(a, b), c(c, d))
#p-value = 0.0001441
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
# 1.137788 1.501256
#sample estimates:
#odds ratio 
#     1.308 


out<-commonFlyIds[commonFlyIds %in% deFlyIds & commonFlyIds %in% snpFlyIds]
write.table(out,'FlybaseDEtimeAndSigSnpFdrEclHawLt0.01FdrEclAppleLt0.01.txt',quote=F,row.names=F,sep="\t")
#all kinds of interesting enrichment in this list


## representation of SNPs relative to DE network statistics
deNet<-read.table('/media/raglandlab/ExtraDrive1/RpomDiapauseRNAseqTraj_GJR/pomResultsOutliersRemoved/DEtimeBothRaces.NetworkStats.txt',stringsAsFactors=F,sep="\t",row.names=NULL,header=T)
flyids<-unique(annoList$Flybase_FBgn)
geneids<-unique(annoList$gene_id)

ind<-deNet$id %in% geneids | deNet$id %in% flyids

median(deNet$hubScore[ind])
median(deNet$hubScore)
ks.test(deNet$hubScore,deNet$hubScore[ind])
#no clear enrichement of snps in more hub-like genes

median(deNet$degree[ind])
median(deNet$degree)
ks.test(deNet$degree,deNet$degree[ind])
#also no clear enrichment of snps in genes with higher degree


##### Next, proceed with correlation between allele freq difs and network statistics

a<-merge( geneList[,c(6,15)],annoList[,c(3:4,6)] )
a<-a[!duplicated(a$snpId),]
a<-a[a$snpId %in% snps$snpId,]
ind<-match(a$snpId,snps$snpId)
fdif<-rowMeans( cbind(snps$pdifEclApple[ind], snps$pdifEclHaw[ind]) )
dat<-data.frame(a,fdif=fdif,degree=rep(NA,nrow(a)),hubScore=rep(NA,nrow(a)))


for (i in 1:nrow(dat)) {
    if (dat$gene_id[i] %in% deNet$id) {
        dat[i,6:7]<-deNet[match(dat$gene_id[i],deNet$id) ,2:3]
    }
    if (sum(dat$Flybase_FBgn[i] %in% deNet$id)==1) {
         dat[i,6:7]<-deNet[match(dat$Flybase_FBgn[i],deNet$id) ,2:3]
    }
}

dat<-dat[!duplicated(dat$gene_id),]
cor.test(abs(dat$fdif),dat$degree)
cor.test(abs(dat$fdif),dat$hubScore)
# no correlation b/t allele freq dif and network stats


physInt<-read.table('physical_interactions_fb_2018_02.tsv',header=F,stringsAsFactors=F,row.names=NULL,sep="\t",quote="\"")
ids<-unique(physInt[,1])
degrees<-vector(length=length(ids))
for (i in 1:length(ids)) {
    degrees[i]<-sum(physInt[,1] %in% ids[i] )
}
flyids<-unique(annoList$Flybase_FBgn)


ind<-ids %in% flyids

mean(degrees)
mean(degrees[ind])
ks.test(degrees,degrees[ind])
#no enrichment of snps sig and de within degree as calculated from the physical interactions dataset imported above


#######################################################################################



################################################################################### 
################################################################################### 
################################################################################### 


######################################################################################################## 
######### SNP ANALYSIS FOCUSING ON HOST RACE COMPARISONS  ############################################## 
########################################################################################################

#some of this is redundant to code above, placed here just to keep things organized by questions


# r = 0.96 correlation b/t allele frequency difs for haw ecl vs apple ecl
# but, r goes down to ~.65 - .69 when consiering only loci sig in apple (or haw)
ind<-snps$fdrEclApple < 0.05 & snps$fdrEclHaw < 0.05
cor.test(snps$pdifEclApple[ind],snps$pdifEclHaw[ind])

#but, non-sig correlation b/t ecl and host race for sig loci in both apple and haw
snps$meanEclFreqDif<-rowMeans(snps[,29:30])
ind<-snps$fdrEclHaw < 0.05 & snps$fdrHost < 0.05 & snps$fdrEclApple < 0.05
cor.test(snps$pdifHost[ind],snps$meanEclFreqDif[ind])

#slightly positive correlation for apple ecl and host race
ind<-snps$fdrEclApple < 0.05 & snps$fdrHost < 0.05
cor.test(snps$pdifHost[ind],snps$pdifEclApple[ind])

#but, stronger negative correlation b/t host race and haw ecl
ind<-snps$fdrEclHaw < 0.1 & snps$fdrHost < 0.1 
cor.test(snps$pdifHost[ind],snps$pdifEclHaw[ind])


sameSign<-function(x,y) {
   out<-(x*y) > 0
   return(out)
}



ind<-snps$fdrEclHaw < 0.01 & snps$fdrHost < 0.01 
isSame<-sameSign(snps$pdifHost[ind],snps$pdifEclHaw[ind])
#Only 34% of sig loci (Ecl haw, host race) go in expected direction
sum(isSame)/sum(ind)


ind<-snps$fdrEclApple < 0.05 & snps$fdrHost < 0.05 
isSame<-sameSign(snps$pdifHost[ind],snps$pdifEclApple[ind])
#Here, 60% of sig loci (Ecl apple, host race) go in expected direction
sum(isSame)/sum(ind)


ind<-snps$fdrEclHaw < 0.05 & snps$fdrHost < 0.05 & snps$fdrEclApple < 0.05
isSame<-sameSign(snps$pdifHost[ind],snps$meanEclFreqDif[ind])
#Here, 54% of sig loci (Ecl apple and Haw, host race) go in expected direction
sum(isSame)/sum(ind)


#indels and snps 
a<-sum( fdrEclHaw < 0.005 & fdrHost < 0.005 )
b<-sum( !(fdrEclHaw) < 0.005 & fdrHost < 0.005 )
c<-sum( fdrEclHaw < 0.005 & !(fdrHost < 0.005) )
d<-sum( !(fdrEclHaw) < 0.005 & !(fdrHost < 0.005) )

fisher.test(rbind(c(a,b),c(c,d)))
#strong overrep of ecl associated snps in host associated snps
#	Fisher's Exact Test for Count Data
#data:  rbind(c(a, b), c(c, d))
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
# 19.26105 37.82126
#sample estimates:
#odds ratio 
#  27.33944 
#NOTE: OR = 2.0, p < 2.2e-16 with ind<-fdrEclHaw < 0.005 & fdrHost < 0.005


#           sigBoth sigEclApple  sigEclHaw nsigBoth
#sigHost      a        b            c         d
#nsigHost     e        f            g         h


#test whether proportion of snps contributing to host race differences is different for Ecl-associated loci in both apple and haw
a<-sum(fdrEclHaw < 0.05 & fdrEclApple < 0.05 & fdrHost < 0.05)
b<-sum(fdrEclHaw > 0.05 & fdrEclApple < 0.05 & fdrHost < 0.05)
c<-sum(fdrEclHaw < 0.05 & fdrEclApple > 0.05 & fdrHost < 0.05)
d<-sum(fdrEclHaw > 0.05 & fdrEclApple > 0.05 & fdrHost < 0.05)

e<-sum(fdrEclHaw < 0.05 & fdrEclApple < 0.05 & fdrHost > 0.05)
f<-sum(fdrEclHaw > 0.05 & fdrEclApple < 0.05 & fdrHost > 0.05)
g<-sum(fdrEclHaw < 0.05 & fdrEclApple > 0.05 & fdrHost > 0.05)
h<-sum(fdrEclHaw > 0.05 & fdrEclApple > 0.05 & fdrHost > 0.05)

library(mosaic)
xchisq.test(cbind(c(a,b,c,d),c(e,f,g,h)))
# all sig-sig loci overrepresented, especially sigBoth (no clear pattern that, e.g. ecl-associated loci in apple tend to be more host-associated)






source('queryDb.r')

ind<-fdrEclHaw < 0.1 & fdrHost < 0.1
snpIds<-c(indels$snpId,snps$snpId)
ids<-snpIds[ind]
#query<-sprintf("select * from annotation where snpId = %s", id)
queryString<-"select * from annotation where snpId = %s"
geneList<-queryDb(ids,queryString)
locs<-geneList$loc[geneList$effect!='intergenic_region']
queryString<-"select * from feature_alias where loc ='%s'"
annoList<-queryDb(locs,queryString)



load('/media/raglandlab/ExtraDrive1/RpomDiapauseRNAseqTraj_GJR/pomResultsOutliersRemoved/clus.Rdat')
names(appledat)[11]<-'Lowest_FDR_from2M.apple.TimeSeries'
pomAll<-merge(appledat[,c(1:6,11)],hawdat[,c(1:5,10)])
names(pomAll)[2]<-'flyid'
pomAll$flyid[is.na(pomAll$flyid)]<-'noHit'
deFlyIds<-unique(pomAll$flyid[pomAll$Lowest_FDR_from2M.apple.TimeSeries < 0.05 & pomAll$Lowest_FDR_from.TimeSeries.haw < 0.05])
snpFlyIds<-unique(annoList$Flybase_FBgn)
snpFlyIds<-snpFlyIds[!is.na(snpFlyIds)]
commonFlyIds<-allFlyIds[allFlyIds %in% pomAll$flyid]

#              DE     notDE
#sigSNP         a        b
#notsigSNP      c        d
a<-sum( commonFlyIds %in% deFlyIds & commonFlyIds %in% snpFlyIds  )
b<-sum( !(commonFlyIds %in% deFlyIds) & commonFlyIds %in% snpFlyIds  )
c<-sum( commonFlyIds %in% deFlyIds & !(commonFlyIds %in% snpFlyIds)  )
d<-sum( !(commonFlyIds %in% deFlyIds) & !(commonFlyIds %in% snpFlyIds)  )
fisher.test(rbind(c(a,b),c(c,d)))
#1.5x more loci than expected by chance (30% more)
#	Fisher's Exact Test for Count Data

#data:  rbind(c(a, b), c(c, d))
#p-value = 3.67e-07
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
# 1.286513 1.761454
#sample estimates:
#odds ratio 
#  1.506959 

source('/media/raglandlab/ExtraDrive1/RpomDiapauseRNAseqTraj_GJR/pomResultsOutliersRemoved/applyEnrichTests.r')
source('/media/raglandlab/ExtraDrive1/RpomDiapauseRNAseqTraj_GJR/pomResultsOutliersRemoved/applyEnrichTestsCustomOnly.r')

flyIds<-annoList$Flybase_FBgn
flyIds[is.na(flyIds)]<-'noHit'
a<-applyEnrichTestsCustomOnly(flyIds,allFlyIds)
# suh and wnt sig at FDR 0.05949280 and 0.06836221, respectively
a<-applyEnrichTests(flyIds,allFlyIds)
#nothing in Kegg

#XXX
out<-snpFlyIds[snpFlyIds %in% deFlyIds]
#write.table(out,'SnpsAndIndelsInGenesFdrEclHawLt0.1FdrHostLt0.1.txt',sep="\t",quote=F,row.names=F)  
#Genes DE b/t races and different among host races and associated with ecl enriched for:
#  immunoglobulin subtype, zinc finger, Wnt signaling, chromatin regulator, nucleus, developmental protein, coiled coil, alternative splicing
# in wnt, some same dir genes DE, some 'opposite' genes DE, doesn't seem to form much of a pattern

source('queryDb.r')

ind<-snps$fdrEclHaw < 0.0005 | snps$fdrHost < 0.0005 | snps$fdrEclApple < 0.0005
#ind<-snps$fdrEclHaw < 0.001 | snps$fdrHost < 0.001 
dat<-snps[ind,]
ids<-dat$snpId
#query<-sprintf("select * from annotation where snpId = %s", id)
queryString<-"select * from annotation where snpId = %s"
geneList<-queryDb(ids,queryString)
locs<-geneList$loc[geneList$effect!='intergenic_region']
queryString<-"select * from feature_alias where loc ='%s'"
annoList<-queryDb(locs,queryString)
adat<-geneList[geneList$effect!='intergenic_region',c(6,15)]
adat<-adat[!duplicated(adat$snpId),]
ldat<-annoList[!duplicated(annoList$loc),c(4,6)]
a<-merge(adat,ldat)
b<-merge(dat,a)
isSame<-sameSign(b$pdifHost,b$pdifEclHaw)
b$isSame<-rep(0,nrow(b))
b$isSame[isSame]<-1
#same and opposite sign genes
#write.table(b,'SnpsInGenesFdrEclHawLt0.1FdrHostLt0.1.sameAndOppositeSign.txt',sep="\t",quote=F,row.names=F)                                       
#opposite sign genes enriched for membrane (Receptors) and zinc finger (others indicating transcription), plus alternative splicing
#same sign enriched for rho gtpase, **Wnt**, membrane, alternative splicing, neg reg of wnt, developmental protein

#re-polarize to most frequent snp in haw
ind<-dat$urbana_hawave_maj/(dat$urbana_hawave_min+dat$urbana_hawave_maj) < dat$urbana_appleave_maj/(dat$urbana_appleave_min+dat$urbana_appleave_maj)
dat$pdifHost[ind]=-dat$pdifHost[ind]
dat$pdifEclHaw[ind]=-dat$pdifEclHaw[ind]
dat$pdifEclApple[ind]=-dat$pdifEclApple[ind]





ids<-unique(geneList$scaffold)
queryString<-"select * from RAD_linkage_matchPool where scaffold = '%s'"
chromList<-queryDb(ids,queryString)
c<-chromList[!duplicated(chromList$scaffold),c(2,6,7)]
d<-merge(dat,c)


mydb = dbConnect(MySQL(), user='raglandlab', password='pomonella', dbname=db)
ldChroms<-dbGetQuery(mydb, 'select   * from     RAD_linkage_matchPool')
lapply(dbListConnections( dbDriver( drv = "MySQL")), dbDisconnect)


source('queryDb.closeSnp.r')

ind<-!is.na(ldChroms$scaffold) & !is.na(ldChroms$position)
ldChroms<-ldChroms[ind,]

ind<-snps$scaffold==scaf & snps$position >= (pos-thresh) & snps$position <= (pos+thresh)

j=1
thresh=200
for (i in 1:nrow(ldChroms)) {
#for (i in 1:50) {
    ind<-snps$scaffold==ldChroms$scaffold[i] & snps$position >= (ldChroms$position[i]-thresh) & snps$position <= (ldChroms$position[i]+thresh)
    if (j==1 & sum(ind) > 0) {
        out<-snps[ind,]
        out$chrom<-ldChroms$chromosome[i]
        out$LDgr<-ldChroms$LDgr[i]
    }
    if (j > 1 & sum(ind) > 0) {
        newRows<-snps[ind,]
        newRows$chrom<-ldChroms$chromosome[i]
        newRows$LDgr<-ldChroms$LDgr[i]
        out<-rbind(out,newRows)
    }
    j=j+1
}

ind<-out$chrom==2 & out$LDgr=='H6'
cor.test(out$pdifHost[ind],out$pdifEclHaw[ind])


closeSnps<-queryDb.closeSnp(ldChroms$scaffold[ind],ldChroms$position[ind],Table='main')


a<-queryDb.closeSnp(d$scaffold,d$position)

d<-d[,-33]
e<-merge(d,a[,c(2:3,7)])




#can also use above to correlate host differences with ecl differences per chromosome, have to adjust 'ind' above to get reasonable numbers
# results qualitatively similar for different ind thresholds, and when using the average freq. diff of eclosion bulks b/t apple and haw
#below is for ind<-snps$fdrEclHaw < 0.0005 | snps$fdrHost < 0.0005 | snps$fdrEclApple < 0.0005
ind<-d$chromosome==1
cor.test(d$pdifHost[ind],d$pdifEclHaw[ind]) # r = -0.2802802 df=7576, p < 2.2e-16
ind<-d$chromosome==2
cor.test(d$pdifHost[ind],d$pdifEclHaw[ind]) # r = -0.05025086  df=726, p = 0.1756
ind<-d$chromosome==3
cor.test(d$pdifHost[ind],d$pdifEclHaw[ind]) # r = 0.1399554  df=7294, p < 2.2e-16
ind<-d$chromosome==4
cor.test(d$pdifHost[ind],d$pdifEclHaw[ind]) # r = -0.07972383  df=278, p = 0.1835
ind<-d$chromosome==5
cor.test(d$pdifHost[ind],d$pdifEclHaw[ind]) # r = -0.0420958  df=1175, p = 0.1489


e<-d[!duplicated(d$loc),]

#for ind<-snps$fdrEclHaw < 0.1 & snps$fdrHost < 0.1
table(e$chromosome[e$isSame==0])
table(e$chromosome[e$isSame==1])
chisq.test(cbind(table(e$chromosome[e$isSame==0]),table(e$chromosome[e$isSame==1])))
#distribution on chroms is really different
#write.table(e,'SnpsInGenesFdrEclHawLt0.1FdrHostLt0.1.sameAndOppositeSign.chromAssign.txt',sep="\t",quote=F,row.names=F) 

table(e$LDgr[e$isSame==0])
table(e$LDgr[e$isSame==1])


## representation of SNPs relative to DE network statistics
deNet<-read.table('/media/raglandlab/ExtraDrive1/RpomDiapauseRNAseqTraj_GJR/pomResultsOutliersRemoved/DEtimeBothRaces.NetworkStats.txt',stringsAsFactors=F,sep="\t",row.names=NULL,header=T)
flyids<-unique(annoList$Flybase_FBgn)
geneids<-unique(annoList$gene_id)

ind<-deNet$id %in% geneids | deNet$id %in% flyids

median(deNet$hubScore[ind])
median(deNet$hubScore)
ks.test(deNet$hubScore,deNet$hubScore[ind])
#no clear enrichement of snps in more hub-like genes

median(deNet$degree[ind])
median(deNet$degree)
ks.test(deNet$degree,deNet$degree[ind])
#no clear enrichement of snps in genes with higher connectivity





##### XX new section applyin Lynch's LR test for frequency differences
source('LRfreqTest.r')

pval<-LRfreqTest(as.numeric(c(dat[3,14:17],48)))

subDat<-snps[,14:17]
subDat$N<-48

pvals<-apply(subDat[1:1000,],1,LRfreqTest)



#### XX permute correlations




permFdif<-function(x) {
    tot1<-sum(x[1:4])
    tot2<-sum(x[5:8])
    f1<-(x[1]+x[3])/tot1
    f2<-(x[5]+x[7])/tot2
    counts1<-tabulate(sample(1:4,tot1,replace=T,prob=rep(c(f1,(1-f1)),2)),nbins=4)
    counts2<-tabulate(sample(1:4,tot2,replace=T,prob=rep(c(f2,(1-f2)),2)),nbins=4)
    fdif<-counts1[1]/(counts1[1]+counts1[2]) - counts1[3]/(counts1[3]+counts1[4])
    fdif[2]<-counts2[1]/(counts2[1]+counts2[2]) - counts2[3]/(counts2[3]+counts2[4])
    return(fdif)
}



for (i in 1:nIter) {
    fdifs<-apply(snps[1:10,6:13],1,permFdif)
    
}

### XX should be able to run these separately in reasonable amounts of time
### list: HostvHawEcl (allSig) HawEclvAppleEcl (allSig) HostvHawEcl (by Chromosome)

nIter=10000
cors<-vector(length=nIter)

ind<-snps$fdrHost < 0.05 & snps$fdrEclHaw < 0.05
dat<-snps[ind,c(6:9,14:17)]

system.time(
for (i in 1:nIter) {
    fdifs<-apply(dat[1:50,],1,permFdif)
    return( cor(fdifs[1,],fdifs[2,]) )
}
)


nIter=100


#permutation test for correlation of allele frequency differences
# permutes reads assuming all come from one sampled population, in proportion to the overall frequencies of the alleles
# test is 2-tailed
permCor<-function(dat,processors=2,nIter=10000,permFdif.=permFdif) {
    corEst<-cor(dat[,9],dat[,10])
    library(doParallel)
    cl <- makeCluster(processors)
    registerDoParallel(cl)
    out<-foreach(i=1:nIter,.combine='c') %dopar% {
        fdifs<-apply(dat,1,permFdif.)
        return( cor(fdifs[1,],fdifs[2,]) )
    }
    stopCluster(cl)
    if (corEst > 0) {return(2*(1-ecdf(out)(corEst)))} else {return(2*ecdf(out)(corEst))}
}


ind<-snps$fdrHost < 0.05 & snps$fdrEclHaw < 0.05
dat<-snps[ind,c(6:9,14:17,28:29)]
permCor(dat)
#output = 0 (i.e., p < 1e-4 of getting a correlation stronger than the estimate by chance alone


ind<-snps$fdrEclHaw < 0.001 & snps$fdrEclApple < 0.001
dat<-snps[ind,c(6:9,14:17,29:30)]
permCor(dat)
#output = 0 (i.e., p < 1e-4 of getting a correlation stronger than the estimate by chance alone



sort(table(snps$scaffold[ind]))

snps$position[ind & snps$scaffold=='NW_016157145.1']

snps$position[ind & snps$scaffold=='NW_016174872.1']

snps$pdifHost[ind & snps$scaffold=='NW_016157420.1']

snps$position[ind & snps$scaffold=='NW_016157420.1']

ind2<-snps$scaffold=='NW_016157420.1' & snps$position <= 122303 & snps$position >= (122303-100)
hist(snps$pdifHost[ind2])




#############


ind<-snps$fdrHost < 0.05 & snps$fdrEclHaw < 0.05
dat<-snps[ind,]

scafPos<-list()
singletonScafs<-vector()
for (i in unique(dat$scaffold)) {
    sind<-dat$scaffold==i
    if (sum(sind) > 1) {
        scafPos[[i]]<-dat$position[sind]
    } else {singletonScafs<-c(singletonScafs,i)}
}

scafRange<-data.frame(scaffold=names(scafPos),range=0)
for (i in 1:nrow(scafRange)) {
    scaffold<-scafRange$scaffold[i]
    scafRange$range[i]<-max(scafPos[[scaffold]]) - min(scafPos[[scaffold]])
}

scafSizes$scaffold<-as.character(scafSizes$scaffold)

range<-max(scafRange$range)
scaf<-sample(scafSizes$scaffold[scafSizes$maxPos > range],1)
subDat<-snps[snps$scaffold==scaf,]

sample(subDat$pos[subDat$pos < ( max(subDat$pos) - range )],1)

getClose<-function(x,y) { x[which.min(abs(x-y))] }


sampleSameQuantiles<-function(realSnps,sampleSnps) {
    getClose<-function(y,x) { x[which.min(abs(x-y))] }
    sampMin<-min(sampleSnps)
    sampMax<-max(sampleSnps)
    ecdfReal<-ecdf(min(realSnps):max(realSnps))
    probs<-ecdfReal(realSnps)
    quants<-quantile(sampMin:sampMax,probs)
    return( sapply(quants,getClose,x=sampleSnps) )
}

i="NW_016157092.1"
scafSnps<-scafPos[[i]]

i=1
range<-scafRange$range[i]
scafSnps<-scafPos[[ scafRange$scaffold[i] ]]

#need to specify ind and export to cluster
scol1<-28 #column with first frequency diference in snps data frame
scol2<-29 #column with secend freq dif
sampleVariants<-function(scafSnps,range) {
     
    scaf<-sample(scafSizes$scaffold[scafSizes$maxPos > range & scafSizes$scaffold %in% scafAvail],1)
    subDat<-snps[snps$scaffold==scaf,c(3,scol1,scol2)]
    if (nrow(subDat)==1 | (max(subDat$position) - min(subDat$position)) < range) {
        out<-data.frame(scaffold="none",rep(0,length(scafSnps)),rep(0,length(scafSnps)))
    } else {
        start<-sample(subDat$position[subDat$position < ( max(subDat$position) - range )],1)
        end<-start + range
        sind<-subDat$position >= start & subDat$position <= end
        if (sum(sind) < length(scafSnps) | (max(subDat$position) - min(subDat$position)) < range) {
            out<-data.frame(scaffold="none",rep(0,length(scafSnps)),rep(0,length(scafSnps)))
        } else {
            sampPos<-sampleSameQuantiles(scafSnps,subDat$position[sind])
            sind<-match(sampPos,subDat$position)
            out<-data.frame( scaffold=scaf,subDat[sind,2], subDat[sind,3] ) 
        }
    }
    return(out)
    
}



permuteCor<-function(x) {

    fun<-function() {
        scafAvail<-uniScafs
        #XX also need to shuffle to permute read mapping to which pop, as above

        
        freqVals<-matrix(nrow=(length(unlist(scafPos)) + length(singletonScafs) ),ncol=2)
        #freqVals<-c(0,0)
        nSnps<-vector()
        j=1
        for (i in 1:length(scafPos)) {
            range<-scafRange$range[i]
            scafSnps<-scafPos[[ scafRange$scaffold[i] ]]
            sampledScaf="none"
            while (sampledScaf=="none") {
                vals<-sampleVariants(scafSnps,range)
                sampledScaf<-vals[1,1]
            }
            scafAvail<-scafAvail[scafAvail!=sampledScaf]
            incr<-nrow(vals)
            freqVals[j:(j+incr - 1),]<-as.matrix(vals[,2:3])
            j=j+incr
            #freqVals<-rbind(freqVals,as.matrix(vals[,2:3]))
            nSnps<-c(nSnps,nrow(vals))
        }
        start<-j
        for (j in start:nrow(freqVals)) {
            scaf<-sample(scafSizes$scaffold[scafSizes$scaffold %in% scafAvail],1)
            subDat<-snps[snps$scaffold==scaf,c(3,scol1,scol2)]
            if (length(subDat$position)==1) {
                pos<-subDat$position
            } else {
                pos<-sample(subDat$position,1)
                   }
            freqVals[j,]<-as.numeric(subDat[subDat$position==pos,2:3])
            scafAvail<-scafAvail[scafAvail!=scaf]
        }
        return(freqVals)
    }
    done<-0
    while (done==0) {
        out<-try(fun(),silent=T)
        if ( !is(out,"try-error") ) {done<-1}
    }
    
    return(cor(out[,1],out[,2]))
    #return(out)
    
}
ptm <- proc.time()
#a<-try(permuteCor(),silent=T)
a<-permuteCor(1)
proc.time() - ptm

if ( is(a,"try-error")


for(i in 1:100) {

    scaf<-sample(scafSizes$scaffold[scafSizes$maxPos > range & scafSizes$scaffold %in% scafAvail],1)
    subDat<-snps[snps$scaffold==scaf,c(3,scol1,scol2)]
    if (nrow(subDat)==1 | (max(subDat$position) - min(subDat$position)) < range) {
        out<-data.frame(scaffold="none",rep(0,length(scafSnps)),rep(0,length(scafSnps)))
    } else {
        start<-sample(subDat$position[subDat$position < ( max(subDat$position) - range )],1)
        end<-start + range
        sind<-subDat$position >= start & subDat$position <= end
        if (sum(sind) < length(scafSnps) | (max(subDat$position) - min(subDat$position)) < range) {
            out<-data.frame(scaffold="none",rep(0,length(scafSnps)),rep(0,length(scafSnps)))
        } else {
            sampPos<-sampleSameQuantiles(scafSnps,subDat$position[sind])
            sind<-match(sampPos,subDat$position)
            out<-data.frame( scaffold=scaf,subDat[sind,2], subDat[sind,3] ) 
        }
    }

}

nnSnps<-vector()
for (i in 1:length(scafPos)) {
nnSnps<-c(nnSnps,(length(scafPos[[i]])))
    
}
