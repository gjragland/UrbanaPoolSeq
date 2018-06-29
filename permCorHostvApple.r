#R
#GJR 6/22/2018
#permute snp positions to estimate approximate null distribution of a correlation between allele frequncy differences
#takes into account physical linkage by sampling physically close SNPs (same scaffold) whose quantiles (base pair position) approximate that of the focal SNPs
#samples random snps from scaffolds for singletons (sinlge snps on a scaffold)
#R CMD BATCH permCor.r

#setwd('/media/raglandlab/ExtraDrive4/Urbana_PoolSeq')
setwd('/scratch/summit/grra3428/poolseq')
load('UrbanPoolseq.Rdat')



ind<-snps$fdrEclApple < 0.05 & snps$fdrHost < 0.05
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






sampleSameQuantiles<-function(realSnps,sampleSnps) {
    getClose<-function(y,x) { x[which.min(abs(x-y))] }
    sampMin<-min(sampleSnps)
    sampMax<-max(sampleSnps)
    ecdfReal<-ecdf(min(realSnps):max(realSnps))
    probs<-ecdfReal(realSnps)
    quants<-quantile(sampMin:sampMax,probs)
    return( sapply(quants,getClose,x=sampleSnps) )
}

sampleVariants<-function(scafSnps,range) {
     
    scaf<-sample(scafSizes$scaffold[scafSizes$maxPos > range & scafSizes$scaffold %in% scafAvail],1)
    subDat<-snps[snps$scaffold==scaf,2:4]
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
        #scafAvail<-uniScafs
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
            subDat<-snps[snps$scaffold==scaf,2:4]
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

    #out<-fun()
    done<-0
    while (done==0) {
        out<-try(fun(),silent=T)
        if ( !is(out,"try-error") ) {done<-1}
    }
    
    return(cor(out[,1],out[,2]))
    #return(out)
    
}

#need to specify ind and export to cluster
scol1<-28 #column with first frequency diference in snps data frame
scol2<-30 #column with secend freq dif
snps<-snps[,c(2:3,scol1,scol2)]




#scafAvail<-uniScafs
#ptm <- proc.time()
#a<-try(permuteCor(),silent=T)
#a<-permuteCor(1)
#proc.time() - ptm

scafAvail<-uniScafs
library(parallel)
nCores=24
nIter=10000
cl <- makeCluster(nCores)
##need to export variable/functions so that they are visible in the virtual machine
clusterExport(cl=cl, varlist=c("sampleVariants","sampleSameQuantiles","scafPos",
                               "scafRange","scafSizes","singletonScafs","snps","uniScafs","permuteCor","scafAvail"),envir=environment())
ptm <- proc.time()
corVec<-unlist(parLapply(cl, 1:nIter, permuteCor ) )
proc.time() - ptm
stopCluster(cl)

write.table(corVec,'corVecHostvAppleEcl.txt',quote=F,row.names=F)
