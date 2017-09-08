## R script to process, analyze, and visualize genetic knockouts in C. albicans experiment 


# load libraries
rm(list = ls())
library(amap)
library(gplots)
library("extrafont")
#font_import(recursive = FALSE)
#loadfonts(device = "pdf")

# load and process index for data - see example index file 
#setwd('/path/to/data/index/')
setwd('/Users/cporter/Google Drive/Collins Lab Docs/Research Projects/Rebecca/publicCode')
geneIndex = read.table('geneIndex.csv', sep=",")
lowerIndex = as.matrix(geneIndex[lower.tri(geneIndex, diag=TRUE)])
upperIndex = as.matrix(geneIndex[upper.tri(geneIndex, diag=TRUE)])
upperIndex_i = upperIndex

# load index to gene name mapping
#geneNames <- read.table('geneNames.csv', sep=",")
geneNames <- read.table('/Users/cporter/Google Drive/Collins Lab Docs/Research Projects/Rebecca/publicCode/geneNames.csv', sep=",")

# reverse the string order of the name for upper genes
strReverse <- function(x)
        sapply(lapply(strsplit(x, NULL), rev), paste, collapse="")
upperIndex = strReverse(upperIndex)

# get index of aligned gene combinations 
sortedLower = sort(lowerIndex, index.return=TRUE)
sortedUpper = sort(upperIndex, index.return=TRUE)

# load and process data: should be a folder with .csv files of 
# OD measurements for each condition; see example file
#setwd('/path/to/OD/data/')
setwd('/Users/cporter/Google Drive/Collins Lab Docs/Research Projects/Rebecca/publicCode/data')
conditions <- list.files()

# define variables 
noReps <- 8 # user defined
noGenes <- nrow(geneNames)
noCond <- length(conditions)
noStrains <- length(lowerIndex)

# initialize matrices
allData <- matrix(data=NA, nrow=noStrains, ncol=noReps*noCond)
names <- matrix(data=NA,nrow=1,ncol=noReps*noCond)
avgData <- matrix(data=NA, nrow=noStrains, ncol=noCond)
avgLower <- matrix(data=NA, nrow=noStrains, ncol=noCond)
avgUpper <- matrix(data=NA, nrow=noStrains, ncol=noCond)

# assign data to appropratiate matrices 
# colums will be conditions, rows will be strains
for (i in 1:noCond){
        # pull out data for current condition 
        currentData <- read.table(conditions[i], sep=",")
        # separate replicates 
        controls <- currentData[,(noGenes+1):ncol(currentData)]
        knockouts <- currentData[,1:noGenes]
        # save only rows from file that contain data 
        knockouts <- knockouts[complete.cases(knockouts),]
        controls <- controls[complete.cases(controls),]
        
        # initialize temporary matrices to store current data 
        allCurrentData <- matrix(data=NA, nrow=noStrains, ncol=noReps)
        currentNames <- matrix(data=NA, nrow=1, ncol=noReps)
        for (j in 1:(noReps/2)){
                # pull out replicate data 
                repData <- knockouts[(j*noGenes-(noGenes-1)):(j*noGenes),]
                # pull out control data and take average
                ctrlData <- controls[(j*2-1):(j*2),]
                ctrlAvg <- mean(as.matrix(ctrlData), na.rm=TRUE)
                # pull out data for forward and reverse KO cases 
                lower <- as.matrix(repData[lower.tri(repData, diag = TRUE)])
                upper <- as.matrix(repData[upper.tri(repData, diag = TRUE)])
                # average strain matches
                orderedLower <- as.matrix(lower[sortedLower$ix])
                orderedUpper <- as.matrix(upper[sortedUpper$ix])
                allCurrentData[,(j*2-1)] <- orderedLower/ctrlAvg
                allCurrentData[,(j*2)] <- orderedUpper/ctrlAvg
                
                # remove .csv from file names for conditions 
                currentNames[(j*2-1):(j*2)] <- gsub(".csv", "", conditions[i])
        }
        # fill final matrices and build average data matrices
        allData[,((i*noReps-(noReps-1)):(i*noReps))] <- allCurrentData
        avgData[,i] <- rowMeans(allCurrentData, na.rm=TRUE)
        names[,((i*noReps-(noReps-1)):(i*noReps))] <- currentNames
        avgLower[,i] <- rowMeans(orderedLower/ctrlAvg, na.rm=TRUE)
        avgUpper[,i] <- rowMeans(orderedUpper/ctrlAvg, na.rm=TRUE)
}

# name rows and columns 
geneList <- sortedLower$x
geneList_i <- geneList
singles <- diag(as.matrix(geneIndex))
singlesIndex <- match(singles, geneList)
geneList <- paste0(paste0('-', geneList),'-')
for (i in 1:nrow(geneNames)){
        geneList <- gsub(paste0('^-', geneNames[i,1]), geneNames[i,2], geneList)
        geneList <- gsub(paste0(geneNames[i,1], '-$'), geneNames[i,2], geneList)
}
geneList <- gsub('-', '', geneList)



# rename columns and rows
avgData <- as.matrix(avgData)
avgLower <- as.matrix(avgLower)
avgUpper <- as.matrix(avgUpper)
rownames(allData)<-geneList
colnames(allData)<-names
rownames(avgData)<-geneList
colnames(avgData)<-gsub(".csv", "", conditions)
rownames(avgLower)<-sortedLower$x
colnames(avgLower)<-gsub(".csv", "", conditions)
rownames(avgUpper)<-upperIndex_i[sortedUpper$ix]
colnames(avgUpper)<-gsub(".csv", "", conditions)

# remove empty columns (for conditions that had fewer reps)
avgData_i <- avgData
allData <- allData[,complete.cases(t(allData))]
avgData <- avgData[,complete.cases(t(avgData))]
completeConditions <- conditions[complete.cases(t(avgData_i))]


# save data to csv files
#setwd('/path/to/save/data/')
setwd('/Users/cporter/Google Drive/Collins Lab Docs/Research Projects/Rebecca/publicCode')
write.csv(allData, "allData.csv")
write.csv(avgData, "avgData.csv")

# make heatmaps of "raw" data in plate shape
triangleHM <- function(triangleData, triangleType, colors, extension){
        tDataCol <- ncol(triangleData)
        for (i in 1:tDataCol){

                triangleGeneIndex <- geneIndex
                rownames(triangleGeneIndex) <- geneNames[,2]
                colnames(triangleGeneIndex) <- geneNames[,2]
                
                # reverse the gene index place holder for half triangle before sorting
                if (triangleType=="half"){
                        
                }
                
                # sort the data
                sortRow <- sort(rownames(triangleGeneIndex), index.return=TRUE)
                sortCol <- sort(colnames(triangleGeneIndex), index.return=TRUE)
                triangleGeneIndex <- triangleGeneIndex[sortRow$ix, sortCol$ix]
                
                # replace upper trinagle with NaN
                if (triangleType=="half"){
                        triangleGeneIndex[upper.tri(triangleGeneIndex)]<-NaN
                }
                triangleGeneIndex <- as.matrix(triangleGeneIndex)
                
                
                # replace letters with data values
                tDataCol <-ncol(triangleData)
                for (j in 1:nrow(triangleData)){
                        triangleGeneIndex <- gsub(paste0(paste0('^', rownames(triangleData)[j]),'$'),triangleData[j,i],triangleGeneIndex)
                }
                
                if (triangleType=="half"){
                        tDataCol <-ncol(triangleData)
                        for (j in 1:nrow(triangleData)){
                                triangleGeneIndex <- gsub(paste0(paste0('^', rownames(avgUpper)[j]),'$'),triangleData[j,i],triangleGeneIndex)
                        }
                }
                
                
                # convert from character to double for heatmap 
                triangleGeneIndex <- apply(triangleGeneIndex, c(1,2), as.numeric)
                
                # generate pdf of heatmap (no clustering)
                pdf(paste0(paste0(paste0(colnames(triangleData)[i], extension), triangleType), '.pdf'),
                    width = 7, height = 7, family="Arial", useDingbats=FALSE)
                heatmap.2(triangleGeneIndex, Rowv=FALSE, Colv=FALSE, trace="none", 
                          dendrogram="none", keysize=1.1, key.title='', key.xlab='',
                          main=colnames(triangleData)[i], col=my_palette, breaks=colors,
                          density.info='none', na.color='gray80')
                dev.off()
                        
        }
}

avgData <- as.matrix(avgData)

# triangle HM of OD data, half
tData <- avgData
rownames(tData) <- geneList_i

# remove outliers if they exist
# in our experiment, tpo3yor1 and yor1tpo3 were experimental outliers
outlier <- c(which(geneList_i=='jg'), 
             which(geneList_i=='gj'))
avgLower[outlier,]<-NaN
avgUpper[outlier,]<-NaN

my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 48)
colors = unique(c(seq(min(tData, na.rm=TRUE),1,length=25),seq(1,max(tData, na.rm=TRUE),length=25)))

# triangle HM of OD data, full
tData <- rbind(avgLower, avgUpper)
triangleHM(tData, "full", colors, "_OD")

# track single gene KO data 
singlesData <- allData[singlesIndex,]
rownames(singlesData) <- singles
singlesAvgData <- as.matrix(avgData[singlesIndex,])
rownames(singlesAvgData) <- singles

# pull out avg data ratio matrix to repopulate with expected growth rates
expectedGrowth <- allData
rownames(expectedGrowth) <- geneList_i

for (r in 1:nrow(expectedGrowth)){
        for (c in 1:ncol(expectedGrowth)){
                strain <- rownames(expectedGrowth)[r]
                splitStrain <- strsplit(strain, "")
                val <- matrix(data=NA,nrow=length(splitStrain[[1]]),ncol=1)
                for (l in 1:length(splitStrain[[1]])){
                        letter <- splitStrain[[1]][l]
                        valIndex <- match(letter, rownames(singlesData))
                        val[l] <- singlesData[valIndex,c]
                }
                expectedVal <- prod(val)
                expectedGrowth[r,c] <- expectedVal

        }
}

# rename rows for expected data matrix
rowsExpectedGrowth <- paste0(paste0('-', rownames(expectedGrowth)), '-')
for (i in 1:nrow(geneNames)){
        rowsExpectedGrowth <- gsub(paste0('^-', geneNames[i,1]),geneNames[i,2],rowsExpectedGrowth)
        rowsExpectedGrowth <- gsub(paste0(geneNames[i,1], '-$'),geneNames[i,2],rowsExpectedGrowth)
}
rownames(expectedGrowth) <- gsub('-', '', rowsExpectedGrowth)

# update strain names for singlesData
rowSD <- paste0(paste0('-', rownames(singlesData)), '-')
rowSDavg <- paste0(paste0('-', rownames(singlesAvgData)), '-')
for (i in 1:nrow(geneNames)){
        rowSD <- gsub(paste0('^-', geneNames[i,1]),geneNames[i,2], rowSD)
        rowSD <- gsub(paste0(geneNames[i,1], '-$'),geneNames[i,2], rowSD)
        rowSDavg <- gsub(paste0('^-', geneNames[i,1]),geneNames[i,2], rowSDavg)
        rowSDavg <- gsub(paste0(geneNames[i,1], '-$'),geneNames[i,2], rowSDavg)
}
rownames(singlesData)<-gsub('-', '', rowSD)
rownames(singlesAvgData)<-gsub('-', '', rowSDavg)

# remove outliers if they exist
# in our experiment, tpo3yor1 and yor1tpo3 were experimental outliers
outlier <- c(which(rownames(allData)=='tpo3yor1'),
             which(rownames(allData)=='yor1tpo3'))
allData <- allData[-outlier,]
avgData <- avgData[-outlier,]
avgUpper <- avgUpper[-outlier,]
avgLower <- avgLower[-outlier,]
expectedGrowth <- expectedGrowth[-outlier,]


## compare expected growth numbers to actual growth values (p-value based)

# initialize matrices 
avgData <- as.matrix(avgData)
avgUpper <- as.matrix(avgUpper)
avgLower <- as.matrix(avgLower)
pVal <- matrix(data=NA,nrow=nrow(allData),ncol=length(completeConditions),
               dimnames=list(rownames(avgData), colnames(avgData)))
pVal_low <- matrix(data=NA,nrow=nrow(allData),ncol=length(completeConditions),
               dimnames=list(rownames(avgData), colnames(avgData)))
pVal_up <- matrix(data=NA,nrow=nrow(allData),ncol=length(completeConditions),
                 dimnames=list(rownames(avgData), colnames(avgData)))
eps <- matrix(data=NA,nrow=nrow(allData),ncol=ncol(allData),
              dimnames=list(rownames(allData), colnames(allData)))
pVal_eps <- matrix(data=NA,nrow=nrow(allData),ncol=length(completeConditions),
                   dimnames=list(rownames(avgData), colnames(avgData)))
epsAvg <- matrix(data=NA,nrow=nrow(allData),ncol=length(completeConditions),
                 dimnames=list(rownames(avgData), colnames(avgData)))
epsUpper <- matrix(data=NA,nrow=nrow(allData),ncol=length(completeConditions),
                   dimnames=list(rownames(avgUpper), colnames(avgUpper)))
epsLower <- matrix(data=NA,nrow=nrow(allData),ncol=length(completeConditions),
                   dimnames=list(rownames(avgLower), colnames(avgLower)))
pVal_epsRecip <- matrix(data=NA,nrow=nrow(allData),ncol=length(completeConditions),
                   dimnames=list(rownames(avgData), colnames(avgData)))

# calculate p-values 
for (i in 1:length(completeConditions)){
        currentCond <- gsub(".csv", "", completeConditions[i])
        currentExpected <- as.matrix(expectedGrowth[,which(colnames(expectedGrowth) %in% currentCond)])
        currentActual <- as.matrix(allData[,which(colnames(allData) %in% currentCond)])
        for (j in 1:nrow(currentExpected)){
                p <- t.test(currentExpected[j,seq(1,ncol(currentExpected),2)],currentActual[j,])
                pLow <- t.test(currentExpected[j,seq(1,ncol(currentExpected),2)],currentActual[j,seq(1,ncol(currentExpected),2)])
                pUp <- t.test(currentExpected[j,seq(1,ncol(currentExpected),2)],currentActual[j,seq(2,ncol(currentExpected),2)])
                pVal[j,i] <- p$p.value
                pVal_low[j,i] <- pLow$p.value
                pVal_up[j,i] <- pUp$p.value
                eps_current <- currentActual[j,]-currentExpected[j,]
                eps[j,(i*noReps-(noReps-1)):(i*noReps)] <- eps_current
                p_eps <- t.test(eps[j,],mu=0)
                pVal_eps[j,i] <- p_eps$p.value
                epsAvg[j,i] <- mean(eps_current)
                epsLower[j,i] <- mean(eps_current[seq(1,ncol(currentExpected),2)])
                epsUpper[j,i] <- mean(eps_current[seq(2,ncol(currentExpected),2)])
                pRecip <- t.test(eps_current[seq(1,ncol(currentExpected),2)], eps_current[seq(2,ncol(currentExpected),2)])
                pVal_epsRecip[j,i] <- pRecip$p.value

        }
}

# eps triangle heatmap "full"
epsAvgT <- rbind(epsLower, epsUpper)
for (i in 1:nrow(geneNames)){
        rownames(epsAvgT) <- gsub(paste0('^', geneNames[i,2]), geneNames[i,1], rownames(epsAvgT))
        rownames(epsAvgT) <- gsub(paste0(geneNames[i,2],'$'), geneNames[i,1], rownames(epsAvgT))
}

colors = unique(seq(-1,1,length=49))

triangleHM(epsAvgT, "full", colors, "_EPS")

pValwSingles <- pVal

# delete single strains from pVal, pVal_eps, epsAvg, eps
rownames(epsLower)<-rownames(eps)
rownames(epsUpper)<-rownames(eps)
pVal <- pVal[-match(rownames(singlesData),rownames(pVal)),]
eps <- eps[-match(rownames(singlesData),rownames(eps)),]
pVal_eps <- pVal_eps[-match(rownames(singlesData),rownames(pVal_eps)),]
epsAvg <- epsAvg[-match(rownames(singlesData),rownames(epsAvg)),]
pVal_low <- pVal_low[-match(rownames(singlesData),rownames(pVal_low)),]
pVal_up <- pVal_up[-match(rownames(singlesData),rownames(pVal_up)),]
pVal_epsRecip <- pVal_epsRecip[-match(rownames(singlesData),rownames(pVal_epsRecip)),]
epsLower <- epsLower[-match(rownames(singlesData),rownames(epsLower)),]
epsUpper <- epsUpper[-match(rownames(singlesData),rownames(epsUpper)),]


# multiple hypothesis correction
pVal <- as.matrix(pVal)
pVal_low <- as.matrix(pVal_low)
pVal_up <- as.matrix(pVal_up)
pVal_epsRecip <- as.matrix((pVal_epsRecip))
pVal_eps <- as.matrix(pVal_eps)

pValAdj <- matrix(data=NA, nrow=nrow(pVal), ncol=ncol(pVal),
                  dimnames=list(rownames(pVal), colnames(pVal)))
pValAdj_low <- matrix(data=NA, nrow=nrow(pVal_low), ncol=ncol(pVal_low),
                  dimnames=list(rownames(pVal_low), colnames(pVal_low)))
pValAdj_up <- matrix(data=NA, nrow=nrow(pVal_up), ncol=ncol(pVal_up),
                  dimnames=list(rownames(pVal_up), colnames(pVal_up)))
pValAdj_recip <- matrix(data=NA, nrow=nrow(pVal_epsRecip), ncol=ncol(pVal_epsRecip),
                        dimnames=list(rownames(pVal_epsRecip), colnames(pVal_epsRecip)))
pValAdj_eps <- matrix(data=NA, nrow=nrow(pVal_eps), ncol=ncol(pVal_eps),
                      dimnames=list(rownames(pVal_eps), colnames(pVal_eps)))

for (i in 1:ncol(pVal)){
        pValAdj[,i] <- p.adjust(pVal[,i], method="BH")
}
for (i in 1:ncol(pVal_low)){
        pValAdj_low[,i] <- p.adjust(pVal_low[,i], method="BH")
}
for (i in 1:ncol(pVal_up)){
        pValAdj_up[,i] <- p.adjust(pVal_up[,i], method="BH")
}
for (i in 1:ncol(pVal_eps)){
        pValAdj_eps[,i] <- p.adjust(pVal_eps[,i], method="BH")
}
for (i in 1:ncol(pVal_epsRecip)){
        pValAdj_recip[,i] <- p.adjust(pVal_epsRecip[,i], method="BH")
}

# Look at some pVal stats for reciprocal knockouts 
bin = as.matrix(pValAdj_recip[,colnames(pValAdj_recip)=="NT"]>0.05)
print(paste("recipNT", sum(bin==TRUE)/(nrow(bin))*100))
print(paste("recipAll", sum(as.matrix(pValAdj_recip>0.05)==TRUE)/(nrow(pValAdj_recip)*ncol(pValAdj_recip))*100))


# postive and negative interaction lists (based on eps)
binary_eps <- matrix(data=0,nrow=nrow(pValAdj_eps),ncol=ncol(pValAdj_eps))
binary_eps[which(pValAdj_eps<0.05 & epsAvg<0)] <- -1 
binary_eps[which(pValAdj_eps<0.05 & epsAvg>0)] <- 1
rownames(binary_eps)<-rownames(pValAdj_eps)
colnames(binary_eps)<-colnames(pValAdj_eps)

# postive and negative interaction lists (based on old p-value calc)
binary <- matrix(data=0,nrow=nrow(pValAdj),ncol=ncol(pValAdj))
binary[which(pValAdj<0.05 & epsAvg<0)] <- -1 
binary[which(pValAdj<0.05 & epsAvg>0)] <- 1
rownames(binary)<-rownames(pValAdj)
colnames(binary)<-colnames(pValAdj)

binary_low <- matrix(data=0,nrow=nrow(pValAdj_low),ncol=ncol(pValAdj_low))
binary_low[which(pValAdj_low<0.05 & epsLower<0)] <- -1 
binary_low[which(pValAdj_low<0.05 & epsLower>0)] <- 1
rownames(binary_low)<-rownames(pValAdj_low)
colnames(binary_low)<-colnames(pValAdj_low)

binary_up <- matrix(data=0,nrow=nrow(pValAdj_up),ncol=ncol(pValAdj_up))
binary_up[which(pValAdj_up<0.05 & epsUpper<0)] <- -1 
binary_up[which(pValAdj_up<0.05 & epsUpper>0)] <- 1
rownames(binary_up)<-rownames(pValAdj_up)
colnames(binary_up)<-colnames(pValAdj_up)

# Look at some pVal stats for reciprocal knockouts 
binary_comp <- (binary_up == binary_low) & (binary_up !=0) & (binary_low !=0)
print(sum(binary_comp==TRUE)/(nrow(binary_comp)*ncol(binary_comp))*100)
print(sum(binary_comp[,colnames(binary_comp)=="NT"]==TRUE)/(nrow(binary_comp))*100)

binary_comp_print <- binary_comp
binary_comp_print[binary_comp_print==TRUE] <- 1
binary_comp_print[binary_comp_print==FALSE] <- 0 

# write FC and pValue data to spreadsheets
write.csv(epsAvg, "epsAvg.csv")
write.csv(epsLower, "eps_low.csv")
write.csv(epsUpper, "eps_up.csv")
write.csv(eps, "eps.csv")
write.csv(pVal_eps, "pVal_eps.csv")
write.csv(pVal, "pVal.csv")
write.csv(pVal_low, "pVal_low.csv")
write.csv(pVal_up, "pVal_up.csv")
write.csv(pValAdj, "pValAdj.csv")
write.csv(pValAdj_low, "pValAdj_low.csv")
write.csv(pValAdj_up, "pValAdj_up.csv")
write.csv(pValAdj_eps, "pValAdj_eps.csv")
write.csv(pValAdj_recip, "pValAdj_recip.csv")
write.csv(binary_eps, "binary_eps.csv")
write.csv(binary, "binary.csv")
write.csv(binary_low, "binary_low.csv")
write.csv(binary_up, "binary_up.csv")
write.csv(expectedGrowth, "expected.csv")
write.csv(binary_comp_print, "binary_comp_print.csv")
write.csv(pVal_epsRecip, "pVal_epsRecip.csv")

my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 48)


# heatmap of average eps values 
colors = unique(seq(-1,1,length=49))
epsAvg <- as.matrix(epsAvg)
pdf(paste0('epsAvg', '.pdf'), width = 5, height = 10, family="Arial", useDingbats=FALSE)
par(mar=c(5, 4, 4, 2))
heatmap.2(epsAvg,scale="none",trace="none",cexRow=0.6,cexCol=0.8,keysize=.3,
          key.title='', col=my_palette, breaks=colors, symkey=FALSE,
          key.xlab='', density.info='none',
          lmat=rbind(rbind(c(0,3),c(2,1)),c(0,4)), lhei=c(2,4,.85), lwid=c(2,4))
dev.off()


# heatmap of average data matrix 
colors = unique(c(seq(min(avgData, na.rm=TRUE),1,length=25),seq(1,max(avgData, na.rm=TRUE),length=25)))
pdf(paste0('avgData', '.pdf'), width = 5, height = 10, family="Arial", useDingbats=FALSE)
par(mar=c(5, 4, 4, 2))
heatmap.2(avgData,scale="none",trace="none",cexRow=0.6,cexCol=0.8,keysize=.3,
          key.title='', col=my_palette, breaks=colors, symkey=FALSE,
          key.xlab='', density.info='none',
          lmat=rbind(rbind(c(0,3),c(2,1)),c(0,4)), lhei=c(2,4,.85), lwid=c(2,4))
dev.off()


