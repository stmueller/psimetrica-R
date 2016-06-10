
source("phonology.R")


if(T)
    {

       # loadMTUDict(recache=T)
        loadMTUDict("../data/abhinavdict.txt",recache=F)
        loadPhonemeTranslation(useStress=F) ##use
        loadFeatureSimilarity.levenshtein(.5)
#        loadFeatureSimilarity(type="psim")
        save.image()
    }


ab <- read.csv("../data/abhinav.csv")
ab$length <- NA
simmat <- matrix(0,nrow(ab),nrow(ab))
maxlength <- matrix(0,nrow(ab),nrow(ab))
##ab contains the transcribed words already; translated into
## the CMU format.

inscost <- .5
stressinsertcost <- 2.0

#print <- TRUE
for(i in 1:nrow(ab))
{
    word1 <-  strsplit(as.character(ab$standardPH[i])," ")[[1]]
        plot(1,1,type="n",ylim=c(0,50))
        text(1,1,paste(i,ab$standardPH[i],cex=3))
        print(paste(i,ab$standardPH[i]))
    pword1 <- translateWord(word1)
    
       for(j in 1:nrow(ab))
        {

            word2 <- strsplit(as.character(ab$standardPH[j])," ")[[1]]
            text(1,50-j,paste(j,ab$standardPH[j],cex=3))
            print(paste(j,ab$standardPH[j]))


            pword2 <- translateWord(word2)
            ab$length[i] <- length(pword1)
            ab$length[j] <- length(pword2)
    
            pmat <- getCompMatrix(pword1,pword2)
            ins1 <-makeInsertCost(pword1,inscost,stressinsertcost)
            ins2 <-makeInsertCost(pword2,inscost,stressinsertcost)
            ##ldist <- levdist(pword1,pword2,pmat,ins1,ins2,print) 
            dist2 <- getDistance_(pword1,pword2,pmat,ins1,ins2,print)
            
            simmat[i,j] <- dist2[[1]]
            maxlength[i,j] <-  max(length(pword1),length(pword2))

    }

}

simmat.normalized <- simmat/maxlength


write.csv(simmat,"simmat.csv")
write.csv(simmat.normalized,"simmat-normed.csv")

sortme(getDistanceDistribution("kiyA",confirm=T ))
