
source("phonology.R")


if(F)
    {

        loadMTUDict(recache=T)
        loadPhonemeTranslation(useStress=F) ##use
#        loadFeatureSimilarity.levenshtein(.5)
        loadFeatureSimilarity(type="psim")
        save.image()
    }


vw <- read.csv("../data/vitzwinkler.csv")
vw$psimR2 <- NA
vw$psimR <- NA
vw$length <- NA
vw$maxlength <- NA

##vw contains the transcribed words already; translated into
## the CMU format.

inscost <- .5
stressinsertcost <- 2.0

print <- TRUE
for(i in 1:nrow(vw))
    {

        word1 <-  strsplit(as.character(vw$standardPH[i])," ")[[1]]
        word2 <- strsplit(as.character(vw$X[i])," ")[[1]]
        pword1 <- translateWord(word1)
        pword2 <- translateWord(word2)
        pmat <- getCompMatrix(pword1,pword2)
        ins1 <-makeInsertCost(pword1,inscost,stressinsertcost)
        ins2 <-makeInsertCost(pword2,inscost,stressinsertcost)
        #ldist <- levdist(pword1,pword2,pmat,ins1,ins2,print) 
        dist2 <- getDistance_(pword1,pword2,pmat,ins1,ins2,print)
        
        vw$psimR[i] <- dist2[[1]]
        vw$length[i] <- dist2[[2]]
        vw$maxlength[i] <- max(length(pword1),length(pword2))
    }

vw$normalized <- vw$psimR2/vw$length
vw$normalized2 <- vw$psimR/vw$length
cor(vw[,c(7,8,9,10,13,14)])

plot(vw$data,vw$vwpred)
cor(vw$data,vw$vwpred)

plot(vw$data,vw$psimR)
cor(vw$data,vw$psimR)

plot(vw$data,vw$psimR2)
cor(vw$data,vw$psimR2)
