#source("lexicalsim.R")
source("phonology.R")
#source("setup-phon.R")

if(F)
    {
#        loadMTUDict("mtushort.txt")
        loadMTUDict(recache=T)
        loadPhonemeTranslation(useStress=T) ##use
        loadFeatureSimilarity(type="simple")
        save.image()
    }
##


getDistance("horse","abet",print=T)
getDistance("horse","abet",print=F)

getDistance2("horse","abet",print=T)
getDistance2("horse","abet",print=F)


getDistance("horse","force",print=T)
getDistance2("horse","force")#compiled version


getDistance("horse","course")
getDistance2("horse","course")



getDistance("horse","choice",print=T)
getDistance2("horse","choice")



getDistance("horse","coltheart")
getDistance2("horse","coltheart")



##This takes way long time:

h <- getDistanceDistribution("horse")

h <- getDistanceDistribution("horse",confirm=T)
sortme(h)    

sortme(getDistanceDistribution("bottle"))
sortme(getDistanceDistribution("insight"))
sortme(getDistanceDistribution("lemon"))
sortme(getDistanceDistribution("coltheart"))

sortme(getDistanceDistribution("pheasant"))



##Here is

ricewords <- read.csv("rice2.csv",header=F)
ricesim <- matrix(NA, nrow(ricewords),nrow(ricewords))
for(i in 1:nrow(ricewords))
    {

        for(j in (i+1):nrow(ricewords))
            {
                sim <- getDistance(ricewords[i,],ricewords[j,])
                ricesim[i,j] <- sim
                ricesim[j,i] <- sim
                print(paste(i,j))
            }
    }
diag(ricesim) <- 0
colnames(ricesim) <- as.vector(unlist(t(ricewords)))
rownames(ricesim) <- as.vector(unlist(t(ricewords)))
write.csv(file="ricesim.csv",ricesim)


sort(ricesim[1,])



##
words <- read.csv("phon_relatedness.csv")
words$sim <- NA
for(i in 1:nrow(words))
    {
        words$sim[i] <- getDistance2(words$Word1[i],words$Word2[i])
    }


write.csv(file="relatedness.csv",words)
