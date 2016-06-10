#source("lexicalsim.R")
source("phonology.R")
#source("setup-phon.R")

if(T)
    {
        loadMTUDict(recache=F)
        loadPhonemeTranslation(useStress=T) ##use
        loadFeatureSimilarity(type="simple") ##basic vowel-consonant differences
##        loadFeatureSimilarity(type="psim")   ##full feature differences
##        loadFeatureSimilarity(type="levenshtein") ##exact matches required

#        ##custom phoneme distance matrix with distinct weightings:
#        loadFeatureSimilarity.psim(inscost = 1.0,c(1,1,1,0,0,0,0,0,0,2,2,2,1))

        save.image()
    }
##



##look at the features
FEATURES


##words are in cmudict or mtudict


dist <- getDistance("horse","abet",print=T)
##first element is total absolute difference
##second is total number of comparisons
##first/second is almost exactly Vitz&Winkler
dist[[1]]/dist[[2]]


dist <- getDistance("horse","force",print=F)
dist
dist[[1]]/dist[[2]]

dist <- getDistance("horse","course")
dist
dist[[1]]/dist[[2]]




dist <- getDistance("horse","choice",print=T)
dist
dist[[1]]/dist[[2]]

dist <- getDistance("horse","choice")
dist
dist[[1]]/dist[[2]]


