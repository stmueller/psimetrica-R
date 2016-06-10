## This file must be loaded befoer the phonological assessment
## functions can be used, but you need to load the functions in
## phonology.R because it uses makeCompmatrix().
## It sets up several global variable/environments which is
## used a hash/lookup table.
##
### Read in the phonological coding library.
### This should have reasonable phonological codings of
### a lot of english words. (159K)

#dyn.unload("../lib/lev_weighted.so")
dyn.load("../lib/lev_weighted.so")

##
levdist <- function(str1, str2,subst,ins1,ins2,print=0)
{
  # check type
#  if (typeof(str1) != "character" && class(str1) != "factor")
#    stop(sprintf("Illegal data type: %s", typeof(str1)))
#  if (class(str1) == "factor")
#    str=as.character(str1)
#  
#  if (typeof(str2) != "character" && class(str2) != "factor")
#    stop(sprintf("Illegal data type: %s", typeof(str2)))
#  if (class(str2) == "factor")
#    str=as.character(str2)
#  if ((is.array(str1) || is.array(str2)) && dim(str1)!=dim(str2))
#    stop ("non-conformable arrays")
#  if (any(is.na(str1),is.na(str2)))
#    out$ans[is.na(str1) | is.na(str2)]=NA
#  
#  if (is.array(str1))
#    return(array(out$ans,dim(str1)))

#  if (is.array(str2))
#    return(array(out$ans,dim(str2)))
#  return(out$ans)


    if(length(str1)==0 || length(str2)==0) return(integer(0))
    l1 <- length(str1)
    l2 <- length(str2)
    
    comp <- as.vector(t(subst))
    #cat("calling lev_weighted with: ", l1, ",",l2,":", comp,"\n and print: ", print , "\n")
   
    out <- .C("lev_weighted",l1,l2, comp,ins1,ins2,print,ans=0)
    
    return (out$ans)
}

sortme <- function(l,n=50){sort(l)[1:n]}


## This will load the MTU dictionary--a cleaned and updated version of
## CMU dictionary (around 0.6) with additional fixes, words, and
##stress characters removed.  Note that it still uses the CMU dictionary
##phoneme codes
loadMTUDict <- function(filename="../data/mtudict.txt",recache=F)
    {

        phono <<- scan(filename,what=character(),sep="-")
        phonopairs <<- strsplit(phono,"  ")
        PHONOWORDS <<- sapply(phonopairs,function(x){x[[1]]})
        
        if(file.exists(paste(filename,".Rimage",sep=""))&recache)
            {
                
                print("loading precached phono map image file")
                load(paste(filename,".Rimage",sep=""))
                print("Done loading  phono map image file")
                
            }else   {
                
                PHONOMAP <<- new.env()
                
                ##this goes through each line in the lexicon and reads it into
                ##the PHONOMAP environment.
                for(i in 1:length(phonopairs))
                    {
                        pair <- phonopairs[i][[1]]
                        word <- pair[1]
                        phones <- strsplit(pair[2]," ")[[1]]
                        
                        
                        ##Now, add all possible pairs
                        if(!is.null(PHONOMAP[[word]]))
                            {
                                phonelist <- append(PHONOMAP[[word]],list(phones))
                                
                            }else{
                                
                                phonelist <- list(phones)
                                
                            }
                        PHONOMAP[[word]] <- phonelist
                        
                        if(i %% 1000==0)print(paste(i,pair))
                    }
                
                
                save(PHONOMAP,file=paste(filename,".Rimage",sep=""))
            }
        
    }



## This loads a CMU-formatted dictionary; hopefully one taken directly
## from the cmu website.
## This lexicon includes stress markers. These can be auto-trimmed out
## or used directly.  Currently, there is no way to transform all stress
## markers to unstressed and use those.

loadCMUDict <- function(filename="../data/cmudict-0.7b")
    {

        phono <<- scan(filename,what=character(),
                       comment.char=";",sep="\n")
        phonopairs <<- strsplit(phono,"  ")
        PHONOWORDS <<- sapply(phonopairs,function(x){x[[1]]})
        
        if(file.exists(paste(filename,".Rimage",sep="")))
            {
                
                print("loading precached phono map image file")
                load(paste(filename,".Rimage",sep=""))
                print("Done loading  phono map image file")
                
            }else   {
                
                PHONOMAP <<- new.env()
                
                ##this goes through each line in the lexicon and reads it into
                ##the PHONOMAP environment.
                for(i in 1:length(phonopairs))
                    {
                        pair <- phonopairs[i][[1]]
                        word <- pair[1]
                        phones <- strsplit(pair[2]," ")[[1]]
                        
                        
                        ##Now, add all possible pairs
                        if(!is.null(PHONOMAP[[word]]))
                            {
                                phonelist <- append(PHONOMAP[[word]],list(phones))
                                
                            }else{
                                
                                phonelist <- list(phones)
                                
                            }
                        PHONOMAP[[word]] <- phonelist
                        
                        if(i %% 1000==0)print(paste(i,pair))
                    }
                
                
                save(PHONOMAP,file=paste(filename,".Rimage",sep=""))
            }
        
    }


## This sets up a phoneme translation table that transforms
## CMU phone set to the psimetrica phone set.  It can,
## based on the option useStress, incorporate or ignore
## stress markers.

loadPhonemeTranslation <- function(useStress=T)
    {



        stressvowels <- list(list("AA",list("ah","stress0")),
                             list("AA0",list("ah","stress0")),
                             list("AA1",list("ah","stress1")),
                             list("AA2",list("ah","stress2")),

            
                             list("AE",list("ae","stress0")),
                             list("AE0",list("ae","stress0")),
                             list("AE1",list("ae","stress1")),
                             list("AE2",list("ae","stress2")),

                             list("AH",list("uu", "stress0")),
                             list("AH0",list("uu", "stress0")),
                             list("AH1",list("uu", "stress1")),
                             list("AH2",list("uu", "stress2")),
                    
                             
                             list("AO",list("aw","stress0")),
                             list("AO0",list("aw","stress0")),
                             list("AO1",list("aw","stress1")),
                             list("AO2",list("aw","stress2")),
                      
                             list("AW",list("ah","u", "stress0")),
                             list("AW0",list("ah","u", "stress0")),
                             list("AW1",list("ah","u", "stress1")),
                             list("AW2",list("ah","u", "stress2")),
                             
                             list("AY",list("ah", "ee","stress0")),
                             list("AY0",list("ah", "ee","stress0")),
                             list("AY1",list("ah", "ee","stress1")),
                             list("AY2",list("ah", "ee","stress2")),
                             
                             list("EH",list("e", "stress0")),
                             list("EH0",list("e", "stress0")),
                             list("EH1",list("e", "stress1")),
                             list("EH2",list("e", "stress2")),
                             


                             list("ER",list("@","stress0","r")),
                             list("ER0",list("@","stress0","r")),
                             list("ER1",list("@","stress1","r")),
                             list("ER2",list("@","stress2","r")),
                             
                             
                             list("EY",list("a_","stress0")),
                             list("EY0",list("a_","stress0")),
                             list("EY1",list("a_","stress1")),
                             list("EY2",list("a_","stress2")),

                             
                             list("IH",list("i","stress0")),
                             list("IH0",list("i","stress1")),
                             list("IH1",list("i","stress2")),
                             list("IH2",list("i","stress3")),
                             
                             list("IY",list("ee","stress0")),
                             list("IY0",list("ee","stress1")),
                             list("IY1",list("ee","stress2")),
                             list("IY2",list("ee","stress3")),

                             list("OW",list("o_", "stress0")),
                             list("OW0",list("o_", "stress0")),
                             list("OW1",list("o_", "stress1")),
                             list("OW2",list("o_", "stress2")),
                           
                             list("OY",list("o_","ee", "stress0")),
                             list("OY0",list("o_","ee","stress0")),
                             list("OY1",list("o_","ee","stress1")),
                             list("OY2",list("o_","ee","stress2")),
                      
                             
                             list("UH",list("oo","stress0")),
                             list("UH0",list("oo","stress0")),
                             list("UH1",list("oo","stress1")),
                             list("UH2",list("oo","stress2")),
                      
                           
                             list("UW",list("u_","stress0")),
                             list("UW0",list("u_","stress0")),
                             list("UW1",list("u_","stress1")),
                             list("UW2",list("u_","stress2"))
                           )

    
        unstressvowels <- list(
            list("AA",list("ah")),
            list("AA0",list("ah")),
            list("AA1",list("ah")),
            list("AA2",list("ah")),
            
            
            list("AE",list("ae")),
            list("AE0",list("ae")),
            list("AE1",list("ae")),
            list("AE2",list("ae")),
            
            list("AH",list("uu")),
            list("AH0",list("uu")),
            list("AH1",list("uu")),
            list("AH2",list("uu")),
                             
            list("AO",list("aw")),
            list("AO0",list("aw")),
            list("AO1",list("aw")),
            list("AO2",list("aw")),
            
            list("AW",list("ah","u")),
            list("AW0",list("ah","u")),
            list("AW1",list("ah","u")),
            list("AW2",list("ah","u")),
                             
            list("AY",list("ah", "ee")),
            list("AY0",list("ah", "ee")),
            list("AY1",list("ah", "ee")),
            list("AY2",list("ah", "ee")),
                             
            list("EH",list("e")),
            list("EH0",list("e")),
            list("EH1",list("e")),
            list("EH2",list("e")),
            
            list("ER",list("@","r")),
            list("ER0",list("@","r")),
            list("ER1",list("@","r")),
            list("ER2",list("@","r")),
            
            
            list("EY",list("a_")),
            list("EY0",list("a_")),
            list("EY1",list("a_")),
            list("EY2",list("a_")),
            
            list("IH",list("i")),
            list("IH0",list("i")),
            list("IH1",list("i")),
            list("IH2",list("i")),
                             
            list("IY",list("ee")),
            list("IY0",list("ee")),
            list("IY1",list("ee")),
            list("IY2",list("ee")),
            
            list("OW",list("o_")),
            list("OW0",list("o_")),
            list("OW1",list("o_")),
            list("OW2",list("o_")),
            
            list("OY",list("o_","ee")),
            list("OY0",list("o_","ee")),
            list("OY1",list("o_","ee")),
            list("OY2",list("o_","ee")),
                      
                             
            list("UH",list("oo")),
            list("UH0",list("oo")),
            list("UH1",list("oo")),
            list("UH2",list("oo")),
                      
            
            list("UW",list("u_")),
            list("UW0",list("u_")),
            list("UW1",list("u_")),
            list("UW2",list("u_")),

            ##some re-translations if they used the PEBL format.
            list("@",  list("@")),
            list("A_", list("a_")),
            list("I", list("i")),
            list("EE", list("ee")),
            list("OO", list("oo")),
            list("O_", list("o_"))
            
        )

        

                           
        consonants <- list(  list("B",list("b")),
                           list("CH",list("ch")),
                           list("D",list("d")),
                           list("DH",list("thz")),
                           list("F",list("f")),
                           list("G",list("g")),
                           list("HH",list("h")),
                           list("H",list("h")),
                           
                           list("JH",list("j")),
                           list("K",list("k")),
                           list("L",list("l")),
                           list("M",list("m")),
                           list("N",list("n")),
                           list("NG",list("ng")),
                           list("P",list("p")),
                           list("R",list("r")),
                           list("S",list("s")),
                           list("SH",list("sh")),
                           list("T",list("t")),
                           list("TH",list("th")),
                           
                           list("V",list("v")),
                           list("W",list("w")),
                           list("Y",list("y")),
                           list("Z",list("z")),
                           list("ZH",list("zh"))
                           
                           )
        
        if(useStress)
            {
                trans <- append(stressvowels,consonants)
            }else{
                trans <- append(unstressvowels,consonants)
            }

        PHON2PHON <<- new.env()

        for(i in 1:length(trans))
            {
                PHON2PHON[[ trans[[i]][[1]] ]] <- trans[[i]][[2]]
            }


    }

## This will create several globals that are used by other functions to
## assess phonological similarity.  These include:
## FEATURES; the feature matrix
## SYMBOLS: set of feature symbols
## COMPCACHE: matrix indexed by SYMBOLS recording pairwise dissimilarity
## values.
##
## FEATURES is only used by makeCompMatrix(); if an alternative similarity
## matrix is desired, be sure to create SYMBOLS and COMPCACHE, but FEATURES
## would not be needed.

loadFeatureSimilarity <- function(type="psim")
    {

        ##other options could be implemented, including
        ##phonetic approach:
        ##        http://aix1.uottawa.ca/~jmielke/research/NELS_similarity_Mielke.pdf

        ##see also Frisch 1996, etc.

        if(type=="psim")
            {
                loadFeatureSimilarity.psim()
            }
        else if(type=="simple")
            {
                loadFeatureSimilarity.vowelconsonant()
            }
        else if(type=="levenshein")
            {
                loadFeatureSimilarity.levenshtein()
            }
    }



loadFeatureSimilarity.psim <- function(filename="../data/phonemes.csv",
                                  inscost = 1.0,
                                  featureweights,
                                  stressweights)
    {
        FEATURES <<- read.csv(filename,header=F)
        ##make sure factor levels are in the same order as the feature file:
        tmp <- c(as.character(FEATURES[,1]),"stress0","stress1","stress2")
        SYMBOLS <<- factor(tmp,levels=tmp)

        ##create a cached feature comparison matrix
        ##this takes less than 10 seconds.
        print("Computing phoneme-to-phoneme dissimilarity matrix.")
        COMPCACHE <<- makeCompMatrix(as.vector(SYMBOLS),as.vector(SYMBOLS),
                                     featureweights,stressweights,printstatus=T)
        print("Done computing phoneme-to-phoneme dissimilarity matrix.")
    }



## This creates a very simple phoneme similarity matrix: 
## identical phonemes have dissimilarity of 0
## within-class (consonant-cons or vowel-vowel) have dissimilarity of 0.25
## partial-class overlap (glide-vowel or glide-consonant) have dissimlarity of 0.5
## cross-class (cons-vowel) have dissimilarity of 1.0

loadFeatureSimilarity.vowelconsonant <- function(withindiss = 0.25,
                                                 glidediss = .5,
                                                 crossdiss = 1.0,
                                                 stressdiss = 2.0)

    {

        print("Computing phoneme-to-phoneme dissimilarity matrix.")

        
        vowels <- c("ah","ae","uu","aw", "a_", "e", "@", 
                    "i","ee","o_","oo","u_","u")
               
        glides <- c("w","y","l","r")
        
        
        consonants <- c("b","ch","d","thz","f","g","h",
                        "j","k","m","n","ng","p","s","sh",
                        "t","th","v","z","zh")
        
        stress <- c("stress0","stress1","stress2")                                   
        

               
        tmp <- c(vowels,glides,consonants,stress)
        SYMBOLS <<- factor(tmp,levels=tmp)
        
        classes <- c(rep(1,length(vowels)),
                     rep(2,length(glides)),
                     rep(3,length(consonants)),
                     rep(4,length(stress)))
        
        COMPCACHE <- matrix(crossdiss,nrow=length(classes),ncol=length(classes))

        ##glide dissimilarities:
        COMPCACHE[classes==2,] <- glidediss
        COMPCACHE[,classes==2] <- glidediss


        ##stress dissimilairty
        COMPCACHE[classes==4,] <- stressdiss
        COMPCACHE[,classes==4] <- stressdiss

        ##withinclass dissimilarity
        COMPCACHE[outer(classes, classes, "==")] <- withindiss

        
        diag(COMPCACHE) <- 0
        COMPCACHE <<- COMPCACHE  ##make it global
        
        ##create a cached feature comparison matrix
        ##this takes less than 10 seconds.



        
        print("Done computing phoneme-to-phoneme dissimilarity matrix.")
    }


loadFeatureSimilarity.levenshtein <- function(filename="../data/phonemes.csv",
                                  mismatchcost = 1.0)
{

        FEATURES <<- read.csv("../data/phonemes.csv",header=F)
        ##make sure factor levels are in the same order as the feature file:
        tmp <- c(as.character(FEATURES[,1]),"stress0","stress1","stress2")
        SYMBOLS <<- factor(tmp,levels=tmp)

        ##create a cached feature comparison matrix
        ##this takes less than 10 seconds.
        print("Computing phoneme-to-phoneme dissimilarity matrix.")
        cc <- matrix(mismatchcost,nrow=length(SYMBOLS),ncol=length(SYMBOLS))
        diag(cc)<-0
        
        COMPCACHE <<- cc

}

getWord <- function(word)
    {
        words <- PHONOMAP[[ (word)]]
        if(is.null(words))
            {
                warning(paste("Word [",word,"] does not exist"))
                words <- list(NA)
            }
        if(length(words)>1)
            {
                print("more than one pronunciation found")
            }
        return(words)
    }


translateWord <- function(pword)
{
    if(!is.vector(pword))
        {
            print(pword)
            if(is.na(pword))
                return(NA)
        }

    
    newword <- c()
    for(i in pword)
        {
            
            newphones <- unlist(PHON2PHON[[i]])
            print(newphones)
            if(is.null(newphones))
            {
                print(paste("untranslateable phoneme",i,"in word " ,
                            paste(pword,collapse="-")))
            }
            
            newword <- c(newword, newphones)

        }
 return (newword)
}


## This makes a comparison matrix, 
##weighing each feature equally.
##also, it ignores vowel/syllable stress
##pword1 and pword2 need to be vectors containing
## have character strings matching the first column
## of FEATURES.

##featureweights should be a weight matrix summing to 1.0 to weight the fatures
##stressweights should be a 4x4 matrix showing the relative distance

##Insertion/deletion of phoneme cost is recorded in a vector.
##stress markers 
makeInsertCost <- function(pword,inscost,stresscost)
    {

        cost <- rep(inscost,length(pword))
        cost[pword=="stress0"] <- stresscost
        cost[pword=="stress1"] <- stresscost
        cost[pword=="stress2"] <- stresscost

        return(cost)
    }


makeCompMatrix <- function(pword1,pword2,
                           featureweights,
                           stressweights,printstatus=F)
{

#    print(pword1)
#    print(pword2)
    ##right now, this is hard-coded to use the FEATURES
    ##matrix, which has 14 columns (the identifier + 13 binary features.

    if(missing(featureweights))
        {
            featureweights <- rep(1,13)/13
        }

    if(missing(stressweights))
        {
            stressweights <- rbind(c(0,   2,  2, 2),  ##first row is comparison of stress to non-stress--should be large.
                                   c(2,   0, .5, 1.0),  ##second row is comparison of stress0 to others
                                   c(2,  .5,  0, .5),   ##third row is comparison of stress1 to others
                                   c(2, 1.0, .5, 0)   ##fourth row is comparison of stress2 to others

                                   ##fifth row/column is the comparison to anything but a stress symbol.
                                   )
        }
    mat <- matrix(NA, nrow=length(pword1),ncol=length(pword2))

    if(printstatus)
        {
            cat(paste(rep("-",length(pword1)),collapse=""),"\n")
        }
    for(i in 1:length(pword1))
        {

            cat(".")

            #Find the row of FEATURES which match SYMBOLS.
            frow <- which(SYMBOLS==pword1[i])

            for(j in 1:length(pword2))
                {


                    fcol <- which(SYMBOLS==pword2[j])

                    
                    if(fcol>nrow(FEATURES)  || ( frow>nrow(FEATURES)))
                        {
                            stressvali <- which(pword1[i]==c("stress0","stress1","stress2"))
                            if(length(stressvali)==0)
                                {
                                    stressvali <- 0
                                }


                            
                            stressvalj <- which(pword2[j]==c("stress0","stress1","stress2"))
                            if(length(stressvalj)==0)
                                {
                                    stressvalj <- 0
                                }
                            
                            mat[i,j] <-stressweights[stressvali+1,stressvalj+1]
                        
                        }else{
                            ##normal feature-feature comparison:

                            f1 <- FEATURES[frow,2:14]
                            f2 <- FEATURES[fcol,2:14]
                            ##get rid of any where both are null (x)
                            filter <- !(f1=="x" & f1==f2)
                            match <- (f1==f2)[filter]
                            mat[i,j] <- 1-weighted.mean(match,featureweights[filter])
                        }
                }
        }

    cat("\n")

    return (mat)
}

##this will extract comp matrix from COMPCACHE; much faster.
getCompMatrix <- function( pword1, pword2)
{
    f1 <- factor(pword1,levels=SYMBOLS)
    f2 <- factor(pword2,levels=SYMBOLS)
    COMPCACHE[f1,f2]
}

##internal getDistance_ function; requires recoded words.
getDistance_ <- function(pword1,pword2,pmat,ins1,ins2,print=F)
    {
        if(print)
            {
                print(ins1)
                print(ins2)
            }
        
        ##pword1 is the left (rows)
        ##pword2 is the right (columns)
        planargraph <- matrix(NA, nrow=length(pword1)+1,ncol=length(pword2)+1)
        stepgraph <- matrix(NA, nrow=length(pword1)+1,ncol=length(pword2)+1)
        
        if(length(ins1)==1)
            {
                planargraph[,1] <- 0:(length(pword1))*ins1

            }else{
                planargraph[,1] <- c(0,cumsum(ins1))

            }
        stepgraph[,1] <- 0:(length(pword1))

        if(length(ins2)==1)
            {
                planargraph[1,] <- 0:(length(pword2))*ins2
            }else{
                planargraph[1,] <- c(0,cumsum(ins2))
            }
            
        stepgraph[1,] <- 0:(length(pword2))

        
        for(i in 2:(length(pword1)+1))
            {
                for(j in 2:(length(pword2)+1))
                    {
                        ##These are the transition costs
                        up <- planargraph[i-1,j]+ins1[i-1]
                        left <- planargraph[i,j-1]+ins2[j-1]
                        uleft <- planargraph[i-1,j-1] + pmat[i-1,j-1]

                        ##for serious debug-fu:
                        ##cat("------------------\n")
                        ##cat("using: ", i,j, "row ",i-1, "col",j-1,"\n")
                        
                        ##cat(ins2, "for value: " , j-1 , "\n")
                        ##cat("LEFT:",planargraph[i,j-1], " + " ,ins2[j-1] ,"=",left, "\n")
                        ##cat("UP:",planargraph[i-1,j], " + " ,ins1[i-1] ,"=",up, "\n")
                        ##cat("UL:",planargraph[i-1,j-1], " + " ,pmat[i-1,j-1] ,"=",uleft , "\n")
                        ##print(paste(i,j,up,left,uleft))


                        ##the best cost is the cheapest + current mismatch
                        whichisit <- which.min(c(up,left,uleft))
                        ups <- c(-1,0,-1)
                        lefts <- c(0,-1,-1)
                        
                        stepgraph[i,j] <- stepgraph[i+ups[whichisit],j+lefts[whichisit]]+1
                        planargraph[i,j] <- min(up,left,uleft) 

                    }
            }
        if(print)
            {
                print(pmat)
                print(stepgraph)
                print(round(planargraph,3))
            }

        return(list(planargraph[length(pword1)+1,length(pword2)+1],
                    stepgraph[length(pword1)+1,length(pword2)+1]))
                    
    }


getDistance <- function(word1,word2,pmat=NULL,print=F)
    {
        
        pword1  <- translateWord(getWord(word1)[[1]])
        pword2  <- translateWord(getWord(word2)[[1]])

        
        if(print)
            {
                print(pword1)
                print(pword2)
            }
        if(is.null(pmat))
            pmat <- getCompMatrix(pword1,pword2)


        ins1 <-makeInsertCost(pword1,1.0,2.0)
        ins2 <-makeInsertCost(pword2,1.0,2.0)

        getDistance_(pword1,pword2,pmat,ins1,ins2,print)        

    }




## This uses a compiled version of partial levenshtein distance.
getDistance2 <- function(word1,word2,inscost=1.0,stressinsertcost=2.0,pmat=NULL,print=F)
    {
        pword1  <- translateWord(getWord(word1)[[1]])
        pword2  <- translateWord(getWord(word2)[[1]])
        printInt = as.integer(print)
            
        
        if(!(is.list(pword1) && is.list(pword2)))
            {
                if(is.na(pword1) || is.na(pword2))
                    {
                        return(NA)
                    }
            }

        if(print)
            {
                print(pword1)
                print(pword2)

            }
        
        if(is.null(pmat))
            {
                pmat <- getCompMatrix(pword1,pword2)
            }

        ins1 <-makeInsertCost(pword1,inscost,stressinsertcost)
        ins2 <-makeInsertCost(pword2,inscost,stressinsertcost)
        levdist(pword1,pword2,pmat,ins1,ins2,printInt)        
    }





getDistanceDistribution <- function(word,inscost=1.0,stressinsertcost=2.0,
                                    confirm=F)
    {

        if(!confirm)
            {
                stop("getDistanceDistribution() can use substantial resources\nto calculate phonological similarity between a probe and all\nwords in the lexicon (e.g., 2-10 minutes).  Please run with\nargument 'confirm=T' to really execute this command.")
            }
        
        print <- 0

        pword1  <- translateWord(getWord(word)[[1]])
        ins1 <-makeInsertCost(pword1,inscost,stressinsertcost)
        out <- rep(NA,length(PHONOWORDS))
        for(i in 1:length(PHONOWORDS))
            {
                pword2  <- translateWord(getWord(PHONOWORDS[i])[[1]])
                pmat <- getCompMatrix(pword1,pword2)

                ins2 <-makeInsertCost(pword2,inscost,stressinsertcost)
  
                out[i] <- levdist(pword1,pword2,pmat,ins1,ins2,print)        
                if(i%%1000==0){cat(i,PHONOWORDS[i],"  ",pword2,"\n")}
            }

        names(out) <- PHONOWORDS
        return(out)
       
    }

