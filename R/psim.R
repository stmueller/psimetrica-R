##This is not the original version of PSIMETRICA, but
##rather a simple edit distance with differential costs.


dyn.load("../src/lev_weighted.so")
##
levdist <- function(str1, str2,subst,inscost)
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
    ##cat("calling lev_weighted with: ", l1, ",",l2,":", inscost,"\n")
    out <- .C("lev_weighted",l1,l2,inscost,
              comp,ans=0)
    
    return (out$ans)
}




### Read in the phonological coding library.
### This should have reasonable phonological codings of
### a lot of english words. (159K)

phono <- scan("mtudict.txt",what=character(),sep="-")
phonopairs <- strsplit(phono,"  ")

PHONOWORDS <- sapply(phonopairs,function(x){x[[1]]})
PHONOMAP <- new.env()



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

     if(i %% 100==0)print(paste(i,pair))
   }


PHON2PHON <- new.env()

trans <- list(list("AA",list("ah")),
               list("AE",list("ae")),
               list("AH",list("uu")),
               list("AO",list("aw")),
               list("AW",list("ah","u")),
               list("AY",list("ah", "ee")),
               list("B",list("b")),
               list("CH",list("ch")),
               list("D",list("d")),
               list("DH",list("thz")),
               list("EH",list("e")),
               list("ER",list("@","r")),
               list("EY",list("a_")),
               list("F",list("f")),
               list("G",list("g")),
               list("HH",list("h")),
               list("IH",list("i")),
               list("IY",list("ee")),
               list("JH",list("j")),
               list("K",list("k")),
               list("L",list("l")),
               list("M",list("m")),
               list("N",list("n")),
               list("NG",list("ng")),
               list("OW",list("o_")),
               list("OY",list("o_","ee")),
               list("P",list("p")),
               list("R",list("r")),
               list("S",list("s")),
               list("SH",list("sh")),
               list("T",list("t")),
               list("TH",list("th")),
               list("UH",list("oo")),
               list("UW",list("u_")),
               list("V",list("v")),
               list("W",list("w")),
               list("Y",list("y")),
               list("Z",list("z")),
               list("ZH",list("zh")))

               

for(i in 1:length(trans))
    {
        PHON2PHON[[ trans[[i]][[1]] ]] <- trans[[i]][[2]]
    }


features <- read.csv("phonemes.csv",header=F)
##create a cached feature comparison matrix
##this takes less than 10 seconds.
COMPCACHE <<- makeCompMatrix(as.list(as.character(features[,1])),
               as.list(as.character(features[,1])))

                     ## Phoneme Example Translation
#        ------- ------- -----------
#        AA	odd     AA D ah
##      AE	at	AE T
#        AH	hut	HH AH T
#        AO	ought	AO T
#        AW	cow	K AW
#        AY	hide	HH AY D
#        B 	be	B IY
#        CH	cheese	CH IY Z
#        D 	dee	D IY
#        DH	thee	DH IY
#        EH	Ed	EH D
#        ER	hurt	HH ER T
#        EY	ate	EY T
#        F 	fee	F IY
#        G 	green	G R IY N
#        HH	he	HH IY
#        IH	it	IH T
#        IY	eat	IY T
#        JH	gee	JH IY
#        K 	key	K IY
#        L 	lee	L IY
#        M 	me	M IY
#       N 	knee	N IY
#        NG	ping	P IH NG
#        OW	oat	OW T
#        OY	toy	T OY
#        P 	pee	P IY
#        R 	read	R IY D
#        S 	sea	S IY
#        SH	she	SH IY
#        T 	tea	T IY
#        TH	theta	TH EY T AH
#        UH	hood	HH UH D
#        UW	two	T UW
#        V 	vee	V IY
#        W 	we	W IY
#        Y 	yield	Y IY L D
#        Z 	zee	Z IY
#        ZH	seizure	S IY ZH ER

##Now read in the phonological symbols and their feature codings.

## ##phonemes <- c(((i1  (+ - + - - - - - + x x x x))
##                         (i2  (+ - + + - - - - + x x x x))
##                         (u1  (+ - + + - - - + + x x x x))
##                         (e1  (+ - - - - - - - + x x x x))
##                         (o1  (+ - - + - - - + + x x x x))
##                         (ae1 (+ - - - + - - - + x x x x))
##                         (a+  (+ - - + + - - - + x x x x))
##                         (-ee (+ - - - + - - + + x x x x))
##                         (-e  (+ - - + + - - + + x x x x))
##                         (i3  (+ - + - - - - - - x x x x))
##                         (u2  (+ - + + - - - + - x x x x))
##                         (e2  (+ - - - - - - - - x x x x))
##                         (u3  (+ - - + - - - - - x x x x))
##                         (o2  (+ - - + - - - + - x x x x))
##                         (ae2 (+ - - - + - - - - x x x x))
##                         (-e2 (+ - - + + - - + - x x x x))
##                         (y   (- - + - - - - - - x x x x))
##                         (w   (- - + + - - - + - x x x x))
##                         (eps (- - - - - - - - - x x x x))
##                         (r   (+ + - - - - + x x + + - -))
##                         (l   (+ + - - - + + x x + + - -))
##                         (p   (- + - - - + - x x - - - -))
##                         (b   (- + - - - + - x x + - - -))
##                         (f   (- + - - - + - x x - + - +))
##                         (v   (- + - - - + - x x + + - +))
##                         (m   (- + - - - + - x x + - + -))
##                         (t   (- + - - - + + x x - - - -))
##                         (d   (- + - - - + + x x + - - -))
##                         (th  (- + - - - + + x x - + - -))
##                         (thz (- + - - - + + x x + + - -))
##                         (n   (- + - - - + + x x + - + -))
##                         (s   (- + - - - + + x x - + - +))
##                         (z   (- + - - - + + x x + + - +))
##                         (c   (- + - - - + + x x - - - +))
##                         (cc  (- + + - - - + x x - - - +))
##                         (jj  (- + + - - - + x x + - - +))
##                         (ss  (- + + - - - + x x - + - +))
##                         (zz  (- + + - - - + x x + + - +))
##                         (k   (- + + + - - - - x - - - -))
##                         (g   (- + + + - - - - x + - - -))
##                         (x   (- + + + - - - - x - + - -))
##                         (ng  (- + + + - - - x x + - + -))
##                         (h   (- - - - + - - x x - + - -))
##                         (kw  (- + + + - - - + x - - - -))
##                         (gw  (- + + + - - - + x + - - -))
##                         (xw  (- + + + - - - + x - + - -))


##                         ))


## ;;;==========================================================================
## ;;;        *phonemes* contains the sub-symbolic description of phonemes
## ;;;==========================================================================
## ;;;This data is derived from chart on pages 176-177 of
## ;;;_The Sound Pattern of English_, Noam Chomsky, 1968.
## ;;;The first 2-3 letter code is the phoneme label, and the rest
## ;;;of the vector represents the speech features that characterize
## ;;;that phoneme.  They are, in order:  vocalic, consonantal, high, back,low,
## ;;;anterior, coronal, round, tense, voice, continuant, nasal, strident.
## ;;;As Chomsky notes, anterior and coronal replace labial, dental, paleto-alveolar, and velar
## ;;;The Xs represent unspecified values.  For the moment, they remain unspecified,
## ;;;but could be determined given phonologic rules and specific words.

## (defvar *phonemes-nc* '((i1  (+ - + - - - - - + x x x x))
##                         (i2  (+ - + + - - - - + x x x x))
##                         (u1  (+ - + + - - - + + x x x x))
##                         (e1  (+ - - - - - - - + x x x x))
##                         (o1  (+ - - + - - - + + x x x x))
##                         (ae1 (+ - - - + - - - + x x x x))
##                         (a+  (+ - - + + - - - + x x x x))
##                         (-ee (+ - - - + - - + + x x x x))
##                         (-e  (+ - - + + - - + + x x x x))
##                         (i3  (+ - + - - - - - - x x x x))
##                         (u2  (+ - + + - - - + - x x x x))
##                         (e2  (+ - - - - - - - - x x x x))
##                         (u3  (+ - - + - - - - - x x x x))
##                         (o2  (+ - - + - - - + - x x x x))
##                         (ae2 (+ - - - + - - - - x x x x))
##                         (-e2 (+ - - + + - - + - x x x x))
##                         (y   (- - + - - - - - - x x x x))
##                         (w   (- - + + - - - + - x x x x))
##                         (eps (- - - - - - - - - x x x x))
##                         (r   (+ + - - - - + x x + + - -))
##                         (l   (+ + - - - + + x x + + - -))
##                         (p   (- + - - - + - x x - - - -))
##                         (b   (- + - - - + - x x + - - -))
##                         (f   (- + - - - + - x x - + - +))
##                         (v   (- + - - - + - x x + + - +))
##                         (m   (- + - - - + - x x + - + -))
##                         (t   (- + - - - + + x x - - - -))
##                         (d   (- + - - - + + x x + - - -))
##                         (th  (- + - - - + + x x - + - -))
##                         (thz (- + - - - + + x x + + - -))
##                         (n   (- + - - - + + x x + - + -))
##                         (s   (- + - - - + + x x - + - +))
##                         (z   (- + - - - + + x x + + - +))
##                         (c   (- + - - - + + x x - - - +))
##                         (cc  (- + + - - - + x x - - - +))
##                         (jj  (- + + - - - + x x + - - +))
##                         (ss  (- + + - - - + x x - + - +))
##                         (zz  (- + + - - - + x x + + - +))
##                         (k   (- + + + - - - - x - - - -))
##                         (g   (- + + + - - - - x + - - -))
##                         (x   (- + + + - - - - x - + - -))
##                         (ng  (- + + + - - - x x + - + -))
##                         (h   (- - - - + - - x x - + - -))
##                         (kw  (- + + + - - - + x - - - -))
##                         (gw  (- + + + - - - + x + - - -))
##                         (xw  (- + + + - - - + x - + - -))


##                         ))
#
#                    symbol  feature-vector                IPA symbol......MRC.......Example
#(defvar *phonemes*    '((ee  (+ - + - - - - - + x x x x))  ;lower-case i...i.........steam
#                        (i   (+ - + - - - - - - x x x x))  ;block i........I.........pits
#                        (a_  (+ - - - - - - - + x x x x))  ;e..............eI........late
#                        (e   (+ - - - - - - - - x x x x))  ;epsilon........e.........mess
#                        (ae  (+ - - - + - - - + x x x x))  ;ae.............&.........fat
#                        (u   (+ - + + - - - - - x x x x))  ;ibar........... .........just
#                        (@   (+ - - + - - - - - x x x x))  ;inverse e...... .........burn, eatEn, Above
#                        (uu  (+ - - + + - - - - x x x x))  ;inverse v......V.........rUbber,abOve
#                        (ah  (+ - - + + - - - + x x x x))  ;a..............A.........hAlf
#                        (u_  (+ - + + - - - + + x x x x))  ;u..............u.........crUise
#                        (oo  (+ - + + - - - + - x x x x))  ;block U........U.........cOOk
#                        (o_  (+ - - + - - - + + x x x x))  ;o..............O.........cOAt
#                        (aw  (+ - - + + - - + - x x x x))  ;backwards c....0 (zero)..bAWl
#                        (y   (- - + - - - - - - x x x x))  ;j..............j.........Yell
#                        (w   (- - + + - - - + - x x x x))  ;w..............w.........Wool
#                        (r   (+ + - - - - + x x + + - -))  ;r..............r.........Rat
#                        (l   (+ + - - - + + x x + + - -))  ;l..............l.........Loin
#                        (h   (- - - - + - - x x - + - -))  ;h..............h.........Hot
#                        (p   (- + - - - + - x x - - - -))  ;p..............p.........Pot
#                        (b   (- + - - - + - x x + - - -))  ;b..............b.........But
#                        (t   (- + - - - + + x x - - - -))  ;t..............t.........Top
#                        (d   (- + - - - + + x x + - - -))  ;d..............d.........Dog
#                        (ch  (- + + - - - + x x - - - +))  ;tSigma.........tS........CHurCH
#                        (j   (- + + - - - + x x + - - +))  ;d3.............dZ........Jowls
#                        (k   (- + + + - - - - x - - - -))  ;k..............k.........Cramp
#                        (g   (- + + + - - - - x + - - -))  ;g..............g.........Goon
#                        (f   (- + - - - + - x x - + - +))  ;f..............f.........Feeler
#                        (v   (- + - - - + - x x + + - +))  ;v..............v.........Vain
#                        (th  (- + - - - + + x x - + - -))  ;theta..........T.........THought
#                        (thz (- + - - - + + x x + + - -))  ;delta..........D.........THis
#                        (s   (- + - - - + + x x - + - +))  ;s..............s.........Sign
#                        (z   (- + - - - + + x x + + - +))  ;z..............z.........Zircon
#                        (sh  (- + + - - - + x x - + - +))  ;sigma..........S.........SHeet
                                        #                        (zh  (- + + - - - + x x + + - +))  ;3..............Z.........viSIon
#                        (m   (- + - - - + - x x + - + -))  ;m..............m.........MiMe
#                        (n   (- + - - - + + x x + - + -))  ;n..............n.........Numb
#                        (ng  (- + + + - - - x x + - + -))  ;eta............9.........youNG
#                        (X   (1 1 1 1 1 1 1 1 1 1 1 1 1))  ;default placeholder
#                        (XX  (0 0 0 0 0 0 0 0 0 0 0 0 0))  ;default nil phoneme
#                        ;total:  37 phonemes +2 placeholder
#                        ))


getWord <- function(word)
    {
        words <- PHONOMAP[[ toupper(word)]]
        if(is.null(words))
            {
                warning(paste("Word [",word,"] does not exist"))
                words <- list(NA)
            }
        if(length(words)>1)
            {
                warning("more than one pronunciation found")
            }
        return(words)
    }


translateWord <- function(pword)
{
    if(!is.list(pword))
        if(is.na(pword))
            return(NA)
    
    newword <- c()
    for(i in pword)
        {
            
            newphones <- unlist(PHON2PHON[[i]])
            newword <- c(newword, newphones)

        }
 return (newword)
}


## This makes a comparison matrix, 
##weighing each feature equally.
##also, it ignores vowel/syllable stress


makeCompMatrix <- function(pword1,pword2)
{

    mat <- matrix(NA, nrow=length(pword1),ncol=length(pword2))
    for(i in 1:length(pword1))
        {

            frow <- which(features$V1==pword1[i])
            for(j in 1:length(pword2))
                {
                    fcol <- which(features$V1==pword2[j])
                    f1 <- features[frow,2:14]
                    f2 <- features[fcol,2:14]
                    ##get rid of any where both are null (x)
                    filter <- !(f1=="x" & f1==f2)
                    match <- (f1==f2)[filter]
                    mat[i,j] <- 1-mean(match)
                }
        }
    
return (mat)
}

##this will extract comp matrix from COMPCACHE; much faster.
getCompMatrix <- function( pword1, pword2)
{
    f1 <- factor(pword1,levels=features$V1)
    f2 <- factor(pword2,levels=features$V1)
    COMPCACHE[f1,f2]
}

##internal getDistance_ function; requires recoded words.
getDistance_ <- function(pword1,pword2,pmat,inscost=1,print=F)
    {
        ##pword1 is the left (rows)
        ##pword2 is the right (columns)
        planargraph <- matrix(NA, nrow=length(pword1)+1,ncol=length(pword2)+1)

        planargraph[1,] <- 0:(length(pword2))*inscost
        planargraph[,1] <- 0:(length(pword1))*inscost


        for(i in 2:(length(pword1)+1))
            {
                for(j in 2:(length(pword2)+1))
                    {

                        ##print(paste(i,j))

                        ##These are the transition costs
                        up <- planargraph[i-1,j]+inscost
                        left <- planargraph[i,j-1]+inscost
                        uleft <- planargraph[i-1,j-1]

                        ##the best cost is the cheapest + current mismatch
                        planargraph[i,j] <- min(up,left,uleft) + pmat[i-1,j-1]

                    }
            }
        if(print)
            {
                print(pmat)
                print(round(planargraph,3))
            }

        return(planargraph[length(pword1)+1,length(pword2)+1])
    }

getDistance <- function(word1,word2,inscost=1.0,pmat=NULL,print=F)
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

        getDistance_(pword1,pword2,pmat,inscost,print)        
    }


getDistance2 <- function(word1,word2,inscost=1.0,pmat=NULL,print=F)
    {
        pword1  <- translateWord(getWord(word1)[[1]])
        pword2  <- translateWord(getWord(word2)[[1]])


        print(pword1)
        print(pword2)
        
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
        
        levdist(pword1,pword2,pmat,inscost)        
    }





getDistanceDistribution <- function(word,inscost=1.0)
    {
        #stop("this takes too long.")
        pword1  <- translateWord(getWord(word)[[1]])

        out <- rep(NA,length(PHONOWORDS))
        for(i in 1:length(PHONOWORDS))
            {
                pword2  <- translateWord(getWord(PHONOWORDS[i])[[1]])
                pmat <- getCompMatrix(pword1,pword2)
                out[i] <- levdist(pword1,pword2,pmat,inscost)        
                if(i%%1000==0){cat(i,PHONOWORDS[i],"  ",pword2,"\n")}
            }

        names(out) <- PHONOWORDS
        return(out)
       
    }


sortme <- function(l,n=50){sort(l)[1:n]}



getDistance2("horse","force")
getDistance2("horse","bigot")
getDistance("horse","bigot")

getDistance2("horse","course")
getDistance("horse","course")


getDistance2("horse","choice")
getDistance("horse","choice",print=T)




##This takes way long time:

h <- getDistanceDistribution("horse")
sortme(h)    


sortme(getDistanceDistribution("bottle"))
sortme(getDistanceDistribution("insight"))
sortme(getDistanceDistribution("lemon"))



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
