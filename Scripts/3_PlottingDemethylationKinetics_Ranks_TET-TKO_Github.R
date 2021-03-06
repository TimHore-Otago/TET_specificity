#######################################################

#######Plot Demethy kinetics for ranked motifs#########
#######         a.k.a. the chicken plot       #########

#######################################################
#load dplyr package to allow select function
library(dplyr)


setwd("~/Dropbox/DemethylationKinetics_Scripts/Oscar/Datasets")

data <- read.delim("ESC_TET_TKO_TET3_meth.csv_DemethVeloc.txt")
data <- data[order(data$slope),]

#load in the methylation data of interest
#the . is a wildcard representing A, B or C replicate
zero.h <- t(select(data, matches("ESC_TET_TKO_TET3_Dox_0h_._Perc_met")))
six.h <- t(select(data, matches("ESC_TET_TKO_TET3_Dox_6h_._Perc_met")))
twelve.h <- t(select(data, matches("ESC_TET_TKO_TET3_Dox_12h_._Perc_met")))
eighteen.h <- t(select(data, matches("ESC_TET_TKO_TET3_Dox_18h_._Perc_met")))
twenty4.h <- t(select(data, matches("ESC_TET_TKO_TET3_Dox_24h_._Perc_met")))
thirty.h <- t(select(data, matches("ESC_TET_TKO_TET3_Dox_30h_._Perc_met")))
thirty6.h <- t(select(data, matches("ESC_TET_TKO_TET3_Dox_36h_._Perc_met")))
forty2.h <- t(select(data, matches("ESC_TET_TKO_TET3_Dox_42h_._Perc_met")))
forty8.h <- t(select(data, matches("ESC_TET_TKO_TET3_Dox_48h_._Perc_met")))
fifty4.h <- t(select(data, matches("ESC_TET_TKO_TET3_Dox_54h_._Perc_met")))
sixty.h <- t(select(data, matches("ESC_TET_TKO_TET3_Dox_60h_._Perc_met")))
sixty6.h <- t(select(data, matches("ESC_TET_TKO_TET3_Dox_66h_._Perc_met")))
seventy2.h <- t(select(data, matches("ESC_TET_TKO_TET3_Dox_72h_._Perc_met")))

####Either plot all motifs, or uncomment particular motifs of interest. Also be careful to choose error bars or not (see "arrows" below)
for (Motif.col in c(1:256)){  #THis gives the full chicken plot.....(i.e. all 256 motifs)
#for (Motif.col in c(1,50,100,150,200,250)){  #THis is for the 1,50,100,150,200,250 rankedset 
#for (Motif.col in c(2,13)){  #These are MCGW pairs in the top 20
#for (Motif.col in c(3,6)){  #These are MCGW pairs in the top 20    
#for (Motif.col in c(7,105)){  #THis pair of complementary motifs have the largest difference in rate
#for (Motif.col in c(15,118)){  #THis pair of complementary motifs have the 2nd largest difference in rate
#for (Motif.col in c(1:10)){  #THis is for the 1,50,100,150,200,250 rankedset  
    
    command.a <- paste("Motif1 <- zero.h[,",Motif.col,"]")
    eval(parse(text=command.a))
    command.a <- paste("Motif1 <- cbind(Motif1, six.h[,",Motif.col,"])")
    eval(parse(text=command.a))
    command.a <- paste("Motif1 <- cbind(Motif1, twelve.h[,",Motif.col,"])")
    eval(parse(text=command.a))
    command.a <- paste("Motif1 <- cbind(Motif1, eighteen.h[,",Motif.col,"])")
    eval(parse(text=command.a))
    command.a <- paste("Motif1 <- cbind(Motif1, twenty4.h[,",Motif.col,"])")
    eval(parse(text=command.a))
    command.a <- paste("Motif1 <- cbind(Motif1, thirty.h[,",Motif.col,"])")
    eval(parse(text=command.a))
    command.a <- paste("Motif1 <- cbind(Motif1, thirty6.h[,",Motif.col,"])")
    eval(parse(text=command.a))
    command.a <- paste("Motif1 <- cbind(Motif1, forty2.h[,",Motif.col,"])")
    eval(parse(text=command.a))
    command.a <- paste("Motif1 <- cbind(Motif1, forty8.h[,",Motif.col,"])")
    eval(parse(text=command.a))
    command.a <- paste("Motif1 <- cbind(Motif1, fifty4.h[,",Motif.col,"])")
    eval(parse(text=command.a))
    command.a <- paste("Motif1 <- cbind(Motif1, sixty.h[,",Motif.col,"])")
    eval(parse(text=command.a))
    command.a <- paste("Motif1 <- cbind(Motif1, sixty6.h[,",Motif.col,"])")
    eval(parse(text=command.a))
    command.a <- paste("Motif1 <- cbind(Motif1, seventy2.h[,",Motif.col,"])")
    eval(parse(text=command.a))

    #add column names to Motif1 table
    colnames(Motif1)<- c("0 h", "6 h", "12 h", "18 h", "24 h" , "30 h", "36 h", "42 h", "48 h", "54 h", "60 h", "66 h", "72 h")
    h <- c(0,6,12,18,24,30,36,42,48,54,60,66,72)
    
    
    #Convert the table from factors to characters to numbers
    Motif1 <- as.data.frame(Motif1, stringsAsFactors = FALSE)
    Motif1 <- as.data.frame(sapply(Motif1, as.numeric))    

    #Create difference data
    Loss.pp <- Motif1 - mean(Motif1$`0 h`) 
    Loss.pp
    
    #Plot line graph with error bars
    mean.data = colMeans(Loss.pp, na.rm=T)
    no.of.samples = apply(!is.na(Loss.pp), 2, sum)
    st.dev = sapply(Loss.pp, sd, na.rm=T)
    par(new=TRUE)
    midpoints = plot(x=h, y=mean.data, type="l", ylim = c(-30,5), las=1, frame.plot=F, ylab = "CG methylation (%)")
    ##either comment in or out if you want to plot error bars
    #arrows (x0 = h, y0 = mean.data - st.dev, x1 = h, y1 = mean.data + st.dev, angle = 90, length = 0.02, code = 3)
    
    # Add motif labels if you want them
    #command.b <- paste("text(x=70, y=mean.data[13], labels = data[",Motif.col,",1], cex=1)")
    #eval(parse(text=command.b))

  }

