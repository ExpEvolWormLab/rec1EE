# Change chromosomes id Roman <-> num


#setwd("/Users/tomparee/Desktop/rec1EE_github/")

#path = "./Data/rec1EE/" #path on computer perso

recpos = 719556 # position of the deletion in rec-1 mutant on WS245 genome

numtorom = function(vector_chromid){
  vector_chromid = as.character(as.roman(vector_chromid))
  vector_chromid[vector_chromid=="VI"]="X"
  vector_chromid = factor(vector_chromid, levels = c("I","II","III","IV","V","X"))
  return(vector_chromid)
}

romtonum = function(vector_chromid){
  vector_chromid[vector_chromid=="I"]=1
  vector_chromid[vector_chromid=="II"]=2
  vector_chromid[vector_chromid=="III"]=3
  vector_chromid[vector_chromid=="IV"]=4
  vector_chromid[vector_chromid=="V"]=5
  vector_chromid[vector_chromid=="X"]=6
  return(vector_chromid)
}


get.win = function(size, winsize, minsize){
  
  if(winsize < size){
    
    nblocks = size/winsize
    fullblocks = floor(nblocks)
    
    ## split column numbers into 'nblocks' groups
    SPLIT <- split(1:(fullblocks*winsize), rep(1:fullblocks, each = winsize))
    if(nblocks>fullblocks){ 
      dblock = (nblocks-fullblocks)*winsize
      SPLIT[[length(SPLIT)+1]] = 1:dblock + max(fullblocks*winsize)
    }
    
    if(length(SPLIT[[length(SPLIT)]]) < minsize & length(SPLIT) > 1){
      SPLIT[[length(SPLIT)-1]] = c(SPLIT[[length(SPLIT)]],SPLIT[[length(SPLIT)-1]])
      SPLIT = SPLIT[1:(length(SPLIT)-1)]
    }
    
    
  }else{
    
    SPLIT = list(1:size)
  }
  
  return(SPLIT)
}


se = function(x){sd(x)/sqrt(sum(!is.na(x)))}



assign.info.melt = function(melted.freq.matrix, haplo = F){
  require(data.table)
  names(melted.freq.matrix) = c('marker', 'sample', 'freq')
  melted.freq.matrix$chrom = tstrsplit(melted.freq.matrix$marker, '_')[[1]]
  melted.freq.matrix$pos = as.numeric(tstrsplit(melted.freq.matrix$marker, '_')[[2]])
  
  if(grepl("_G", melted.freq.matrix$sample[1])){
    melted.freq.matrix$generation = as.numeric(tstrsplit(melted.freq.matrix$sample, '_G')[[2]])
    melted.freq.matrix$pop = tstrsplit(melted.freq.matrix$sample, '_G')[[1]]
  }else{names(melted.freq.matrix)[2]='pop'}
  
  melted.freq.matrix = assign.infopop(melted.freq.matrix)
  
  if(haplo == T){names(melted.freq.matrix)[1] = 'haplo'
  melted.freq.matrix = melted.freq.matrix[,-which(names(melted.freq.matrix)=="pos")] 
  }
  
  return(melted.freq.matrix)
}

#
arraytomatrix = function(arr,dimtomerge=NULL){
  require(abind)
  if(is.null(dimtomerge)) dimtomerge= 2:length(dim(arr))
  dims=dim(arr)[dimtomerge]
  dims=lapply(dims,function(x){1:x})
  
  dcomb = expand.grid(dims[1:(length(dims)-1)])
  tardim = dimtomerge[1:(length(dimtomerge)-1)]
  dn = dimnames(arr)
  
  mx = NULL
  for(i in 1:nrow(dcomb)){
    tar = dcomb[i,]
    tarl=list()
    for(n in 1:length(tar)){tarl[n]=dn[[tardim[n]]][unlist(tar[n])]}
    x=asub(arr,tarl,tardim)
    colnames(x) = paste0(colnames(x),"_",paste(unlist(tarl), collapse='_'))
    mx=cbind(mx,as.matrix(x))
  }
  return(mx)
}

# Asign populations informations
assign.infopop.Rhom = function(data, col.pop.name='pop'){
  x = which(colnames(data)==col.pop.name)
  data$rec[grepl('LR',data[,x])]='mut'
  data$rec[grepl('HR',data[,x])]='wt'
  data$env[grepl('S',data[,x])]='salt'
  data$env[!grepl('S',data[,x])]='domestication'
  data$rec = factor(data$rec, levels = c("wt", "mut"))
  data$env = factor(data$env, levels = c("domestication", "salt"))
  return(data)
}




#Just a function which do all combination of a vector c(1,2,3) => 1-1,1-2,1-3, 2-2,2-3, 3-3

vector.combination <- function(t, with.itself=T, param.name=c('r1','r2')){
  if(with.itself){
    r1=NULL
    r2 = NULL
    for(n in 1:length(t)){
      r1 = c(r1, rep(t[n],((length(t)+1)-(1:length(t)))[n]))
      r2=c(r2, t[n:length(t)])
    }
  }else{
    r1= (length(t))-(1:length(t))
    r1=NULL
    r2 = NULL
    for(n in 1:length(t)){
      r1 = c(r1, rep(t[n],(length(t)-(1:length(t)))[n]))
      r2=c(r2, t[(n+1):length(t)])
    }
    r2=r2[1:length(r1)]
  }
  
  r = data.frame(r1,r2)
  colnames(r)=param.name
  
  return(r)
}




do.bins = function(snps, binlength=100000, bp=T){
  binsnp=function(x,nsnp){
    #nsnp =300
    x$bin=NA
    for(n in 1:(floor(nrow(x)/nsnp))){
      x$bin[(((n-1)*nsnp)+1):(n*nsnp)]=n
    }
    x$bin[which(is.na(x$bin))]=max(x$bin, na.rm=T)
    return(x)
  }
  
  binbp=function(x,nbp){
    require(birk)
    size=max(x$POS, na.rm=T)
    x$bin=NA
    for(n in 1:(floor(size/nbp))){
      start=(n-1)*nbp
      end=n*nbp
      x$bin[x$POS>=start & x$POS<end]=n
    }
    x$bin[which(is.na(x$bin))]=max(x$bin, na.rm=T)
    return(x)
  }
  
  
  if(bp){x <- do.call(rbind, lapply(split(snps, snps$chrom), function(i) {i=binbp(i,binlength); i}))
  }else{
    x <- do.call(rbind, lapply(split(snps, snps$chrom), function(i) {i=binsnp(i,binlength); i}))
  }
  
  x
  
}







##########################################################
###### GGPLOT

# Some color for colorblind
cbs = c(rgb(0,0.45,0.7),
        rgb(0.9, 0.6, 0),
        rgb(0,0.6,0.5),
        rgb(0.8, 0.4, 0),
        rgb(0.35, 0.7, 0.9),
        rgb(0.8, 0.6, 0.7),
        rgb(0.95, 0.9, 0.25),
        rgb(0,0,0,),
        "gray45")

blueWT = cbs[1]
orangeMut = cbs[2]


# function to have axis is base 10 exponant 
# from: https://stackoverflow.com/questions/10762287/how-can-i-format-axis-labels-with-exponents-with-ggplot2-and-scales
scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}

# ex: +scale_y_continuous(label=scientific_10)



# modify from: https://rdrr.io/github/HanjoStudy/quotidieR/man/theme_Publication.html
theme_Publication3 <- function(base_size=7) {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size)
    + theme(plot.title = element_text(face = "bold",
                                      hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "plain"),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.title = element_text(face="plain"),
            plot.margin=unit(c(t=5,r=5,b=5,l=5),"pt"),
            strip.background=element_rect(colour=NA,fill=NA),
            strip.text = element_text(face="plain")
    ))
  
}



# from: https://stackoverflow.com/questions/44688623/adding-custom-images-to-ggplot-facets
annotation_custom2 <- function (grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, data) 
{
  layer(data = data, stat = StatIdentity, position = PositionIdentity, 
        geom = ggplot2:::GeomCustomAnn,
        inherit.aes = TRUE, params = list(grob = grob, 
                                          xmin = xmin, xmax = xmax, 
                                          ymin = ymin, ymax = ymax))
}


