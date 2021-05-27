#' Easy RT-qPCR analysis
#'
#' This function will make RT-qPCR anlaysis easy by asking for step-by-step
#' information, including input file, reference genes, and reference conditions.
#'
#' @param infile Character indicating input file (raw from qPCR maching). Leave NULL to choose interactively.
#' @param title Character identifying the title of the plot. Leave NULL to choose interactively.
#' @param refGene Character vector indicating gene(s) to use for reference expression. Leave NULL to choose interactively.
#' @param delim Character indicating the delimiter in the 'Sample' column to separate conditions (entered first) and replicates (entered second). Leave NULL to choose interactively.
#' @param refCond Character indicating the control condition to use for normalization. Leave NULL to choose interactively.
#' @param method Character indicating analysis method. Currently only delta delta Ct method available.
#' @param std Numeric indicating standard deviation threshold before evaluating outliers. Leave NULL to choose interactively.
#' @param avg Character indicating the averaging method ('geoMean', or 'mean'). Leave NULL to choose interactively.
#' @param writeOut Character ("Y"/"N") indicating if results are to be written out to csv. Leave NULL to choose interactively.
#' @param col Character indicating the RColorBrewer palette for colour scheme. Default is "Dark2".
#' @param showStat Boolean indicating if statistics are to be shown. Default is TRUE.
#' @param showEB Boolean indicating if error bars are to be showns. Default is TRUE.
#' @return data.frame object with relative expression values and a bar plot depicting relative expression
#' @export
#'
easyRT<-function(infile=NULL, title=NULL, refGene=NULL, delim=NULL, refCond=NULL,
                 method="ddCt", std=NULL, avg=NULL, writeOut=NULL, bioRad=NULL, col="Dark2",
                 showStat=T, showEB=T){
  require(dplyr, quietly=T)
  #Check for the input file:
  if(is.null(infile)){
    infile<-file.choose()
  }
  dat<-NULL
  form<-bioRad
  if(is.null(bioRad)){
    form<-as.character(readline(prompt="Is this bioRad data without header? [Y/N]"))
  }
  if(form=="Y" || form=="y"){
    dat<-bioRadImport(infile)
  } else {
    dat<-read.delim(infile, skip=10)[,c(1:6)] #Change according to file format
  }
  #Change empt characters to NA
  dat[(dat=="")]<-NA
  dat<-dat[complete.cases(dat),]
  #Ensure everything is of the right class
  dat$Sample<-factor(dat$Sample)
  dat$Detector<-factor(dat$Detector)
  dat$Ct<-as.numeric(dat$Ct)
  #Group data into triplicates to check standard deviation
  if(is.null(std)){
    std<-as.numeric(readline(prompt="Please enter a standard deviation threshold of Ct to filter outlying wells (recommended 0.3):"))
  }
  message("Filtering wells based on a standard deviation of Ct of: ", std)
  y<-plyr::ddply(dat, c("Sample","Detector"), .fun=sdFilter, "Ct", std)
  message(nrow(dat)-nrow(y)," wells removed to reduce standard deviation...")
  if(is.null(delim)){
    print(head(levels(as.factor(y$Sample))))
    delim<-as.character(readline(prompt="Please enter the delimiter separating Condition from Replicate (if space, leave blank):"))
    if(delim==""){
      delim<- " "
    }
  }
  message("Using delimiter: '", delim,"' to assess conditions from 'Sample' column...")
  y$Condition<-factor(sapply(strsplit(as.character(y$Sample), split=delim, fixed=T),'[[',1))
  #Next, we want to get the reference conditions and genes
  if(is.null(refCond)){
    conds<-levels(y$Condition)
    message("The following conditions were found: ")
    print(paste0(1:length(conds),". ",conds))
    refCond<-conds[as.numeric(readline(prompt="Enter the number of the reference condition:"))]
  }
  if(is.null(refGene)){
    num<-as.numeric(readline(prompt="How many reference genes were used? (must be at least 1):"))
    conds<-levels(dat$Detector)
    message("The following primer targets were found: ")
    print(paste0(1:length(conds),". ",conds))
    refGene<-conds[as.numeric(readline(prompt="Enter the number of the reference target:"))]
    if(num > 1){
      for(i in 2:num){
        refGene<-c(refGene, conds[as.numeric(readline(prompt="Enter the number of the next reference target:"))])
      }
    }
  }
  if(is.null(avg)){
    options<-c("geoMean","mean")
    message("The following options for averaging Ct values:")
    print(paste0(1:length(options),". ",options))
    avg<-options[as.numeric(readline(prompt="Enter the number indicating the averaging method to use:"))]
  }

  #Now we have the controls, let's run the ddCT
  res<-ddCT(y, refGene, refCond, avg, delim)
  x<-res$res %>% dplyr::group_by(Detector)
  x$Condition<-factor(x$Condition)
  x$Condition<-relevel(x$Condition, refCond)
  #Run pairwise comparisons using forStat things (need to log2 transform)
  forStat<-res$forStat
  forStat$Condition<-relevel(forStat$Condition, refCond)
  forStat$comp<-paste(forStat$Detector, forStat$Condition, sep="_")
  forStat$comp<-factor(forStat$comp)
  if(length(levels(forStat$Condition))>2){
    forStat$comp<-relevel(forStat$comp, grep(refCond, levels(forStat$comp)))
  }
  forStat$dCt<- log2(forStat$dCt)*-1
  comps<-list()
  #If there is only two conditions (i.e. one comparison per gene):
  if(length(levels(forStat$Condition))<3){
    i<-1
    index<-1
    while(i < length(levels(forStat$comp))){
      comps[[index]]<-c(as.character(levels(forStat$comp)[i]), as.character(levels(forStat$comp)[i+1]))
      i<-i+2
      index<-index+1
    }
  } else {
    #Run pairwise comparisons for each condition to refCond
    i<-1
    index<-1
    while(i < length(levels(forStat$comp))){
      for(k in 1:length(levels(forStat$Condition))){
        if(k<length(levels(forStat$Condition))){
          comps[[index]]<-c(as.character(levels(forStat$comp)[i]), as.character(levels(forStat$comp)[i+k]))
          index<-index+1
        }
      }
      i<-i+length(levels(forStat$comp))
    }
  }
  pwc<-forStat %>% rstatix::pairwise_t_test(dCt ~ comp, comparisons=comps, p.adjust.method="BH")
  #pwc<- pwc %>% rstatix::add_xy_position(x="comp")
  pwc <- pwc %>% tibble::add_column(Detector=levels(factor(x$Detector)))
  #get y-values for pwc
  i<-1
  yvals<-c()
  while(i < length(x$ddCt)){
    yvals<-c(yvals, max(c(x$ddCt[i], x$ddCt[i+1]))*1.1)
    i<-i+2
  }
  pwc <- pwc %>% tibble::add_column(y=yvals)
  #Now we can run pairwise comparisons, normalize to control, and plot
  if(is.null(title)){
    title<-as.character(readline(prompt="Enter a title for your plot: "))
  }
  p<-ggplot2::ggplot(x, ggplot2::aes(x=Detector,y=ddCt, fill=Condition))+
    ggplot2::geom_bar(stat="identity", color="black", position=ggplot2::position_dodge())+
    ggplot2::labs(title=title, y="Relative Gene Expression", x="Gene") +
    ggplot2::theme_minimal() + ggplot2::theme(axis.text=ggplot2::element_text(angle=60, hjust=1))+
    ggplot2::scale_fill_brewer(palette=col)
  if(isTRUE(showEB)){
    p <- p + ggplot2::geom_errorbar(ggplot2::aes(ymin=ddCt-sd, ymax=ddCt+sd), width=.2, position=ggplot2::position_dodge(0.9))
  }
  if(isTRUE(showStat)){
    p<- p + ggpubr::stat_pvalue_manual(pwc, label="p.adj.signif",
                               tip.length=0, step.increase=0,
                               x="Detector", y="y", inherit.aes=F)
  }
  print(p + ggplot2::geom_hline(yintercept=1, linetype="dashed"))
  #return(res)
  if(is.null(writeOut)){
    writeOut<-as.character(readline(prompt="Would you like to write the results out to csv? [Y/N]:"))
  }
  if(writeOut=="Y" || writeOut=="y"){
    filename<-as.character(readline(prompt="Please enter a file name for the results (ending in .csv):"))
    fe<-file.exists(filename)
    while(isTRUE(fe)){
      message("File: ",filename," exists.")
      ow<-as.character(readline(prompt="Overwrite? [Y/N]"))
      if(ow=="Y" || ow=="y"){
        message(filename, " will be overwritten...")
        fe<-F
      } else {
        filename<-as.character(readline(prompt="Please enter a different filename (.csv):"))
        fe<-file.exists(filename)
      }
    }
    message("Writing results to ", filename,"...")
    write.csv(x, filename, row.names=F, quote=F)
  }
}


sdFilter<-function(dat, col,std=0.3){
  if(length(dat[[col]][complete.cases(dat[[col]])])>2){
    dev<-sd(dat[[col]], na.rm=T)
    #print(head(dat))
    #print(dev)
    if(dev>std){
      #message("Standard deviation of: ",dev,". Removing Largest outlier")
      #Calculate the differences of min and max from median
      dmin<-diff(c(median(dat[[col]]), min(dat[[col]])))
      dmax<-diff(c(max(dat[[col]]), median(dat[[col]])))
      if(dmin<dmax){
        #min is closer to median, so remove max
        dat<-dat[which(dat[[col]]!=min(dat[[col]])),]
      } else {
        #max is closer to median, so remove min
        dat<-dat[which(dat[[col]]!=max(dat[[col]])),]
      }
    }
  } else {
    dat<-dat[complete.cases(dat),]
  }
  return(dat)
}

geoMean<-function(a){
 prod(a)^(1/length(a))
}

ddCT<-function(dat, refGene, refCond, avg, delim){
  #Start by averaging Ct for each technical replicate
  meanCt<-function(x, col, avg){
    res<-NULL
    if(avg=="geoMean"){
      res<-data.frame(mean=geoMean(x[[col]]),
                      sd=sd(x[[col]]))
    }else{
      res<-data.frame(mean=mean(x[[col]]),
                      sd=sd(x[[col]]))
    }
    return(res)
  }
  dCt<-function(x, gene="Detector", mean="mean", sd="sd",refGene){
    #print(x[[mean]])
    #print(x[[mean]][which(x[[gene]] %in% refGene)])
    #dmean subtracts avg mean of test gene in rep from control gene in rep
    res<-NULL
    if(length(refGene)<2){
      message("Only 1 reference gene...")
      res<-data.frame(Detector=x[[gene]],
                      dCt=2^-(x[[mean]]-x[[mean]][which(x[[gene]] %in% refGene)]),
                      #sd=(x[[sd]]^2 + (x[[sd]][which(x[[gene]] %in% refGene)])^2)^0.5,
                      Condition=x$Condition,
                      rep=x$rep)
      #any dCt of 0, turn to NA and remove
    } else {
      message(length(refGene)," reference genes...")
      gm<-geoMean(x[[mean]][which(x[[gene]] %in% refGene)])
      res<-data.frame(Detector=x[[gene]],
                      dCt=2^-(x[[mean]]-gm),
                      Condition=x$Condition,
                      rep=x$rep)
      #Now set the reference genes to 1, so we can elminate them in the next line
      res$dCt[which(res$Detector %in% refGene)]<-1
    }
    res$dCt[res$dCt==1]<-NA
    res<-res[complete.cases(res),]
    #res$sd<-sd(res$dCt)
    return(res)
  }
  delta<-function(x, dCt="dCt", col="Condition", refCond, conditions, avg){
    #Average control dCt values and subract average from dCt of sample
    ddCt<-NULL
    sd<-NULL
    for(i in 1:length(conditions)){
      cmean<-NULL
      tmean<-NULL
      std<-NULL
      if(avg=="geoMean"){
        cmean<-geoMean(x[[dCt]][which(x[[col]] %in% refCond)])
        tmean<-geoMean(x[[dCt]][which(x[[col]] %in% conditions[i])])
        std<-sd(x[[dCt]][which(x[[col]] %in% conditions[i])])
      } else {
        cmean<-mean(x[[dCt]][which(x[[col]] %in% refCond)])
        tmean<-mean(x[[dCt]][which(x[[col]] %in% conditions[i])])
        std<-sd(x[[dCt]][which(x[[col]] %in% conditions[i])])
      }
      ddCt<-c(ddCt, tmean/cmean)
      sd<-c(sd, std/cmean/sqrt(3))
    }
    res<-data.frame(#Sample=x$Sample,
                    ddCt=ddCt,
                    sd=sd,
                    Condition=conditions)#,
                    #rep=x$rep)
    return(res)

  }
  calcExp<-function(x, col, avg){
    res<-NULL
    if(avg=="geoMean"){
      res<-data.frame(meanExp=geoMean(2^-x[[col]]),
                      sd=sd(x))
    } else {
      res<-data.frame(meanExp=mean(2^-x[[col]]),
                      sd=x$sd)
    }
    return(res)
  }
  step1<-plyr::ddply(dat, c("Sample","Detector"), .fun=meanCt, "Ct", avg)
  step1$Condition<-factor(sapply(strsplit(as.character(step1$Sample), split=delim, fixed=T),'[[',1))
  step1$rep<-factor(sapply(strsplit(as.character(step1$Sample), split=delim, fixed=T),'[[',2))
  #print(step1)
  step2<-plyr::ddply(step1, "Sample", .fun=dCt, "Detector","mean","sd",refGene)
  #print(step2)
  step3<-plyr::ddply(step2, "Detector", .fun=delta, "dCt","Condition", refCond, as.character(levels(step2$Condition)), avg)
  #res<-plyr::ddply(step3, c("Condition","Detector"), .fun=calcExp, "ddCt", avg)
  return(list(forStat=step2, res=step3, step1=step1))
}

bioRadImport<-function(filename){
  #Import bioRad formatted data and reformat it to be consistent with the format of easyRT function
  dat<-read.delim(filename)[,c(1,2,5,3,4,7)]
  colnames(dat)<-c("Position","Flag","Sample","Detector","Task","Ct")
  return(dat)
}
