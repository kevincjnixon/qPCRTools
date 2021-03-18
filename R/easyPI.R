#' Easy percent input analysis
#'
#' This function will make percent input anlaysis easy by asking for step-by-step
#' information, including input file, reference genes, and reference conditions.
#'
#' @param infile Character indicating input file (raw from qPCR maching). Leave NULL to choose interactively.
#' @param title Character identifying the title of the plot. Leave NULL to choose interactively.
#' @param dilution Numeric indicating the diltion of the input as a percentage. e.g. 1:10 = 10 (10 percent)
#' @param delim Character indicating the delimiter in the 'Sample' column to separate conditions (entered first) and replicates (entered second). Leave NULL to choose interactively.
#' @param input Character indicating the input condition
#' @param refCond Character indicating the control condition to use for normalization. Leave NULL to choose interactively.
#' @param method Character indicating analysis method. Currently only percent input method available.
#' @param std Numeric indicating standard deviation threshold before evaluating outliers. Leave NULL to choose interactively.
#' @param avg Character indicating the averaging method ('geoMean', or 'mean'). Leave NULL to choose interactively.
#' @param writeOut Character ("Y"/"N") indicating if results are to be written out to csv. Leave NULL to choose interactively.
#' @param col Character indicating the RColorBrewer palette for colour scheme. Default is "Dark2".
#' @param returnDat Character ("Y"/"N") indicating if results should be returned to console. Leave NULL to choose interactively.
#' @param showStat Boolean indicating if statistics are to be shown. Default is TRUE.
#' @param showEB Boolean indicating if error bars are to be showns. Default is TRUE.
#' @return data.frame object with relative expression values and a bar plot depicting relative expression
#' @export
#'
easyPI<-function(infile=NULL, title=NULL, dilution=NULL, delim=NULL, input=NULL, refCond=NULL,
                 method="PI", std=NULL, avg=NULL, writeOut=NULL, col="Dark2", returnDat=NULL,
                 showStat=T, showEB=T){
  require(dplyr, quietly=T)
  #Check for the input file:
  if(is.null(infile)){
    infile<-file.choose()
  }
  dat<-NULL
  form<-as.character(readline(prompt="Is this bioRad data without header? [Y/N]"))
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

  if(is.null(dilution)){
    dilution<-as.numeric(readline(prompt="Please enter the dilution of the input (as a percentage):"))
  }
  dilution<-log2(dilution)

  if(is.null(delim)){
    print(head(levels(as.factor(y$Sample))))
    delim<-as.character(readline(prompt="Please enter the delimiter separating Condition from Replicate (if space, leave blank):"))
    if(delim==""){
      delim<- " "
    }
  }
  message("Using delimiter: '", delim,"' to assess conditions from 'Sample' column...")
  options<-unlist(strsplit(as.character(y$Sample[1]), split=delim, fixed=T))
  print(paste(1:length(options), options, sep="   "))
  condInd<-as.numeric(readline(prompt="Enter the number indicating the sample condition:"))
  IPInd<-as.numeric(readline(prompt="Enter the number indicating the precipitated/input condition:"))
  y$Condition<-factor(sapply(strsplit(as.character(y$Sample), split=delim, fixed=T),'[[',condInd))
  y$IP<-factor(sapply(strsplit(as.character(y$Sample), split=delim, fixed=T),'[[',IPInd))
  #Next, we want to get the reference conditions and genes
  if(is.null(refCond)){
    conds<-levels(y$Condition)
    message("The following conditions were found: ")
    print(paste0(1:length(conds),". ",conds))
    refCond<-conds[as.numeric(readline(prompt="Enter the number of the reference condition:"))]
  }
  if(is.null(input)){
    ips<-levels(y$IP)
    message("The following IPs were found:")
    print(paste0(1:length(ips),". ", ips))
    input<-ips[as.numeric(readline("Enter the number that indicates the input conditon:"))]
  }
  if(is.null(avg)){
    options<-c("geoMean","mean")
    message("The following options for averaging Ct values:")
    print(paste0(1:length(options),". ",options))
    avg<-options[as.numeric(readline(prompt="Enter the number indicating the averaging method to use:"))]
  }

  #Now we have the controls, let's run the PI
  res<-PI(y, input, dilution, refCond, avg, delim, condInd, IPInd)
  x<-res$res %>% dplyr::group_by(Detector)
  x$X<-factor(paste(x$Detector,x$IP, sep=" "))
  x$Condition<-factor(x$Condition)
  x$Condition<-relevel(x$Condition, refCond)
  #Run pairwise comparisons using forStat things (need to log2 transform)
  #forStat<-res$forStat
  #forStat$Condition<-relevel(forStat$Condition, refCond)
  #forStat$comp<-paste(forStat$Detector, forStat$Condition, sep="_")
  #forStat$comp<-factor(forStat$comp)
  #if(length(levels(forStat$Condition))>2){
  #  forStat$comp<-relevel(forStat$comp, grep(refCond, levels(forStat$comp)))
  #}
  #forStat$dCt<- log2(forStat$dCt)*-1
  #comps<-list()
  #If there is only two conditions (i.e. one comparison per gene):
  #if(length(levels(forStat$Condition))<3){
  #  i<-1
  #  index<-1
  #  while(i < length(levels(forStat$comp))){
  #    comps[[index]]<-c(as.character(levels(forStat$comp)[i]), as.character(levels(forStat$comp)[i+1]))
  #    i<-i+2
  #    index<-index+1
  #  }
  #} else {
    #Run pairwise comparisons for each condition to refCond
  #  i<-1
  #  index<-1
  #  while(i < length(levels(forStat$comp))){
  #    for(k in 1:length(levels(forStat$Condition))){
  #      if(k<length(levels(forStat$Condition))){
  #        comps[[index]]<-c(as.character(levels(forStat$comp)[i]), as.character(levels(forStat$comp)[i+k]))
  #        index<-index+1
  #      }
  #    }
  #   i<-i+length(levels(forStat$comp))
  #  }
  #}
  #pwc<-forStat %>% rstatix::pairwise_t_test(dCt ~ comp, comparisons=comps, p.adjust.method="BH")
  #pwc<- pwc %>% rstatix::add_xy_position(x="comp")
  #pwc <- pwc %>% tibble::add_column(Detector=levels(factor(x$Detector)))
  #get y-values for pwc
  # i<-1
  # yvals<-c()
  # while(i < length(x$ddCt)){
  #   yvals<-c(yvals, max(c(x$ddCt[i], x$ddCt[i+1]))*1.1)
  #   i<-i+2
  # }
  # pwc <- pwc %>% tibble::add_column(y=yvals)
  #Now we can run pairwise comparisons, normalize to control, and plot
  if(is.null(title)){
    title<-as.character(readline(prompt="Enter a title for your plot: "))
  }
  p<-ggplot2::ggplot(x, ggplot2::aes(x=X,y=PI, fill=Condition))+
    ggplot2::geom_bar(stat="identity", color="black", position=ggplot2::position_dodge())+
    ggplot2::labs(title=title, y="Percent Input (%)", x="Gene") +
    ggplot2::theme_minimal() + ggplot2::theme(axis.text=ggplot2::element_text(angle=60, hjust=1))+
    ggplot2::scale_fill_brewer(palette=col)
  #if(isTRUE(showEB)){
  #  p <- p + ggplot2::geom_errorbar(ggplot2::aes(ymin=ddCt-sd, ymax=ddCt+sd), width=.2, position=ggplot2::position_dodge(0.9))
  #}
  #if(isTRUE(showStat)){
  #  p<- p + ggpubr::stat_pvalue_manual(pwc, label="p.adj.signif",
  #                                     tip.length=0, step.increase=0,
  #                                     x="Detector", y="y", inherit.aes=F)
  #}
  print(p)
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
  if(is.null(returnDat)){
    returnDat<-as.character(readline(prompt="Would you like to print results to the console? [Y/N]:"))
  }
  if(returnDat=="Y" || returnDat=="y"){
    return(x)
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

PI<-function(dat, input, dilution, refCond, avg, delim, condInd, IPInd){
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
  getPI<-function(x, gene="Detector",IP="IP", mean="mean", sd="sd",input, dilution){
    #print(x[[mean]])
    #print(x[[mean]][which(x[[gene]] %in% refGene)])
    #dmean subtracts avg mean of test gene in rep from control gene in rep
    res<-NULL
    dip<-x[[mean]][which(x[[IP]] %in% input)]-dilution
    PI<-c()
    eb<-c()
    for(i in 1:length(x[[mean]])){
      PI<-c(PI, 100*(2^((dip-x[[mean]][i]))))
    }
    res<-data.frame(Detector=x[[gene]],
                    PI=PI,
                    #sd=(x[[sd]]^2 + (x[[sd]][which(x[[gene]] %in% refGene)])^2)^0.5,
                    Condition=x$Condition,
                    IP=x$IP)
    res$PI[res$IP==input]<-NA
    res<-res[complete.cases(res),]
    #res$sd<-sd(res$dCt)
    return(res)
  }

  step1<-plyr::ddply(dat, c("Sample","Detector"), .fun=meanCt, "Ct", avg)
  step1$Condition<-factor(sapply(strsplit(as.character(step1$Sample), split=delim, fixed=T),'[[',condInd))
  step1$IP<-factor(sapply(strsplit(as.character(step1$Sample), split=delim, fixed=T),'[[',IPInd))
  #print(step1)
  step2<-plyr::ddply(step1, c("Condition","Detector"), .fun=getPI,"Detector", "IP","mean","sd",input, dilution)
  return(list(res=step2,step1=step1))
}

bioRadImport<-function(filename){
  #Import bioRad formatted data and reformat it to be consistent with the format of easyRT function
  dat<-read.delim(filename)[,c(1,2,5,3,4,7)]
  colnames(dat)<-c("Position","Flag","Sample","Detector","Task","Ct")
  return(dat)
}
