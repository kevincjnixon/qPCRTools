#' easyRT for Shiny App
RTshiny<-function(infile=NULL, title=NULL, refGene=NULL, delim=NULL, refCond=NULL,
                  method="ddCt", std=NULL, avg=NULL, col="Dark2",
                  showStat=T, showEB=T, bioRad=NULL, returnDat=F){
  require(dplyr, quietly=T)
  #Check for the input file:
  if(is.null(infile)){
    #infile<-file.choose()
  }
  dat<-NULL
  form<-bioRad
  if(is.null(bioRad)){
    #form<-as.character(readline(prompt="Is this bioRad data without header? [Y/N]"))
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
    #std<-as.numeric(readline(prompt="Please enter a standard deviation threshold of Ct to filter outlying wells (recommended 0.3):"))
  }
  #message("Filtering wells based on a standard deviation of Ct of: ", std)
  y<-plyr::ddply(dat, c("Sample","Detector"), .fun=sdFilter, "Ct", std)
  #message(nrow(dat)-nrow(y)," wells removed to reduce standard deviation...")
  if(is.null(delim)){
    #print(head(levels(as.factor(y$Sample))))
    #delim<-as.character(readline(prompt="Please enter the delimiter separating Condition from Replicate (if space, leave blank):"))
    if(delim==""){
      delim<- " "
    }
  }
  #message("Using delimiter: '", delim,"' to assess conditions from 'Sample' column...")
  y$Condition<-factor(sapply(strsplit(as.character(y$Sample), split=delim, fixed=T),'[[',1))
  #Next, we want to get the reference conditions and genes
  if(is.null(refCond)){
    conds<-levels(y$Condition)
    #message("The following conditions were found: ")
    #print(paste0(1:length(conds),". ",conds))
    #refCond<-conds[as.numeric(readline(prompt="Enter the number of the reference condition:"))]
  }
  if(is.null(refGene)){
    #num<-as.numeric(readline(prompt="How many reference genes were used? (must be at least 1):"))
    conds<-levels(dat$Detector)
    #message("The following primer targets were found: ")
    #print(paste0(1:length(conds),". ",conds))
    #refGene<-conds[as.numeric(readline(prompt="Enter the number of the reference target:"))]
    if(num > 1){
      for(i in 2:num){
        #refGene<-c(refGene, conds[as.numeric(readline(prompt="Enter the number of the next reference target:"))])
      }
    }
  }
  if(is.null(avg)){
    options<-c("geoMean","mean")
    #message("The following options for averaging Ct values:")
    #print(paste0(1:length(options),". ",options))
    #avg<-options[as.numeric(readline(prompt="Enter the number indicating the averaging method to use:"))]
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
    #title<-as.character(readline(prompt="Enter a title for your plot: "))
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
  if(isTRUE(returnDat)){
    summ<-as.data.frame(x)
    summ$pvalue<-rep(NA, nrow(summ))
    summ$padj<-rep(NA, nrow(summ))
    for(i in 1:nrow(pwc)){
      summ$pvalue[i*2]<-pwc$p[i]
      summ$padj[i*2]<-pwc$p.adj[i]
    }
    return(summ)
  }
}

PIshiny<-function(infile=NULL, title=NULL, dilution=NULL, delim=NULL, input=NULL, refCond=NULL,
                 method="PI", std=NULL, avg=NULL, writeOut=NULL, col="Dark2", returnDat=NULL,
                 showStat=T, showEB=T, bioRad=NULL, condFirst=T){
  require(dplyr, quietly=T)
  #Check for the input file:
  if(is.null(infile)){
    #infile<-file.choose()
  }
  dat<-NULL
  form<-bioRad
  if(is.null(bioRad)){
  #  form<-as.character(readline(prompt="Is this bioRad data without header? [Y/N]"))
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
    #std<-as.numeric(readline(prompt="Please enter a standard deviation threshold of Ct to filter outlying wells (recommended 0.3):"))
  }
  #message("Filtering wells based on a standard deviation of Ct of: ", std)
  y<-plyr::ddply(dat, c("Sample","Detector"), .fun=sdFilter, "Ct", std)
  #message(nrow(dat)-nrow(y)," wells removed to reduce standard deviation...")

  if(is.null(dilution)){
    #dilution<-as.numeric(readline(prompt="Please enter the dilution of the input (as a percentage):"))
  }
  dilution<-log2(dilution)

  if(is.null(delim)){
    #print(head(levels(as.factor(y$Sample))))
    #delim<-as.character(readline(prompt="Please enter the delimiter separating Condition from Replicate (if space, leave blank):"))
    if(delim==""){
      delim<- " "
    }
  }
  #message("Using delimiter: '", delim,"' to assess conditions from 'Sample' column...")
  options<-unlist(strsplit(as.character(y$Sample[1]), split=delim, fixed=T))
  #print(paste(1:length(options), options, sep="   "))
  condInd<-1
  IPInd<-2
  if(isFALSE(condFirst)){
    condInd<-2
    IPInd<-1
  }
  #condInd<-as.numeric(readline(prompt="Enter the number indicating the sample condition:"))
  #IPInd<-as.numeric(readline(prompt="Enter the number indicating the precipitated/input condition:"))
  y$Condition<-factor(sapply(strsplit(as.character(y$Sample), split=delim, fixed=T),'[[',condInd))
  y$IP<-factor(sapply(strsplit(as.character(y$Sample), split=delim, fixed=T),'[[',IPInd))
  #Next, we want to get the reference conditions and genes
  if(is.null(refCond)){
    conds<-levels(y$Condition)
    #message("The following conditions were found: ")
    #print(paste0(1:length(conds),". ",conds))
    #refCond<-conds[as.numeric(readline(prompt="Enter the number of the reference condition:"))]
  }
  if(is.null(input)){
    ips<-levels(y$IP)
    #message("The following IPs were found:")
    #print(paste0(1:length(ips),". ", ips))
    #input<-ips[as.numeric(readline("Enter the number that indicates the input conditon:"))]
  }
  if(is.null(avg)){
    options<-c("geoMean","mean")
    #message("The following options for averaging Ct values:")
    #print(paste0(1:length(options),". ",options))
    #avg<-options[as.numeric(readline(prompt="Enter the number indicating the averaging method to use:"))]
  }

  #Now we have the controls, let's run the PI
  res<-PI(y, input, dilution, refCond, avg, delim, condInd, IPInd)
  x<-res$res %>% dplyr::group_by(Detector)
  x$X<-factor(paste(x$Detector,x$IP, sep=" "))
  x$Condition<-factor(x$Condition)
  x$Condition<-relevel(x$Condition, refCond)
  print(head(x))

  if(is.null(title)){
    #title<-as.character(readline(prompt="Enter a title for your plot: "))
  }
  p<-ggplot2::ggplot(x, ggplot2::aes(x=X,y=PI, fill=Condition))+
    ggplot2::geom_bar(stat="identity", color="black", position=ggplot2::position_dodge())+
    ggplot2::labs(title=title, y="Percent Input (%)", x="Gene") +
    ggplot2::theme_minimal() + ggplot2::theme(axis.text=ggplot2::element_text(angle=60, hjust=1))+
    ggplot2::scale_fill_brewer(palette=col)

  print(p)
  summ<-as.data.frame(x)
  #summ$pvalue<-rep(NA, nrow(summ))
  #summ$padj<-rep(NA, nrow(summ))
  #for(i in 1:nrow(pwc)){
  #  summ$pvalue[i*2]<-pwc$p[i]
  #  summ$padj[i*2]<-pwc$p.adj[i]
  #}
  return(summ)
}




#' User Interface Shiny App for easyRT delta-delta Ct function
#'
#' @return A Shiny app that makes ddCt analysis easier for RT-qPCR data
#' @export
RTUI<-function(){
  ui<-shiny::bootstrapPage(
    # App title ----
    shiny::titlePanel("ddCt Analysis"),

    # Sidebar layout with input and output definitions ----
    shiny::sidebarLayout(

      # Sidebar panel for inputs ----
      shiny::sidebarPanel(
        shiny::fileInput(inputId="infile",
                         label="Select File:",
                         accept=".txt"),
        shiny::radioButtons(inputId="bioRad",
                            label="Is this bioRad data without header?",
                            choices=c("Yes","No"),
                            selected="No"),
        shiny::textInput(inputId="title",
                         label="Plot Title:",
                         value="ddCt",
                         placeholder = "ddCt"),
        # Input: Numeric input for the log2FoldChange Threshold ----
        shiny::numericInput(inputId = "std",
                            label = "Standard Deviation Threshold:",
                            min = 0,
                            max = NA,
                            value = 0.3),
        # Radio buttons to determine if adjusted or raw p-value should be used
        shiny::h4("Delimiter:"),
        shiny::helpText("The character separating the condition",
                        "from the replicate in the 'Sample' column",
                        "(e.g. space, -, or _). Defaults to space."),
        shiny::textInput(inputId= "delim",
                         label=NULL,
                         value=" "),


        #Input: Numeric input for the
        shiny::uiOutput('condOpt'),
        shiny::uiOutput('geneOpt'),

        shiny::radioButtons(inputId="avg",
                            label="Averaging method?",
                            choices=c("mean","geoMean"),
                            selected="mean"),
        shiny::radioButtons(inputId="stat",
                            label="Show statistics?",
                            choices=c("Yes","No"),
                            selected="Yes"),
        shiny::radioButtons(inputId="eb",
                            label="Show Error Bars?",
                            choices=c("Yes","No"),
                            selected="Yes"),
        #shiny::actionButton(inputId="button", label="Generate Plot!"),
        shiny::radioButtons(inputId="showTab",
                            label="Show Table:",
                            choices=c("Raw data","Processed Data"),
                            selected="Raw data")
      ),

      # Main panel for displaying outputs ----
      shiny::mainPanel(
        shiny::fluidRow(
          # Output: ddCt bar graph ----
          shiny::column(12, shiny::plotOutput(outputId = "ddCt")),
          shiny::column(12, shiny::dataTableOutput(outputId = "tab"))
        )
      )
    )
  )

  server<-function(input, output){
    filename<-shiny::reactive(
      return(input$infile$datapath)
    )

    data<-reactive({
      fn<-filename()
      x<-NULL
      if(input$bioRad=="Yes"){
        x<-qPCRTools::bioRadImport(fn)
      } else {
        x<-read.delim(fn, skip=10)[,c(1:6)]
      }
      x[x==""]<-NA
      x<-x[complete.cases(x),]
      return(x)
    })
    getCond<-shiny::reactive({
      y<-data()
      #print(head(y))
      x<-as.character(levels(factor(y$Sample)))
      x<-unique(sapply(strsplit(x, split=input$delim, fixed=T),'[[',1))
      #print(x)
      return(x)
    })
    output$condOpt<-shiny::renderUI({
      opt<-getCond()
      selectInput(inputId="refCond",
                  label="Reference Condition:",
                  multiple=F,
                  choices=opt,
                  selected=opt[1])
    })
    getGenes<-shiny::reactive({
      y<-data()
      x<-as.character(levels(factor(y$Detector)))
      return(x)
    })
    output$geneOpt<-shiny::renderUI({
      opt<-getGenes()
      selectInput(inputId="refGene",
                  label="Reference Gene(s):",
                  multiple=T,
                  choices=opt,
                  selected=opt[1])
    })
    output$tab<-shiny::renderDataTable({
      if(input$showTab=="Raw data"){
        data()
      } else {
        br<-"N"
        if(input$bioRad == "Yes"){
          br<-"Y"
        }
        stat<-T
        eb<-T
        if(input$stat=="No"){
          stat<-F
        }
        if(input$eb=="No"){
          eb<-F
        }
        fn<-filename()
        RTshiny(infile=fn, title=input$title, refGene=input$refGene,
                delim=input$delim, refCond=input$refCond, std=input$std, avg=input$avg,
                returnDat=T, bioRad=br, showStat=stat, showEB=eb)
      }
    })
    output$ddCt<-shiny::renderPlot({
      br<-"N"
      if(input$bioRad == "Yes"){
        br<-"Y"
      }
      stat<-T
      eb<-T
      if(input$stat=="No"){
        stat<-F
      }
      if(input$eb=="No"){
        eb<-F
      }
      fn<-filename()
      RTshiny(infile=fn, title=input$title, refGene=input$refGene,
              delim=input$delim, refCond=input$refCond, std=input$std, avg=input$avg,
              returnDat=F, bioRad=br, showStat=stat, showEB=eb)
    })
  }
  shiny::shinyApp(ui = ui, server = server)
}

#' User Interface Shiny App for easyPI percent input function
#'
#' @return A Shiny app that makes percent input analysis easier for ChIP-qPCR data
#' @export

PIUI<-function(){
ui<-shiny::bootstrapPage(
  # App title ----
  shiny::titlePanel("Percent Input Analysis"),

  # Sidebar layout with input and output definitions ----
  shiny::sidebarLayout(

    # Sidebar panel for inputs ----
    shiny::sidebarPanel(
      shiny::fileInput(inputId="infile",
                       label="Select File:",
                       accept=".txt"),
      shiny::radioButtons(inputId="bioRad",
                          label="Is this bioRad data without header?",
                          choices=c("Yes","No"),
                          selected="No"),
      shiny::textInput(inputId="title",
                       label="Plot Title:",
                       value="ddCt",
                       placeholder = "Percent Input"),
      shiny::h4("Dilution factor:"),
      shiny::helpText("Enter dilution as whole number percentage.",
                      "(e.g. if it's 10%, enter 10)"),
      shiny::numericInput(inputId="dilution",
                          label=NULL,
                          min=0,
                          max=100,
                          value=10),
      # Input: Numeric input for the log2FoldChange Threshold ----
      shiny::numericInput(inputId = "std",
                          label = "Standard Deviation Threshold:",
                          min = 0,
                          max = NA,
                          value = 0.3),
      # Radio buttons to determine if adjusted or raw p-value should be used
      shiny::h4("Delimiter:"),
      shiny::helpText("The character separating the condition",
                      "from the replicate in the 'Sample' column",
                      "(e.g. space, -, or _). Defaults to space."),
      shiny::textInput(inputId= "delim",
                       label=NULL,
                       value=" "),
      shiny::h4("Condition First?"),
      shiny::helpText("Does the condition come before the immunoprecipitated condition",
                      "in the \"Sample\" column?"),
      shiny::radioButtons(inputId="condFirst",
                           label=NULL,
                           choices=c("Yes","No"),
                           selected="Yes"),
      #Input: Numeric input for the
      shiny::uiOutput('condOpt'),
      shiny::uiOutput('IPOpt'),

      shiny::radioButtons(inputId="avg",
                          label="Averaging method?",
                          choices=c("mean","geoMean"),
                          selected="mean"),
      shiny::radioButtons(inputId="stat",
                          label="Show statistics?",
                          choices=c("Yes","No"),
                          selected="Yes"),
      shiny::radioButtons(inputId="eb",
                          label="Show Error Bars?",
                          choices=c("Yes","No"),
                          selected="Yes"),
      #shiny::actionButton(inputId="button", label="Generate Plot!"),
      shiny::radioButtons(inputId="showTab",
                          label="Show Table:",
                          choices=c("Raw data", "Processed data"),
                          selected="Raw data")
    ),

    # Main panel for displaying outputs ----
    shiny::mainPanel(
      shiny::fluidRow(
        # Output: ddCt bar graph ----
        shiny::column(12, shiny::plotOutput(outputId = "piFig")),
        shiny::column(12, shiny::dataTableOutput(outputId = "tab"))
      )
    )
  )
)

server<-function(input, output){
  filename<-shiny::reactive(
    return(input$infile$datapath)
  )

  data<-reactive({
    fn<-filename()
    x<-NULL
    if(input$bioRad=="Yes"){
      x<-qPCRTools::bioRadImport(fn)
    } else {
      x<-read.delim(fn, skip=10)[,c(1:6)]
    }
    x[x==""]<-NA
    x<-x[complete.cases(x),]
    return(x)
  })
  getCond<-shiny::reactive({
    y<-data()
    ind<-1
    if(input$condFirst=="No"){
      ind<-2
    }
    #print(head(y))
    x<-as.character(levels(factor(y$Sample)))
    x<-unique(sapply(strsplit(x, split=input$delim, fixed=T),'[[',ind))
    #print(x)
    return(x)
  })
  output$condOpt<-shiny::renderUI({
    opt<-getCond()
    selectInput(inputId="refCond",
                label="Reference Condition:",
                multiple=F,
                choices=opt,
                selected=opt[1])
  })
  getIP<-shiny::reactive({
    y<-data()
    ind<-2
    if(input$condFirst=="No"){
      ind<-1
    }
    #print(head(y))
    x<-as.character(levels(factor(y$Sample)))
    x<-unique(sapply(strsplit(x, split=input$delim, fixed=T),'[[',ind))
    #print(x)
    return(x)
  })
  output$IPOpt<-shiny::renderUI({
    opt<-getIP()
    selectInput(inputId="IPref",
                label="Input:",
                multiple=T,
                choices=opt,
                selected=opt[1])
  })
  output$tab<-shiny::renderDataTable({
    if(input$showTab=="Raw data"){
      data()
    } else {
      br<-"N"
      if(input$bioRad == "Yes"){
        br<-"Y"
      }
      stat<-T
      eb<-T
      if(input$stat=="No"){
        stat<-F
      }
      if(input$eb=="No"){
        eb<-F
      }
      cf<-T
      if(input$condFirst=="No"){
        cf<-F
      }
      fn<-filename()
      print(c(fn, input$title, input$IPref, input$dilution, input$delim, input$refCond, input$std, input$avg,
              br, stat, eb, cf))
      PIshiny(infile=fn, title=input$title, input=input$IPref, dilution=input$dilution,
              delim=input$delim, refCond=input$refCond, std=input$std, avg=input$avg,
              writeOut="N", bioRad=br, showStat=stat, showEB=eb, condFirst=cf)
    }
  })
  output$piFig<-shiny::renderPlot({
    br<-"N"
    if(input$bioRad == "Yes"){
      br<-"Y"
    }
    stat<-T
    eb<-T
    if(input$stat=="No"){
      stat<-F
    }
    if(input$eb=="No"){
      eb<-F
    }
    cf<-T
    if(input$condFirst=="No"){
      cf<-F
    }
    fn<-filename()
    print(c(fn, input$title, input$IPref, input$dilution, input$delim, input$refCond, input$std, input$avg,
            br, stat, eb, cf))
    PIshiny(infile=fn, title=input$title, input=input$IPref, dilution=input$dilution,
            delim=input$delim, refCond=input$refCond, std=input$std, avg=input$avg,
            writeOut="N", bioRad=br, showStat=stat, showEB=eb, condFirst=cf)
  })
}
shiny::shinyApp(ui = ui, server = server)
}
