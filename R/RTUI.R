#' easyRT for Shiny App
RTshiny<-function(infile=NULL, title=NULL, refGene=NULL, delim=NULL, refCond=NULL,
                  method="ddCt", std=NULL, avg=NULL, writeOut=NULL, col="Dark2",
                  showStat=T, showEB=T, bioRad=NULL){
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
    #writeOut<-as.character(readline(prompt="Would you like to write the results out to csv? [Y/N]:"))
  }
  if(writeOut=="Y" || writeOut=="y"){
    #filename<-as.character(readline(prompt="Please enter a file name for the results (ending in .csv):"))
    fe<-file.exists(filename)
    while(isTRUE(fe)){
      #message("File: ",filename," exists.")
      #ow<-as.character(readline(prompt="Overwrite? [Y/N]"))
      if(ow=="Y" || ow=="y"){
        #message(filename, " will be overwritten...")
        fe<-F
      } else {
        #filename<-as.character(readline(prompt="Please enter a different filename (.csv):"))
        fe<-file.exists(filename)
      }
    }
    #message("Writing results to ", filename,"...")
    write.csv(x, filename, row.names=F, quote=F)
  }
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
        shiny::actionButton(inputId="button", label="Generate Plot!")
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
      data()
    })
    shiny::observeEvent(input$button, {

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
                writeOut="N", bioRad=br, showStat=stat, showEB=eb)
      })
    })
  }
  shiny::shinyApp(ui = ui, server = server)
}
