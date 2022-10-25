library(shiny)
library(magrittr)
library(dplyr)
library(BiocManager)

shinyServer(function(input, output, session) {

  options(shiny.maxRequestSize=30*1024^2)
  ranges <- reactiveValues(x = NULL, y = NULL)

  observeEvent(input$double_click, {
    brush <- input$brush_plot
    if (!is.null(brush)) {
      ranges$x <- c(brush$xmin, brush$xmax)
      ranges$y <- c(brush$ymin, brush$ymax)

    } else {
      ranges$x <- NULL
      ranges$y <- NULL
    }
  })

  output$plot_data = renderPlot({
    tryCatch(
      {
        LimaBeanGBS = vcfR::read.vcfR(input$filevcf$datapath)
      },
      error = function(e){
        stop('Upload VCF file')
      }
    )

    tryCatch(
      {
        data = read.table(input$filetxt$datapath, header = TRUE)
      },
      error = function(e){
        stop('Upload TXT file')
      }
    )

    Genepool <- as.character(data$Genepool)
    LimaBeanData2 <- vcfR::vcfR2genind(LimaBeanGBS,
                                       pop= Genepool,
                                       NA.char= "NA")

    EuclideanDistance <- dist(LimaBeanData2,
                              method = "euclidean",
                              diag = FALSE,
                              upper = FALSE,
                              p=2)
    NJtree <- ape::nj(EuclideanDistance)

    mycol = c("light blue",
              "gray",
              "green",
              "red",
              "blue",
              "light green",
              "pink")[LimaBeanData2$pop]

    plot(NJtree, tip.color=mycol, type="fan",
         x.lim = ranges$x,
         y.lim = ranges$y)
  })

  ranges1 <- reactiveValues(x = NULL, y = NULL)

  output$plot_data2 = renderPlot({
    tryCatch(
      {
        LimaBeanGBS = vcfR::read.vcfR(input$filevcf$datapath)
      },
      error = function(e){
        stop('Upload VCF file')
      }
    )

    tryCatch(
      {
        data = read.table(input$filetxt$datapath, header = TRUE)
      },
      error = function(e){
        stop('Upload TXT file')
      }
    )

    Genepool <- as.character(data$Genepool)
    LimaBeanData2 <- vcfR::vcfR2genind(LimaBeanGBS,
                                       pop= Genepool,
                                       NA.char= "NA")

    LimaBeanData3 <- vcfR::vcfR2genlight(LimaBeanGBS) # convert the vcf file into a genlight object
    adegenet::ploidy(LimaBeanData3) <- 2
    adegenet::pop(LimaBeanData3) <- Genepool

    UPGMAtree <- poppr::aboot(LimaBeanData3, tree = "upgma", distance = poppr::bitwise.dist, sample = 100, showtree = F, cutoff = 50, quiet = T)

    mycol = c("light blue",
              "gray",
              "green",
              "red",
              "blue",
              "light green",
              "pink")[LimaBeanData2$pop]

    plot(UPGMAtree, tip.color=mycol, cex = 1.2,
         x.lim = ranges1$x,
         y.lim = ranges1$y)
  })

  observeEvent(input$double_click2, {
    brush <- input$brush_plot
    if (!is.null(brush)) {
      ranges1$x <- c(brush$xmin, brush$xmax)
      ranges1$y <- c(brush$ymin, brush$ymax)

    } else {
      ranges1$x <- NULL
      ranges1$y <- NULL
    }
  })

  output$distance_matrix1 <- renderPlot({
    tryCatch(
      {
        LimaBeanGBS = vcfR::read.vcfR(input$filevcf$datapath)
      },
      error = function(e){
        stop(safeError(e))
      }
    )

    tryCatch(
      {
        data = read.table(input$filetxt$datapath, header = TRUE)
      },
      error = function(e){
        stop(safeError(e))
      }
    )
    Genepool <- as.character(data$Genepool)
    LimaBeanData2 <- vcfR::vcfR2genind(LimaBeanGBS,
                                       pop = Genepool,
                                       NA.char= "NA")

    adegenet::strata(LimaBeanData2) <- data.frame(adegenet::pop(LimaBeanData2)) #to add the strata of a genind object
    adegenet::nameStrata(LimaBeanData2) <- ~Genepool # to assign a name to the strata of a genind object
    #Building of a UPGMA topology with bootstrap support using the function aboot of the poppr package
    set.seed(999)
    poppr::aboot(LimaBeanData2, strata = Genepool, tree= "upgma", distance="nei.dist", sample=1000, cutoff=0, showtree=TRUE,quiet=FALSE)

  })

  output$distance_matrix2 <- renderPlot({
    tryCatch(
      {
        LimaBeanGBS = vcfR::read.vcfR(input$filevcf$datapath)
      },
      error = function(e){
        stop(safeError(e))
      }
    )

    tryCatch(
      {
        data = read.table(input$filetxt$datapath, header = TRUE)
      },
      error = function(e){
        stop(safeError(e))
      }
    )
    Genepool <- as.character(data$Genepool)
    LimaBeanData2 <- vcfR::vcfR2genind(LimaBeanGBS, pop= Genepool, NA.char= "NA")

    adegenet::strata(LimaBeanData2) <- data.frame(adegenet::pop(LimaBeanData2)) #to add the strata of a genind object
    adegenet::nameStrata(LimaBeanData2) <- ~Genepool # to assign a name to the strata of a genind object
    #Building of a UPGMA topology with bootstrap support using the function aboot of the poppr package
    set.seed(999)
    #Building of a NJ topology with bootstrap support using the function aboot of the poppr package
    poppr::aboot(LimaBeanData2, strata = Genepool, tree= "nj", distance="nei.dist", sample=1000, cutoff=0, showtree=TRUE,quiet=FALSE)
  })

  ######### PCA

  output$plotPCA <- renderPlot({
    tryCatch(
      {
        LimaBeanGBS = vcfR::read.vcfR(input$filevcf$datapath)
      },
      error = function(e){
        stop('Upload VCF file')
      }
    )

    tryCatch(
      {
        data = read.table(input$filetxt$datapath, header = TRUE)
      },
      error = function(e){
        stop('Upload TXT file')
      }
    )

    Genepool = as.character(data$Genepool)
    LimaBeanData3 <- vcfR::vcfR2genlight(LimaBeanGBS)
    adegenet::ploidy(LimaBeanData3) <- 2
    adegenet::pop(LimaBeanData3) <- Genepool

    LimaBeanPCA <- adegenet::glPca(LimaBeanData3, nf = 3)
    # to carry out a PCA on a genlight object. With the argument nf as NULL,
    #you are asked interactively for the number of principal components to be
    #retained. For this data, three axes were retained.

    # to create a customized PCA plot with the package ggplot2

    Limapcascores <- as.data.frame(LimaBeanPCA$scores)
    Limapcascores$pop <- adegenet::pop(LimaBeanData3)


    set.seed(5)
    colors <- RColorBrewer::brewer.pal(n = adegenet::nPop(LimaBeanData3), name = "Set1")
    p <- ggplot2::ggplot(Limapcascores, ggplot2::aes(x=PC1, y=PC2, colour=pop))
    p <- p + ggplot2::geom_hline(yintercept = 0)
    p <- p + ggplot2::geom_vline(xintercept = 0)
    p <- p + ggplot2::geom_point(size=3)
    p <- p + ggplot2::theme_bw()
    p <- p + ggplot2::scale_color_manual(values=colors ) +
      ggplot2::xlab(sprintf("PC1 %f percent", 100*LimaBeanPCA$eig[1]/sum(LimaBeanPCA$eig))) +
      ggplot2::ylab(sprintf("PC2 %f percent", 100*LimaBeanPCA$eig[2]/sum(LimaBeanPCA$eig)))
    p

  })

  ############### DPCA

  output$plotDPCA <- renderPlot({
    tryCatch(
      {
        LimaBeanGBS <- vcfR::read.vcfR(input$filevcf$datapath)
      },
      error = function(e){
        stop('Upload VCF file') #safeError(e)
      }
    )

    tryCatch(
      {
        PopData = read.table(input$filetxt$datapath, header = T)
      },
      error = function(e){
        stop('Upload VCF file')
      }
    )
    Genepool <- as.character(PopData$Genepool)

    require("adegenet")

    LimaBeanData2 <- vcfR::vcfR2genind(LimaBeanGBS, pop= Genepool, NA.char= "NA")
    grp <- find.clusters(LimaBeanData2, max.n.clust = 30, n.pca = input$npca, n.clust = input$nda)
    dapc <- dapc(LimaBeanData2, grp$grp, n.pca = input$npca, n.da = input$nda)
    scatter(dapc)
  })

  output$table1D <- renderDataTable({
    tryCatch(
      {
        LimaBeanGBS <- vcfR::read.vcfR(input$filevcf$datapath)
      },
      error = function(e){
        stop('Upload VCF file') #safeError(e)
      }
    )

    tryCatch(
      {
        PopData = read.table(input$filetxt$datapath, header = T)
      },
      error = function(e){
        stop('Upload VCF file')
      }
    )
    Genepool <- as.character(PopData$Genepool)

    require("adegenet")

    LimaBeanData2 <- vcfR::vcfR2genind(LimaBeanGBS, pop= Genepool, NA.char= "NA")
    grp <- find.clusters(LimaBeanData2, max.n.clust = 30, n.pca = input$npca1, n.clust = input$nclust)
    table(pop(LimaBeanData2), grp$grp)
  })

  output$plotDPCA2 <- renderPlot({
    tryCatch(
      {
        LimaBeanGBS <- vcfR::read.vcfR(input$filevcf$datapath)
      },
      error = function(e){
        stop('Upload VCF file') #safeError(e)
      }
    )

    tryCatch(
      {
        PopData = read.table(input$filetxt$datapath, header = T)
      },
      error = function(e){
        stop('Upload TXT file')
      }
    )
    Genepool <- as.character(PopData$Genepool)

    require("adegenet")

    LimaBeanData2 <- vcfR::vcfR2genind(LimaBeanGBS, pop= Genepool, NA.char= "NA")
    grp <- find.clusters(LimaBeanData2, max.n.clust = 30, n.pca = input$npca1, n.clust = input$nclust)
    dapc <- dapc(LimaBeanData2, grp$grp, n.pca = input$npca, n.da = input$nda)
    myCol2 <- c("pink","red","blue","light blue", "green") # to assign colors to each of the five clusters

    scatter(dapc, scree.da=FALSE, bg="white", pch=20,  cell=0, cstar=0, col=myCol2, solid=1.0,
            cex=3,clab=0, leg=TRUE, posi.leg= "bottomleft", scree.pca=TRUE, posi.pca = "topright", ratio.pca=0.3)
  })

  output$plotDPCA3 <- renderPlot({
    tryCatch(
      {
        LimaBeanGBS <- vcfR::read.vcfR(input$filevcf$datapath)
      },
      error = function(e){
        stop('Upload VCF file') #safeError(e)
      }
    )

    tryCatch(
      {
        PopData = read.table(input$filetxt$datapath, header = T)
      },
      error = function(e){
        stop('Upload TXT file')
      }
    )
    Genepool <- as.character(PopData$Genepool)

    require("adegenet")

    LimaBeanData2 <- vcfR::vcfR2genind(LimaBeanGBS, pop= Genepool, NA.char= "NA")
    grp <- find.clusters(LimaBeanData2, max.n.clust = 30, n.pca = input$npca1, n.clust = input$nclust)
    dapc <- dapc(LimaBeanData2, grp$grp, n.pca = input$npca, n.da = input$nda)
    compoplot.dapc(dapc)
  })


  ############### BDSpl

  output$plot_d1 <- renderPlot({

    tryCatch(
      {
        LimaBeanGBS <- vcfR::read.vcfR(input$filevcf$datapath)
      },
      error = function(e){
        stop('Upload VCF file') #safeError(e)
      }
    )

    tryCatch(
      {
        PopData = read.table(input$filetxt$datapath, header = T)
      },
      error = function(e){
        stop('Upload TXT file')
      }
    )
    require("adegenet")
    Genepool <- as.character(PopData$Genepool)
    LimaBeanData2 <- vcfR::vcfR2genind(LimaBeanGBS, pop= Genepool, NA.char= "NA")
    diversity <- summary(LimaBeanData2)
    plot(diversity$Hobs, xlab="loci number", ylab="Hobs", main="Hobs per locus")
  })

  output$plot_d2 <- renderPlot({

    tryCatch(
      {
        LimaBeanGBS <- vcfR::read.vcfR(input$filevcf$datapath)
      },
      error = function(e){
        stop('Upload VCF file') #safeError(e)
      }
    )

    tryCatch(
      {
        PopData = read.table(input$filetxt$datapath, header = T)
      },
      error = function(e){
        stop('Upload TXT file')
      }
    )
    require("adegenet")
    Genepool <- as.character(PopData$Genepool)
    LimaBeanData2 <- vcfR::vcfR2genind(LimaBeanGBS, pop= Genepool, NA.char= "NA")
    diversity <- summary(LimaBeanData2)
    plot(diversity$Hexp, xlab="loci number", ylab="Hexp", main="Hexp per locus")
  })

  output$plot_d3 <- renderPlot({
    tryCatch(
      {
        LimaBeanGBS <- vcfR::read.vcfR(input$filevcf$datapath)
      },
      error = function(e){
        stop('Upload VCF file') #safeError(e)
      }
    )

    tryCatch(
      {
        PopData = read.table(input$filetxt$datapath, header = T)
      },
      error = function(e){
        stop('Upload TXT file')
      }
    )
    require("adegenet")
    Genepool <- as.character(PopData$Genepool)
    LimaBeanData2 <- vcfR::vcfR2genind(LimaBeanGBS, pop= Genepool, NA.char= "NA")
    diversity <- summary(LimaBeanData2)
    plot(diversity$Hobs, diversity$Hexp, xlab="Hobs", ylab="Hexp", main="Hexp as a function of Hobs per locus")
  })

  ############################################## BDSpc

  data1 <- reactive({
    tryCatch(
      {
        LimaBeanGBS <- vcfR::read.vcfR(input$filevcf$datapath)
      },
      error = function(e){
        stop("Upload VCF file")
      }
    )

    tryCatch(
      {
        data = read.table(input$filetxt$datapath, header = TRUE)
      },
      error = function(e){
        stop("Upload TXT file")
      }
    )
    Genepool <- as.character(data$Genepool)
    LimaBeanData3 <- vcfR::vcfR2genlight(LimaBeanGBS)
    adegenet::ploidy(LimaBeanData3) <- 2
    adegenet::pop(LimaBeanData3) <- Genepool

    LimaBeanData2 <- vcfR::vcfR2genind(LimaBeanGBS,
                                       pop= Genepool,
                                       NA.char= "NA")
  })

  output$summary1 <- renderPrint({

    LimaBeanData2 = data1()

    diversity.clusters <- hierfstat::basic.stats(LimaBeanData2, diploid=TRUE, digits=2)
    print(summary(diversity.clusters$Hs))
  })

  output$summary2 <- renderPrint({
    LimaBeanData2 = data1()

    diversity.clusters <- hierfstat::basic.stats(LimaBeanData2, diploid=TRUE, digits=2)
    print(summary(diversity.clusters$Ho))
  })

  output$summary3 <- renderPrint({
    LimaBeanData2 = data1()

    diversity.clusters <- hierfstat::basic.stats(LimaBeanData2, diploid=TRUE, digits=2)
    print(summary(diversity.clusters$Fis))
  })

  output$boxplot <- renderPlot({
    LimaBeanData2 = data1()

    diversity.clusters <- hierfstat::basic.stats(LimaBeanData2, diploid=TRUE, digits=2)
    boxplot(diversity.clusters$Hs, ylab="Hs")
  })

  ################################################## Results ML

  data3 <- reactive({
    tryCatch(
      {
        LimaBeanGBS = vcfR::read.vcfR(input$filevcf$datapath)
      },
      error = function(e){
        stop('Upload VCF file')
      }
    )

    tryCatch(
      {
        data = read.table(input$filetxt$datapath, header = TRUE)
      },
      error = function(e){
        stop('Upload TXT file')
      }
    )
    Genepool <- as.character(data$Genepool)
    LimaBeanData3 <- vcfR::vcfR2genlight(LimaBeanGBS)
    adegenet::ploidy(LimaBeanData3) <- 2
    adegenet::pop(LimaBeanData3) <- Genepool

    LimaBeanData2 <- vcfR::vcfR2genind(LimaBeanGBS,
                                       pop= Genepool,
                                       NA.char= "NA")
  })

  data2 <- reactive({
    tryCatch(
      {
        LimaBeanGBS = vcfR::read.vcfR(input$filevcf$datapath)
      },
      error = function(e){
        stop('Upload VCF file')
      }
    )

    tryCatch(
      {
        data = read.table(input$filetxt$datapath, header = TRUE)
      },
      error = function(e){
        stop('Upload TXT file')
      }
    )
    Genepool <- as.character(data$Genepool)
    LimaBeanData3 <- vcfR::vcfR2genlight(LimaBeanGBS)
    adegenet::ploidy(LimaBeanData3) <- 2
    adegenet::pop(LimaBeanData3) <- Genepool
    LimaBeanPCA <- adegenet::glPca(LimaBeanData3, nf = 3)
    Limapcascores <- as.data.frame(LimaBeanPCA$scores)
    Limapcascores$pop <- adegenet::pop(LimaBeanData3)
    return(list(Limapcascores, LimaBeanPCA))
  })

  output$boxplot1 <- renderPlot({

    LimaBeanData2 <- data3()

    diversity.clusters <- hierfstat::basic.stats(LimaBeanData2, diploid=TRUE, digits=2)
    boxplot(diversity.clusters$perloc[,2:3])
  })

  output$summary4 <- renderPrint({
    LimaBeanData2 <- data3()
    global.Fst.weir_cock <- hierfstat::wc(LimaBeanData2)
    print(global.Fst.weir_cock)
  })

  output$summary5 <- renderPrint({
    LimaBeanData2 <- data3()
    pairwise.fst.genepools <- hierfstat::genet.dist(LimaBeanData2, method="WC84")
    print(pairwise.fst.genepools)
  })

  output$plot1 <- renderPlot({
    Limapcascores <- data2()[[1]]
    factoextra::fviz_nbclust(as.data.frame(Limapcascores[,-4]), FUNcluster = kmeans) +
      ggplot2::theme_minimal()
  })

  output$plot2 <- renderPlot({
    Limapcascores <- data2()[[1]]
    LimaBeanPCA <- data2()[[2]]
    PVCA <-as.data.frame(Limapcascores[,1:3])
    set.seed(2022)
    putput <- kmeans(x = as.matrix.data.frame(Limapcascores[,-4]), centers = 6, nstart = 4)
    CLUSTER <- as.data.frame(putput$cluster)
    colnames(CLUSTER) <- "CLUSTER"

    final_cluster_data <- cbind(PVCA$PC1, PVCA$PC2, CLUSTER)
    rownames(final_cluster_data) <- rownames(Limapcascores)
    final_cluster_data <- as.data.frame(final_cluster_data)
    final_cluster_data$CLUSTER <- as.factor(final_cluster_data$CLUSTER)
    colnames(final_cluster_data) <- c("V1","V2","CLUSTER")

    my_pal <- c("darkgreen","darkblue", "orangered","darkred","lightslateblue","orange","purple4","darkred","green","red","pink","yellow","black","deeppink4","darkturquoise","khaki3")
    my_fill <- c("darkgreen","darkblue", "orangered","darkred","lightslateblue","orange","purple4","darkred","green","red","pink","yellow","black","deeppink4","darkturquoise","khaki3")
    p1 <- ggplot2::ggplot(final_cluster_data, ggplot2::aes(x = V1, y = V2,color = CLUSTER))
    p1 <- p1 + ggplot2::geom_point(size = 3, ggplot2::aes(fill = CLUSTER),alpha =0.5)
    p1 <- p1 + ggplot2::geom_hline(yintercept = 0)
    p1 <- p1 + ggplot2::geom_vline(xintercept = 0)
    p1 <- p1 + ggplot2::geom_point(size=3)
    p1 <- p1 + ggplot2::theme_bw()
    p1 <- p1 + ggplot2::scale_color_manual(values=c(my_pal))
    p1 <- p1 + ggplot2::scale_fill_manual(values=c(paste(my_fill)))+
      ggplot2::xlab(sprintf("PC1 %f percent", 100*LimaBeanPCA$eig[1]/sum(LimaBeanPCA$eig))) +
      ggplot2::ylab(sprintf("PC2 %f percent", 100*LimaBeanPCA$eig[2]/sum(LimaBeanPCA$eig)))
    p1
  })

  output$plot3 <- renderPlot({
    tryCatch(
      {
        LimaBeanGBS = vcfR::read.vcfR(input$filevcf$datapath)
      },
      error = function(e){
        stop('Upload VCF file')
      }
    )

    tryCatch(
      {
        data = read.table(input$filetxt$datapath, header = TRUE)
      },
      error = function(e){
        stop('Upload TXT file')
      }
    )
    Genepool <- as.character(data$Genepool)
    LimaBeanData3 <- vcfR::vcfR2genlight(LimaBeanGBS)
    adegenet::ploidy(LimaBeanData3) <- 2
    adegenet::pop(LimaBeanData3) <- Genepool
    Limapcascores <- data2()[[1]]
    PVCA <-as.data.frame(Limapcascores[,1:3])
    set.seed(2022)
    putput <- kmeans(x = as.matrix.data.frame(Limapcascores[,-4]), centers = 6, nstart = 4)
    CLUSTER <- as.data.frame(putput$cluster)
    colnames(CLUSTER) <- "CLUSTER"

    final_cluster_data <- cbind(PVCA$PC1, PVCA$PC2, CLUSTER)
    rownames(final_cluster_data) <- rownames(Limapcascores)
    final_cluster_data <- as.data.frame(final_cluster_data)
    final_cluster_data$CLUSTER <- as.factor(final_cluster_data$CLUSTER)
    colnames(final_cluster_data) <- c("V1","V2","CLUSTER")
    # metodo distancia
    my_pal <- c("darkgreen","darkblue", "orangered","darkred","lightslateblue","orange","purple4","darkred","green","red","pink","yellow","black","deeppink4","darkturquoise","khaki3")
    adegenet::pop(LimaBeanData3) <- final_cluster_data$CLUSTER
    tree <- poppr::aboot(LimaBeanData3, tree = "upgma", distance = nei.dist, sample = 100, showtree = F, cutoff = 50, quiet = T)
    ape::plot.phylo(tree, cex = 0.3, font = 2, adj = 0, tip.color =  my_pal[adegenet::pop(LimaBeanData3)])
    ape::nodelabels(tree$node.label, adj = c(1.3, -0.5), frame = "n", cex = 0.3,font = 3, xpd = TRUE)
    axis(side = 1)
    title(xlab = "Genetic distance (proportion of loci that are different)")
  })

  model <- reactive({
    tryCatch(
      {
        # Import data from module import_data
        racesR <- read.csv(input$filecsv$datapath)
      },
      error = function(e){
        stop('Upload CSV file')
      }
    )

    ## Name data ####

    data <- racesR

    ### Convert double data to factor ###
    indx <- sapply(data, is.double)
    data[indx] <- lapply(data[indx], function(data) as.factor(as.double(data)))

    ### Convert factor data to double ###

    data$P100S <- as.numeric(data$P100S)
    data$DAF <- as.numeric(data$DAF)
    data$DAM <- as.numeric(data$DAM)
    data$LONG_SEM <- as.numeric(data$LONG_SEM)
    data$ANCHO_SEM <- as.numeric(data$ANCHO_SEM)
    data$GRUESO_SEM <- as.numeric(data$GRUESO_SEM)
    data$Races <- as.factor(data$Races)

    ######### ML ######################

    Datos <- as.data.frame(data[,-1])
    missingcols <- sapply(Datos, function(x) { any(is.na(x)) })
    tcontrol <- caret::trainControl(method="repeatedcv", number=50, repeats=3)
    set.seed(123)

    # replace data by keeping only those variables that don't have missing data
    data <- Datos[ , !missingcols]

    data %>% dplyr::group_by(Races) %>% dplyr::summarise(n = dplyr::n())

    # create training and test sets
    inTrain <- caret::createDataPartition(y = data$Races, p = 0.6, list = FALSE)

    # subset
    training <- data[inTrain, ]
    testing <- data[-inTrain, ]

    training$Races <- as.factor(training$Races)
    testing$Races <- as.factor(testing$Races)

    metric <- "Accuracy"
    mtry <- sqrt(ncol(training))
    modFitN <- caret::train(Races~., method = "rf",
                            data = training[], trControl = tcontrol,
                            metric=metric, tuneLength=15)

    modFit.rfN <- randomForest::randomForest(Races ~., data = training[],
                                             mtry= modFitN$bestTune$mtry)

    return(list(modFitN, modFit.rfN))
  })

  output$modelF <- renderPrint({
    modFitN <- model()[[1]]

    print(modFitN)

  })

  output$modFit <- renderPlot({

    modFitN <- model()[[1]]

    plot(modFitN)

  })

  output$con_matrix <- renderPlot({

    modFit.rfN <- model()[[2]]

    ### Confusion matrix ###

    TablaE <- as.table(modFit.rfN$confusion)
    TablaE <- as.data.frame(TablaE)

    TablaE <- TablaE[1:36, ]

    ### Training confusion ###

    ggplot2::ggplot(TablaE, aes(x=Var1, y=Var2, fill=Freq)) +
      ggplot2::geom_tile(color="black") +
      ggplot2::theme_bw() +
      ggplot2::coord_equal() +
      ggplot2::scale_fill_distiller(palette="Greys", direction=1) +
      ggplot2::guides(fill=F) + # removing legend for `fill`
      ggplot2::labs(title = "Training Confusion Matrix") + # using a title instead
      ggplot2::geom_text(aes(label=Freq), color="black") +
      ggplot2::theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.x=element_blank(), axis.title.y=element_blank())

  })

  model1 <- reactive({
    tryCatch(
      {
        # Import data from module import_data
        racesR <- read.csv(input$filecsv1$datapath)
      },
      error = function(e){
        stop('Upload CSV file')
      }
    )

    ## Name data ####

    data <- racesR

    ### Convert double data to factor ###
    indx <- sapply(data, is.double)
    data[indx] <- lapply(data[indx], function(data) as.factor(as.double(data)))

    ### Convert factor data to double ###

    data$P100S <- as.numeric(data$P100S)
    data$DAF <- as.numeric(data$DAF)
    data$DAM <- as.numeric(data$DAM)
    data$LONG_SEM <- as.numeric(data$LONG_SEM)
    data$ANCHO_SEM <- as.numeric(data$ANCHO_SEM)
    data$GRUESO_SEM <- as.numeric(data$GRUESO_SEM)
    data$Races <- as.factor(data$Races)

    ######### ML ######################

    Datos <- as.data.frame(data[,-1])
    missingcols <- sapply(Datos, function(x) { any(is.na(x)) })
    tcontrol <- caret::trainControl(method="repeatedcv", number=50, repeats=3)
    set.seed(123)

    # replace data by keeping only those variables that don't have missing data
    data <- Datos[ , !missingcols]

    data %>% dplyr::group_by(Races) %>% summarise(n = n())

    # create training and test sets
    inTrain <- caret::createDataPartition(y = data$Races, p = 0.6, list = FALSE)

    # subset
    training <- data[inTrain, ]
    testing <- data[-inTrain, ]

    training$Races <- as.factor(training$Races)
    testing$Races <- as.factor(testing$Races)

    metric <- "Accuracy"
    mtry <- sqrt(ncol(training))
    modFitN <- caret::train(Races~., method = "rf",
                            data = training[], trControl = tcontrol,
                            metric=metric, tuneLength=15)

    predictions.rfN<- stats::predict(modFitN, newdata = testing[])

    Datos$Races<- as.factor(Datos$Races)

    modFit.rf <- randomForest::randomForest(Races ~ ., data = Datos, localImp = TRUE, mtry=modFitN$bestTune$mtry)

    VI <- varImp(modFit.rf, scale=T)
    VI$Variable <- row.names(VI)
    return(modFit.rf)
  })

  output$min_depth <- renderPlot({
    modFit.rf <- model1()
    min_depth_frame <- randomForestExplainer::min_depth_distribution(modFit.rf)
    randomForestExplainer::plot_min_depth_distribution(min_depth_frame, mean_sample = "top_trees")
  })

  output$multi <- renderPlot({
    modFit.rf <- model1()
    randomForestExplainer::plot_multi_way_importance(modFit.rf, size_measure = "no_of_nodes")
  })

  output$inter <- renderPlot({
    modFit.rf <- model1()
    importance_frame <- randomForestExplainer::measure_importance(modFit.rf,  measures = NULL)
    vars <- randomForestExplainer::important_variables(importance_frame, k = 10, measures = c("mean_min_depth", "no_of_trees"))
    interactions_frame <- randomForestExplainer::min_depth_interactions(modFit.rf, vars)
    randomForestExplainer::plot_min_depth_interactions(interactions_frame)
  })

  output$plot_dataS <- renderPlot({

    tryCatch(
      {
        Specie = readxl::read_xlsx(input$filexlsx$datapath)
      },
      error = function(e){
        stop('Upload Excel file')
      }
    )

    # Deleting duplicates
    recordsSpecie <- unique(Specie)
    # Deleting occurrence points with NA
    recordsSpecie <- na.omit(Specie)
    # Deleting ID column
    recordsSpecie <- recordsSpecie[,-1]

    # We transform into sf format
    Lunatus <- recordsSpecie %>% sf::st_as_sf(coords = c(1, 2),
                                             crs = "+proj=longlat +ellps=WGS84
                                      +datum=WGS84 +no_defs +towgs84=0,0,0")
    # Generate a buffer
    Hull <- Lunatus %>%
      sf::st_union() %>%
      sf::st_convex_hull()
    Buffer <- Hull %>%
      sf::st_buffer(dist = 1) %>%
      sf::st_as_sf()

    Bioclimatic <- raster::getData("worldclim", res = 2.5, var = "bio", path = tempdir())

    # Crop layers using the buffer
    Bioclimatic <- Bioclimatic %>%
      raster::crop(Buffer) %>%
      raster::trim()

    # Selección del número de background points
    Number_background_points = 1000
    # Run ENMevaluate
    Results <- ENMeval::ENMevaluate(occs =  recordsSpecie, envs = Bioclimatic,
                                    n.bg = Number_background_points,
                                    algorithm = 'maxnet', partitions = 'block',
                                    tune.args = list(fc = c("L","LQ"), rm = 1:2)) #,"LQH","H"

    # Modeling results
    Results@results


    ## Best Model Prediction
    Models <- Results@results
    Models$ID <- 1:nrow(Models)
    Models <- Models %>%
      dplyr::arrange(AICc)
    BestModels <- Results@models[[Models$ID[1]]]
    Prediction <- raster::predict(Bioclimatic, BestModels, type = "cloglog")

    plot(BestModels, c("bio1","bio2","bio3","bio4","bio5","bio6","bio7","bio8","bio9","bio10",
                       "bio11","bio12"), type = "cloglog")

  })

  output$plot_data2S <- renderPlot({
    tryCatch(
      {
        Specie = readxl::read_xlsx(input$filexlsx$datapath)
      },
      error = function(e){
        stop('Upload Excel file')
      }
    )

    # Deleting duplicates
    recordsSpecie <- unique(Specie)
    # Deleting occurrence points with NA
    recordsSpecie <- na.omit(Specie)
    # Deleting ID column
    recordsSpecie <- recordsSpecie[,-1]

    # We transform into sf format
    Lunatus <- recordsSpecie %>% sf::st_as_sf(coords = c(1, 2),
                                             crs = "+proj=longlat +ellps=WGS84
                                      +datum=WGS84 +no_defs +towgs84=0,0,0")
    # Generate a buffer
    Hull <- Lunatus %>%
      sf::st_union() %>%
      sf::st_convex_hull()
    Buffer <- Hull %>%
      sf::st_buffer(dist = 1) %>%
      sf::st_as_sf()

    Bioclimatic <- raster::getData("worldclim", res = 2.5, var = "bio", path = tempdir())

    # Crop layers using the buffer
    Bioclimatic <- Bioclimatic %>%
      raster::crop(Buffer) %>%
      raster::trim()

    # Selección del número de background points
    Number_background_points = 1000
    # Run ENMevaluate
    Results <- ENMeval::ENMevaluate(occs =  recordsSpecie, envs = Bioclimatic,
                                    n.bg = Number_background_points,
                                    algorithm = 'maxnet', partitions = 'block',
                                    tune.args = list(fc = c("L","LQ"), rm = 1:2)) #,"LQH","H"
    # Modeling results
    Results@results

    ## Best Model Prediction
    Models <- Results@results
    Models$ID <- 1:nrow(Models)
    Models <- Models %>%
      dplyr::arrange(AICc)
    BestModels <- Results@models[[Models$ID[1]]]
    Prediction <- raster::predict(Bioclimatic, BestModels, type = "cloglog")

    gdd_ggplot2 <- Prediction$layer %>% as("SpatialPixelsDataFrame") %>%
      as.data.frame()

    # Plotting using ggplot2
    plot_1 <- ggplot2::ggplot() +
      ggplot2::geom_tile(data = gdd_ggplot2,
                         ggplot2::aes(x = x, y = y, fill = layer)) +
      viridis::scale_fill_viridis(option="turbo") +
      ggplot2::xlab("Longitude") +
      ggplot2::ylab("Latitude")+
      ggplot2::theme(panel.grid.major = ggplot2::element_line(size = 0.5,
                                                              linetype = 'solid',
                                                              colour = "gray14"),
                     panel.grid.minor = ggplot2::element_line(size = 0.25,
                                                              linetype = 'solid',
                                                              colour = "gray14"),
                     panel.background = ggplot2::element_rect(fill = "black",
                                                              colour = "black",
                                                              size = 0.5,
                                                              linetype = "solid"))+
      ggplot2::labs(fill = "Probability") +
      ggplot2::labs(title= "Spatial Distribution Modeling") +
      ggplot2::labs(subtitle = "MaxEnt By ENMeval (maxnet)",
                    caption = " Maxent Algorithm (Phillips et al. 2006, Phillips et al. 2008, Phillips et al. 2017)\n using the R-package ENMeval (Kass et al. 2021, https://doi.org/10.1111/2041-210X.13628)")+
      ggplot2::theme(plot.caption.position = "plot",
                     plot.caption = ggplot2::element_text(hjust = 0)) +
      ggplot2::labs(fill = "Probability")+
      ggspatial::annotation_north_arrow(location = "tr",
                                        which_north = "true",
                                        pad_x = ggplot2::unit(0.4, "in"),
                                        pad_y = ggplot2::unit(0.4, "in"),
                                        style = ggspatial::north_arrow_nautical(fill = c("grey40", "gray49"),
                                                                                line_col = "gray30"))
    plot_1
  })

  output$plot_data3S <- renderPlot({
    tryCatch(
      {
        Specie = readxl::read_xlsx(input$filexlsx$datapath)
      },
      error = function(e){
        stop('Upload Excel file')
      }
    )

    # Deleting duplicates
    recordsSpecie <- unique(Specie)
    # Deleting occurrence points with NA
    recordsSpecie <- na.omit(Specie)
    # Deleting ID column
    recordsSpecie <- recordsSpecie[,-1]

    # We transform into sf format
    Lunatus <- recordsSpecie %>% sf::st_as_sf(coords = c(1, 2),
                                             crs = "+proj=longlat +ellps=WGS84
                                      +datum=WGS84 +no_defs +towgs84=0,0,0")
    # Generate a buffer
    Hull <- Lunatus %>%
      sf::st_union() %>%
      sf::st_convex_hull()
    Buffer <- Hull %>%
      sf::st_buffer(dist = 1) %>%
      sf::st_as_sf()

    Bioclimatic <- raster::getData("worldclim", res = 2.5, var = "bio", path = tempdir())

    # Crop layers using the buffer
    Bioclimatic <- Bioclimatic %>%
      raster::crop(Buffer) %>%
      raster::trim()

    # Option 1: Extracting Environmental Data from Occurrences Points of Species
    Occurrences<-raster::extract(Bioclimatic, recordsSpecie)
    # Removing Missing Data
    Occurrences<-na.omit(Occurrences)
    # Removing Repeating Data from Points
    Occurrences<-unique(Occurrences)

    # Selecting pseudo-absence Points
    background_data <- raster::sampleRandom(Bioclimatic,length(recordsSpecie$Longitude),cells=T)
    # Removing Repeating Data from Layer Cells
    background_data <- unique(background_data)
    # Removing the Number of Cells Column
    background_data <- background_data[,-1]
    # Merging Data
    Environmental_values <- data.frame(rbind(Occurrences,background_data))
    # Labeling the presences as '1' and background as '0'
    y <- c(rep(1,nrow(Occurrences)), rep(0,nrow(background_data)))
    Data_01 <- cbind(y, Environmental_values)
    # Table as Data.Frame
    Data_SDM <- as.data.frame(Data_01)
    Data_SDM <- na.omit(Data_SDM)

    # Set Seed
    set.seed(2021)
    training <- caret::createDataPartition(y = Data_SDM$y, p = 0.8, list = FALSE, times = 1)

    # Training
    Data_training <- Data_SDM[training, ]
    Data_training_modeling <- Data_training
    Data_training_modeling$y <- as.factor(Data_training_modeling$y)
    levels(Data_training_modeling$y) <- c("back", "pres")
    # Training Data Ratio
    prop.table(table(Data_training$y))
    #Testing
    Data_test  <- Data_SDM[-training, ]
    Data_test_modeling <- Data_test
    Data_test_modeling$y <- as.factor(Data_test_modeling$y)
    levels(Data_test_modeling$y) <- c("back", "pres")
    # Test Data Ratio
    prop.table(table(Data_test$y))

    # Definition of partitions and repetitions
    partitions  <- 5
    repetitions <- 3

    # Definition of Training
    control_tr <- caret::trainControl(method = "repeatedcv", number = partitions,
                                      repeats = repetitions,
                                      returnResamp = "final", verboseIter = FALSE,
                                      classProbs = TRUE)
    # Model fit
    SDM_glm_model <- caret::train(y ~ ., data = Data_training_modeling,
                                  method = "glm",
                                  preProcess = c('center', 'scale'),
                                  trControl = control_tr,
                                  metric = 'Accuracy')
    # Model output
    #SDM_glm_model

    # Testing
    Prediction_raw <- raster::predict(SDM_glm_model, newdata = Data_test_modeling)
    # Confusio nMatrix
    confusionMatrix <- caret::confusionMatrix(data = Prediction_raw, reference = Data_test_modeling$y,positive = "pres")
    #confusionMatrix

    # Error Test
    error_test <- mean(Prediction_raw != Data_test_modeling$y)

    # Extrapolating Model to others Zones
    mapGLM <- raster::predict(Bioclimatic,SDM_glm_model, progress="text",type = "prob",args=c("extrapolate=T", "doclamp=TRUE"))
    mapGLM <- 1-mapGLM

    # Projections to Data.Frame Format
    mapGLM.ggplot2 <- mapGLM %>%
      as("SpatialPixelsDataFrame") %>%
      as.data.frame()

    # Plotting using ggplot2
    plot_2 <- ggplot2::ggplot() +
      ggplot2::geom_tile(data = mapGLM.ggplot2,
                         ggplot2::aes(x = x, y = y, fill = layer)) +
      viridis::scale_fill_viridis(option="turbo") +
      ggplot2::xlab("Longitude") +
      ggplot2::ylab("Latitude")+
      # annotation_north_arrow(location = "tr", which_north = "true", style = north_arrow_nautical) +
      ggplot2::theme(panel.grid.major = ggplot2::element_line(size = 0.5,
                                                              linetype = 'solid',
                                                              colour = "gray14"),
                     panel.grid.minor = ggplot2::element_line(size = 0.25,
                                                              linetype = 'solid',
                                                              colour = "gray14"),
                     panel.background = ggplot2::element_rect(fill = "black",
                                                              colour = "black",
                                                              size = 0.5,
                                                              linetype = "solid"))+
      ggplot2::labs(fill = "Probability") +
      ggplot2::labs(title= "Spatial Distribution Modeling") +
      ggplot2::labs(subtitle = "Generalized Linear Models By glm() function",
                    caption = " glm() function from the basic R-package using caret (https://www.cienciadedatos.net/machine-learning-r.h)")+
      ggplot2::theme(plot.caption.position = "plot",
                     plot.caption = ggplot2::element_text(hjust = 0)) +
      ggplot2::labs(fill = "Probability") +
      ggspatial::annotation_north_arrow(location = "tr",
                                        which_north = "true",
                                        pad_x = ggplot2::unit(0.4, "in"),
                                        pad_y = ggplot2::unit(0.4, "in"),
                                        style = ggspatial::north_arrow_nautical(fill = c("grey40", "gray49"),
                                                                                line_col = "gray30"))
    plot_2
  })


  values <- reactiveValues(ploting = NULL, plott = NULL)

  observeEvent(input$gob,{
    tryCatch(
      {
        Lima = readxl::read_excel(input$filexlsxG$datapath)
      },
      error = function(e){
        stop('Upload VCF file')
      }
    )

    tryCatch(
      {
        data = input$filetxtG$datapath
      },
      error = function(e){
        stop('Upload TXT file')
      }
    )

    myY = as.data.frame(Lima)
    myG <- as.data.frame(readr::read_delim(data,
                                           delim = "\t", escape_double = FALSE,
                                           col_names = FALSE, trim_ws = TRUE))

    PCA_number <- 8
    kinship <- "Vanraden"
    kinship_cluster <- "ward"
    MAF <- 0.05
    defaultW <- getOption("warn")
    dirname <-  '~/GAPIT/MLM'
    if (!dir.exists(dirname))dir.create(dirname,recursive=TRUE)

    # set working folder path for MLM
    setwd("~/GAPIT/MLM")
    # Run the MLM model
    myGAPIT_MLM <- GAPIT3::GAPIT(
      Y = myY,
      G = myG,
      model = "MLM",
      PCA.total = PCA_number,
      file.output = TRUE,
      Geno.View.output = FALSE,
      kinship.cluster = kinship_cluster,
      SNP.MAF = MAF
    )
    options(warn = defaultW)

    dirname <-  '~/GAPIT/Farm-CPU'
    if (!dir.exists(dirname))dir.create(dirname,recursive=TRUE)

    # set working folder path for Farm-CPU
    setwd('~/GAPIT/Farm-CPU')
    # Run the Farm-CPU model
    myGAPIT_MLM <- GAPIT3::GAPIT(
      Y=myY,
      G=myG,
      model="Farm-CPU",
      PCA.total=PCA_number,
      file.output=T,Geno.View.output=F,
      # kinship.algorithm=kinship,
      kinship.cluster=kinship_cluster,
      SNP.MAF=MAF
    )
    options(warn = defaultW)

    dirname <-  '~/GAPIT/Blink'
    if (!dir.exists(dirname))dir.create(dirname,recursive=TRUE)

    # set working folder path for Blink
    setwd('~/GAPIT/Blink')
    # Run the Blink model
    myGAPIT_MLM <- GAPIT3::GAPIT(
      Y=myY,
      G=myG,
      model="Blink",
      PCA.total=PCA_number,
      file.output=T,Geno.View.output=F,
      # kinship.algorithm=kinship,
      kinship.cluster=kinship_cluster,
      SNP.MAF=MAF
    )
    options(warn = defaultW)
    values$ploting <- 1

  })

  output$ready <- renderPrint({
    if (values$ploting == 1){
      print('Please click in See plot!')
    }else if(is.null(values$ploting)){
      print('Please wait!')
    }
  })

  observeEvent(input$Manhattan, {
    if(values$ploting == 1){
      values$plott = 1
    }
  })

  output$plot_data1G <- renderPlot({
    if(values$plott == 1){
      image <- magick::image_read_pdf("~/GAPIT/MLM/GAPIT.Manhattan.Multiple.Plot.circular.pdf")
      file_i <- magick::image_scale(image, geometry = "700")
      plot(file_i)
    }
    else{
      file_i <- NULL
    }
  })

  output$plot_data2G <- renderPlot({
    if(values$plott == 1){
      image = magick::image_read_pdf("~/GAPIT/MLM/GAPIT.Multiple.QQ.plot.symphysic.pdf")
      file_i <- magick::image_scale(image, geometry = "700")
      plot(file_i)
    }else{file_i <- NULL}
  })

  output$plot_data3G <- renderPlot({
    if(values$plott == 1){
      image <- magick::image_read_pdf("~/GAPIT/Farm-CPU/GAPIT.Manhattan.Multiple.Plot.circular.pdf")
      file_i <- magick::image_scale(image, geometry = "700")
      plot(file_i)
    }
    else{
      file_i <- NULL
    }
  })

  output$plot_data4G <- renderPlot({
    if(values$plott == 1){
      image = magick::image_read_pdf("~/GAPIT/Farm-CPU/GAPIT.Multiple.QQ.plot.symphysic.pdf")
      file_i <- magick::image_scale(image, geometry = "700")
      plot(file_i)
    }else{file_i <- NULL}
  })

  output$plot_data5G <- renderPlot({
    if(values$plott == 1){
      image <- magick::image_read_pdf("~/GAPIT/Blink/GAPIT.Manhattan.Multiple.Plot.circular.pdf")
      file_i <- magick::image_scale(image, geometry = "700")
      plot(file_i)
    }
    else{
      file_i <- NULL
    }
  })

  output$plot_data6G <- renderPlot({
    if(values$plott == 1){
      image = magick::image_read_pdf("~/GAPIT/Blink/GAPIT.Multiple.QQ.plot.symphysic.pdf")
      file_i <- magick::image_scale(image, geometry = "700")
      plot(file_i)
    }else{file_i <- NULL}
  })


  data <- reactive({

    tryCatch(
      {
        Y1 <- read.table(input$filelfmmL$datapath, quote="\"", comment.char="", na.strings="9")
      },
      error = function(e){
        stop('Upload LFMM file')
      }
    )

    pc1 <- prcomp(Y1) # carry out a PCA

  })

  output$plot_dataL <- renderPlot({

    pc1 <- data()

    plot(pc1$sdev[1:15]^2, ylab = "percentage of variance explained", xlab = 'Axes', ) #plot results

  })

  output$plot_data1L <- renderPlot({
    pc1 <- data()

    screeplot(pc1, main = "Screeplot of Genetic Data with Broken Stick", bstick=TRUE, type="barplot")
  })

  output$plot_data2L <- renderPlot({
    tryCatch(
      {
        data1 <- LEA::vcf2lfmm(input$filevcfL$datapath, force = TRUE)
      },
      error = function(e){
        stop('Upload VCF file')
      }
    )

    tryCatch(
      {
        Y1 <- read.table(input$filelfmmL$datapath, quote="\"", comment.char="", na.strings="9")
      },
      error = function(e){
        stop('Upload LFMM file')
      }
    )

    tryCatch(
      {
        bio10 <- read.table(input$filetxtL$datapath, quote="\"", comment.char="")
      },
      error = function(e){
        stop('Upload TXT file')
      }
    )

    mod.lfmm <- lfmm::lfmm_ridge(Y = Y1,
                                 X = bio10,
                                 K = 8)

    pv <- lfmm::lfmm_test(Y = Y1,
                          X = bio10,
                          lfmm = mod.lfmm,
                          calibrate = "gif")
    hist(pv$pvalue[,1], main="Unadjusted p-values", xlab = 'P value')
  })

  output$plot_data3L <- renderPlot({
    tryCatch(
      {
        data1 <- LEA::vcf2lfmm(input$filevcfL$datapath, force = TRUE)
      },
      error = function(e){
        stop('Upload VCF file')
      }
    )

    tryCatch(
      {
        Y1 <- read.table(input$filelfmmL$datapath, quote="\"", comment.char="", na.strings="9")
      },
      error = function(e){
        stop('Upload LFMM file')
      }
    )

    tryCatch(
      {
        bio10 <- read.table(input$filetxtL$datapath, quote="\"", comment.char="")
      },
      error = function(e){
        stop('Upload TXT file')
      }
    )

    mod.lfmm <- lfmm::lfmm_ridge(Y = Y1,
                                 X = bio10,
                                 K = 8)

    pv <- lfmm::lfmm_test(Y = Y1,
                          X = bio10,
                          lfmm = mod.lfmm,
                          calibrate = "gif")

    hist(pv$calibrated.pvalue[,1], main="GIF-adjusted p-values", xlab = 'P value calibrated')
  })

  output$plot_data4L <- renderPlot({
    tryCatch(
      {
        data1 <- LEA::vcf2lfmm(input$filevcfL$datapath, force = TRUE)
      },
      error = function(e){
        stop('Upload VCF file')
      }
    )

    tryCatch(
      {
        Y1 <- read.table(input$filelfmmL$datapath, quote="\"", comment.char="", na.strings="9")
      },
      error = function(e){
        stop('Upload LFMM file')
      }
    )

    tryCatch(
      {
        bio10 <- read.table(input$filetxtL$datapath, quote="\"", comment.char="")
      },
      error = function(e){
        stop('Upload TXT file')
      }
    )
    mod.lfmm <- lfmm::lfmm_ridge(Y = Y1,
                                 X = bio10,
                                 K = 8)

    pv <- lfmm::lfmm_test(Y = Y1,
                          X = bio10,
                          lfmm = mod.lfmm,
                          calibrate = "gif")

    pvaluesUncal <- pv$pvalue
    qqplot(rexp(length(pvaluesUncal), rate = log(10)),
           -log10(pvaluesUncal),
           xlab = "Expected quantile",
           pch = 19, cex = .4)
    abline(0,1)
  })

  output$plot_data5L <- renderPlot({
    tryCatch(
      {
        data1 <- LEA::vcf2lfmm(input$filevcfL$datapath, force = TRUE)
      },
      error = function(e){
        stop('Upload VCF file')
      }
    )

    tryCatch(
      {
        Y1 <- read.table(input$filelfmmL$datapath, quote="\"", comment.char="", na.strings="9")
      },
      error = function(e){
        stop('Upload LFMM file')
      }
    )

    tryCatch(
      {
        bio10 <- read.table(input$filetxtL$datapath, quote="\"", comment.char="")
      },
      error = function(e){
        stop('Upload TXT file')
      }
    )
    mod.lfmm <- lfmm::lfmm_ridge(Y = Y1,
                                 X = bio10,
                                 K = 8)
    pv <- lfmm::lfmm_test(Y = Y1,
                          X = bio10,
                          lfmm = mod.lfmm,
                          calibrate = "gif")

    pvalues <- pv$calibrated.pvalue
    qqplot(rexp(length(pvalues), rate = log(10)),
           -log10(pvalues), xlab = "Expected quantile",
           pch = 19, cex = .4)
    abline(0,1)
  })

  output$table1L <- DT::renderDataTable({
    tryCatch(
      {
        bio10 <- read.table(input$filetxtL$datapath, quote="\"", comment.char="")
      },
      error = function(e){
        stop('Upload TXT file')
      }
    )

    tryCatch(
      {
        Y1 <- read.table(input$filelfmmL$datapath, quote="\"", comment.char="", na.strings="9")
      },
      error = function(e){
        stop('Upload LFMM file')
      }
    )

    tryCatch(
      {
        chr_pos_10668snps <- read.delim(input$filetxt1L$datapath)
      },
      error = function(e){
        stop('Upload TXT file')
      }
    )

    mod.lfmm <- lfmm::lfmm_ridge(Y = Y1,
                                 X = bio10,
                                 K = 8)

    pv <- lfmm::lfmm_test(Y = Y1,
                          X = bio10,
                          lfmm = mod.lfmm,
                          calibrate = "gif")
    qv <- qvalue::qvalue(pv$calibrated.pvalue)$qvalues

    lima_lfmm <- cbind(chr_pos_10668snps, pv$calibrated.pvalue, qv)
    colnames(lima_lfmm)[4:5] <- c("calibrated.pvalues", "pvalues")

    DT::datatable(lima_lfmm,
                  filter = 'top', extensions = c('Buttons', 'Scroller'),
                  options = list(scrollY = 650,
                                 scrollX = 500,
                                 deferRender = TRUE,
                                 scroller = TRUE,
                                 # paging = TRUE,
                                 # pageLength = 25,
                                 buttons = list('excel',
                                                list(extend = 'colvis', targets = 0, visible = FALSE)),
                                 dom = 'lBfrtip',
                                 fixedColumns = TRUE),
                  rownames = FALSE)

  })

  output$plot_data6L <- renderPlot({
    tryCatch(
      {
        bio10 <- read.table(input$filetxtL$datapath, quote="\"", comment.char="")
      },
      error = function(e){
        stop('Upload TXT file')
      }
    )

    tryCatch(
      {
        Y1 <- read.table(input$filelfmmL$datapath, quote="\"", comment.char="", na.strings="9")
      },
      error = function(e){
        stop('Upload LFMM file')
      }
    )

    tryCatch(
      {
        chr_pos_10668snps <- read.delim(input$filetxt1L$datapath)
      },
      error = function(e){
        stop('Upload TXT file')
      }
    )

    mod.lfmm <- lfmm::lfmm_ridge(Y = Y1,
                                 X = bio10,
                                 K = 8)

    pv <- lfmm::lfmm_test(Y = Y1,
                          X = bio10,
                          lfmm = mod.lfmm,
                          calibrate = "gif")
    qv <- qvalue::qvalue(pv$calibrated.pvalue)$qvalues

    lima_lfmm <- cbind(chr_pos_10668snps, pv$calibrated.pvalue, qv)
    colnames(lima_lfmm)[4:5] <- c("calibrated.pvalues", "pvalues")

    qqman::manhattan(lima_lfmm, chr="chrom", bp="bp", snp="snp",
                     p="calibrated.pvalues",suggestiveline = FALSE,
                     genomewideline = -log10(8e-05)) #genomewideline corresponds to qvalues < 0.1
  })

  output$plot_data7L <- renderPlot({
    tryCatch(
      {
        bio10 <- read.table(input$filetxtL$datapath, quote="\"", comment.char="")
      },
      error = function(e){
        stop('Upload TXT file')
      }
    )

    tryCatch(
      {
        Y1 <- read.table(input$filelfmmL$datapath, quote="\"", comment.char="", na.strings="9")
      },
      error = function(e){
        stop('Upload LFMM file')
      }
    )

    tryCatch(
      {
        chr_pos_10668snps <- read.delim(input$filetxt1L$datapath)
      },
      error = function(e){
        stop('Upload TXT file')
      }
    )

    mod.lfmm <- lfmm::lfmm_ridge(Y = Y1,
                                 X = bio10,
                                 K = 8)

    b.values <- lfmm::effect_size(Y1, bio10, mod.lfmm)
    hist(b.values, main = '')
  })

  ########### CPS

  output$plot_dataC <- renderPlot({
    tryCatch(
      {
        Y1 <- read.table(input$filelfmmC$datapath, quote="\"", comment.char="", na.strings="9")
      },
      error = function(e){
        stop('Upload LFMM file')
      }
    )
    tryCatch(
      {
        bio10 <- read.table(input$filetxt2C$datapath, quote="\"", comment.char="")
      },
      error = function(e){
        stop('Upload TXT file')
      }
    )

    tryCatch(
      {
        PopData <- read.table(input$filetxt1C$datapath, header = T)
      },
      error = function(e){
        stop('Upload TXT file')
      }
    )

    pop <- PopData$pop

    #Create dataframe with genetic data, environmental data and SNP matrix
    data<-as.data.frame(cbind(pop,bio10,Y1))

    colnames(data)[1:2]<-c("pop","bio10")

    # First, explore the overall distribution of the phenotype
    hist(data$bio10, xlab = 'Phenotype', main = '')

  })

  output$corrected <- renderPrint({

    tryCatch(
      {
        Y1 <- read.table(input$filelfmmC$datapath, quote="\"", comment.char="", na.strings="9")
      },
      error = function(e){
        stop('Upload LFMM file')
      }
    )
    tryCatch(
      {
        bio10 <- read.table(input$filetxt2C$datapath, quote="\"", comment.char="")
      },
      error = function(e){
        stop('Upload TXT file')
      }
    )

    tryCatch(
      {
        PopData <- read.table(input$filetxt1C$datapath, header = T)
      },
      error = function(e){
        stop('Upload TXT file')
      }
    )

    pop <- PopData$pop

    #Create dataframe with genetic data, environmental data and SNP matrix
    data<-as.data.frame(cbind(pop,bio10,Y1))

    data_corrected <- data

    data_corrected[,2:ncol(data_corrected)] <- NA

    for (i in 2:ncol(data)){
      LM_SNP_i <- lm(data[,i] ~ factor(data$pop))
      data_corrected[,i] <- LM_SNP_i$residuals
      colnames(data_corrected)[i]<-colnames(data)[i]
      #if(i%%50==0) print(i)
    }

    data_corrected[,ncol(data_corrected)]-LM_SNP_i$residuals
  })


  output$plot_data2C <- renderPlot({

    tryCatch(
      {
        Y1 <- read.table(input$filelfmmC$datapath, quote="\"", comment.char="", na.strings="9")
      },
      error = function(e){
        stop('Upload LFMM file')
      }
    )
    tryCatch(
      {
        bio10 <- read.table(input$filetxt2C$datapath, quote="\"", comment.char="")
      },
      error = function(e){
        stop('Upload TXT file')
      }
    )

    tryCatch(
      {
        PopData <- read.table(input$filetxt1C$datapath, header = T)
      },
      error = function(e){
        stop('Upload TXT file')
      }
    )

    pop <- PopData$pop

    #Create dataframe with genetic data, environmental data and SNP matrix
    data <- as.data.frame(cbind(pop,bio10,Y1))

    colnames(data)[1:2] <- c("pop","bio10")

    data_corrected <- data

    data_corrected[,2:ncol(data_corrected)] <- NA

    for (i in 2:ncol(data)){
      LM_SNP_i <- lm(data[,i] ~ factor(data$pop))
      data_corrected[,i] <- LM_SNP_i$residuals
      colnames(data_corrected)[i]<-colnames(data)[i]
      #if(i%%50==0) print(i)
    }

    p <- ncol(data_corrected)-2

    results_optimization <- matrix(data=NA , nrow = 0, ncol = 3)
    for (i in seq(from = 100, to = 1000 , by = 100)){
      print(i)
      for (j in c(sqrt(p), 2*sqrt(p), 0.1*p, 0.2*p, p/3, p)){
        rf_ij <- ranger::ranger(x = data_corrected[,3:ncol(data_corrected)],
                                y = data_corrected$bio10,
                                num.trees=i,
                                mtry=j)
        results_optimization <- rbind(results_optimization, c(i,j,tail(rf_ij$r.squared,1)))
      }
    }

    # Clean up the file format
    results_optimization<-as.data.frame(results_optimization)
    colnames(results_optimization)<-c("ntree", "mtry","PVE")

    #write.csv(results_optimization,file="lima_bean_results.optimization.csv", row.names = FALSE, quote=FALSE)

    # Now plot results to see if there's a plateau


    plot(results_optimization$ntree[results_optimization$mtry == sqrt(p)],results_optimization$PVE[results_optimization$mtry == sqrt(p)], type="l", col="black", xlab="ntree",ylab="PVE", ylim=range(results_optimization$PVE))
    lines(results_optimization$ntree[results_optimization$mtry == 2*sqrt(p)],results_optimization$PVE[results_optimization$mtry == 2*sqrt(p)], col="blue")
    lines(results_optimization$ntree[results_optimization$mtry == 0.1*p],results_optimization$PVE[results_optimization$mtry == 0.1*p], col="green")
    lines(results_optimization$ntree[results_optimization$mtry == 0.2*p],results_optimization$PVE[results_optimization$mtry == 0.2*p], col="purple")
    lines(results_optimization$ntree[results_optimization$mtry == p/3],results_optimization$PVE[results_optimization$mtry == p/3], col="orange")
    lines(results_optimization$ntree[results_optimization$mtry == p],results_optimization$PVE[results_optimization$mtry == p], col="red")
  })

  output$hist1 <- renderPlot({
    tryCatch(
      {
        Y1 <- read.table(input$filelfmmC$datapath, quote="\"", comment.char="", na.strings="9")
      },
      error = function(e){
        stop('Upload LFMM file')
      }
    )
    tryCatch(
      {
        bio10 <- read.table(input$filetxt2C$datapath, quote="\"", comment.char="")
      },
      error = function(e){
        stop('Upload TXT file')
      }
    )

    tryCatch(
      {
        PopData <- read.table(input$filetxt1C$datapath, header = T)
      },
      error = function(e){
        stop('Upload TXT file')
      }
    )

    tryCatch(
      {
        pos <- read.table(input$filetxt3C$datapath, header = T)
      },
      error = function(e){
        stop('Upload TXT file')
      }
    )

    pop <- PopData$pop

    #Create dataframe with genetic data, environmental data and SNP matrix
    data <- as.data.frame(cbind(pop,bio10,Y1))

    colnames(data)[1:2] <- c("pop","bio10")

    data_corrected <- data

    data_corrected[,2:ncol(data_corrected)] <- NA

    for (i in 2:ncol(data)){
      LM_SNP_i <- lm(data[,i] ~ factor(data$pop))
      data_corrected[,i] <- LM_SNP_i$residuals
      colnames(data_corrected)[i]<-colnames(data)[i]
      #if(i%%50==0) print(i)
    }

    #Run 10 random forests with 1000 trees each and the optimized mtry value (0.1)
    r2vim <- Pomona::var.sel.r2vim(x = data_corrected[,3:ncol(data_corrected)],
                                   y=data_corrected$bio10,
                                   no.runs = 10,
                                   factor = 1,
                                   ntree = 1000,
                                   mtry.prop = 0.1,
                                   method = "ranger",
                                   type = "regression")

    #calculate correlation between runs
    runs <- as.data.frame(cbind(r2vim$info$vim.run.1,
                                r2vim$info$vim.run.2,
                                r2vim$info$vim.run.3,
                                r2vim$info$vim.run.4,
                                r2vim$info$vim.run.5,
                                r2vim$info$vim.run.6,
                                r2vim$info$vim.run.7,
                                r2vim$info$vim.run.8,
                                r2vim$info$vim.run.9,
                                r2vim$info$vim.run.10))

    hist(cor(runs), col="blue", xlab = 'Correlation', main= '')

  })
  output$data_table <- DT::renderDataTable({
    tryCatch(
      {
        Y1 <- read.table(input$filelfmmC$datapath, quote="\"", comment.char="", na.strings="9")
      },
      error = function(e){
        stop('Upload LFMM file')
      }
    )
    tryCatch(
      {
        bio10 <- read.table(input$filetxt2C$datapath, quote="\"", comment.char="")
      },
      error = function(e){
        stop('Upload TXT file')
      }
    )

    tryCatch(
      {
        PopData <- read.table(input$filetxt1C$datapath, header = T)
      },
      error = function(e){
        stop('Upload TXT file')
      }
    )

    tryCatch(
      {
        pos <- read.table(input$filetxt3C$datapath, header = T)
      },
      error = function(e){
        stop('Upload TXT file')
      }
    )

    pop <- PopData$pop

    #Create dataframe with genetic data, environmental data and SNP matrix
    data <- as.data.frame(cbind(pop,bio10,Y1))

    colnames(data)[1:2] <- c("pop","bio10")

    data_corrected <- data

    data_corrected[,2:ncol(data_corrected)] <- NA

    for (i in 2:ncol(data)){
      LM_SNP_i <- lm(data[,i] ~ factor(data$pop))
      data_corrected[,i] <- LM_SNP_i$residuals
      colnames(data_corrected)[i]<-colnames(data)[i]
      #if(i%%50==0) print(i)
    }

    #Run 10 random forests with 1000 trees each and the optimized mtry value (0.1)
    r2vim <- Pomona::var.sel.r2vim(x = data_corrected[,3:ncol(data_corrected)],
                                   y=data_corrected$bio10,
                                   no.runs = 10,
                                   factor = 1,
                                   ntree = 1000,
                                   mtry.prop = 0.1,
                                   method = "ranger",
                                   type = "regression")
    # Build a data frame with results
    df<-data.frame(chrom = pos$chrom,
                   bp = pos$bp,
                   min.rel.vim = r2vim$info$rel.vim.min,
                   snp_name = pos$snp)
    DT::datatable(df,
                  filter = 'top',
                  extensions = c('Buttons', 'Scroller'),
                  options = list(scrollY = 650,
                                 scrollX = 500,
                                 deferRender = TRUE,
                                 scroller = TRUE,
                                 # paging = TRUE,
                                 # pageLength = 25,
                                 buttons = list('excel',
                                                list(extend = 'colvis',
                                                     targets = 0,
                                                     visible = FALSE)),
                                 dom = 'lBfrtip',
                                 fixedColumns = TRUE),
                  rownames = FALSE)
  })

  output$manhattan <- renderPlot({
    tryCatch(
      {
        Y1 <- read.table(input$filelfmmC$datapath, quote="\"", comment.char="", na.strings="9")
      },
      error = function(e){
        stop('Upload LFMM file')
      }
    )
    tryCatch(
      {
        bio10 <- read.table(input$filetxt2C$datapath, quote="\"", comment.char="")
      },
      error = function(e){
        stop('Upload TXT file')
      }
    )

    tryCatch(
      {
        PopData <- read.table(input$filetxt1C$datapath, header = T)
      },
      error = function(e){
        stop('Upload TXT file')
      }
    )

    tryCatch(
      {
        pos <- read.table(input$filetxt3C$datapath, header = T)
      },
      error = function(e){
        stop('Upload TXT file')
      }
    )

    pop <- PopData$pop

    #Create dataframe with genetic data, environmental data and SNP matrix
    data <- as.data.frame(cbind(pop,bio10,Y1))

    colnames(data)[1:2] <- c("pop","bio10")

    data_corrected <- data

    data_corrected[,2:ncol(data_corrected)] <- NA

    for (i in 2:ncol(data)){
      LM_SNP_i <- lm(data[,i] ~ factor(data$pop))
      data_corrected[,i] <- LM_SNP_i$residuals
      colnames(data_corrected)[i]<-colnames(data)[i]
      #if(i%%50==0) print(i)
    }

    r2vim <- Pomona::var.sel.r2vim(x = data_corrected[,3:ncol(data_corrected)],
                                   y=data_corrected$bio10,
                                   no.runs = 10,
                                   factor = 1,
                                   ntree = 1000,
                                   mtry.prop = 0.1,
                                   method = "ranger",
                                   type = "regression")
    # Build a data frame with results
    df<-data.frame(chrom = pos$chrom,
                   bp = pos$bp,
                   min.rel.vim = r2vim$info$rel.vim.min,
                   snp_name = pos$snp)

    #plot Manhattan plot with variable importance for each SNP
    plot(x=1:nrow(df),y=r2vim$info$rel.vim.min, xlab = 'SNP', ylab = 'Importance')
    abline(h=1, col="red")
  })

})

## User NBClust for gap method
## Arreglr logo agrosavia
