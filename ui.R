library(shiny)

shinyUI(fluidPage(
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "style_app.css")
  ),
  shinydashboard::dashboardPage(
    skin = 'purple',
    shinydashboard::dashboardHeader(
      title = 'GermVersity',
      titleWidth = 300,
      shinydashboard::dropdownMenu(type = 'messages',
                                   icon = icon('github'),
                                   shinydashboard::messageItem(
                                     from = 'Github',
                                     message = '',
                                     icon = icon('github'),
                                     href = 'https://github.com/GermVersity/GermVersity/'
                                   ))
    ),
    shinydashboard::dashboardSidebar(
      width = 300,
      shinydashboard::sidebarMenu(
        id = 'sidebar',
        style = "position: relative; overflow: visible;",
        tags$div(
          tags$img(height = "280px",
                   alt="Logo",
                   style = 'display: flex; justify-content: center',
                   src="Logo.png")
        ),
        shinydashboard::menuItem("Home",
                                 tabName = 'home',
                                 icon = icon('fas fa-house')),
        shinydashboard::menuItem('Module 1',
                                 tabName = 'module1',
                                 icon = icon('bezier-curve'),
                                 #shinydashboard::menuSubItem('Upload files',
                                 #           tabName = 'upload',
                                 #           icon = icon('upload')),
                                 shinydashboard::menuSubItem('Genotypic module',
                                                             tabName = 'genotypic1',
                                                             icon = icon('dna')),
                                 shinydashboard::menuSubItem('Phenotypic data',
                                                             tabName = 'phenotypic',
                                                             icon = icon('leaf')),
                                 shinydashboard::menuSubItem('Descriptor optimization',
                                                             tabName = 'optimization',
                                                             icon = icon('chart-bar'))
        ),
        shinydashboard::menuItem('Module 2',
                                 tabName = 'module2',
                                 icon = icon('bezier-curve'),
                                 shinydashboard::menuSubItem('Spatial DM',
                                                             tabName = 'spatialdm',
                                                             icon = icon('globe'))
                                 #,icon = icon('globe'))
        ),
        shinydashboard::menuItem('Module 3',
                                 tabName = 'module3',
                                 icon = icon('bezier-curve'),
                                 shinydashboard::menuSubItem('GAPIT',
                                                             tabName = 'gapit',
                                                             icon = icon('archway')),
                                 shinydashboard::menuSubItem('CPS',
                                                             tabName = 'cps',
                                                             icon = icon('bell')),
                                 shinydashboard::menuSubItem('LFMM',
                                                             tabName = 'lfmm',
                                                             icon = icon('brain'))
        ),
        shinydashboard::menuItem('Team',
                                 tabName = 'team',
                                 icon = icon('users')),
        HTML(paste0(
          "<br><br><br><br><br><br><br><br><br>",
          "<br><br><br><br><br><br><br><br><br>",
          "<br><br><br><br><br><br><br><br><br>",
          "<script>",
          "var today = new Date();",
          "var yyyy = today.getFullYear();",
          "</script>",
          "<p style = 'text-align: center;'><small>&copy; - <a href='#' target='_blank'>GermVersity </a> - <script>document.write(yyyy);</script></small></p>")
        )
      )

    ),
    shinydashboard::dashboardBody(
      shinydashboard::tabItems(
        shinydashboard::tabItem(
          tabName = 'home',
          fluidPage(
            sidebarLayout(
              sidebarPanel(
                tags$h3(class = "title",
                        "Application in development..."),
                tags$h3(),
                tags$p(class = "texto_inicial",
                       "This project transforms the way we address environmental
               challenges by more systematically and consistently studying
               the biodiversity of crop gene pools. Platforms for the
               discovery of agrobiodiversity, such as the one proposed here
               using AI, will allow the adaptability of crops not to depend
               exclusively on their intrinsic diversity, probably without
               sufficient previous adaptations to unusual climates.This project
               will harness modern AI technology to assist in the identification,
               characterization and utilization of new adaptive sources among
               thousands of landraces and related wild species, all stored in genebanks.."),
                tags$h4(class = "version",
                        "Version 1.8")
              ),
              mainPanel(
                tags$div(class = "images",
                         tags$img(src = "https://alimentro.agrosavia.co/Content/imagenes/logo-agrosavia.png", width = "230px", height = "50px"),
                         tags$img(src = "https://upload.wikimedia.org/wikipedia/commons/thumb/0/0a/Logotipo_de_la_Universidad_Nacional_de_Colombia.svg/1200px-Logotipo_de_la_Universidad_Nacional_de_Colombia.svg.png", width = "150px", height = "150px"),
                         tags$img(src = "https://upload.wikimedia.org/wikipedia/en/thumb/f/f4/Kansas_Jayhawks_logo.svg/1200px-Kansas_Jayhawks_logo.svg.png", width = "150px", height = "150px"),
                         tags$img(src = "https://climaloca.org/wp-content/uploads/alliance_logo_standard_cropped.png", width = "150px",  height = "70px")),
                tags$h1(class = "title",
                        "Welcome to GermVersity"),
                tags$h1(),
                tags$h4("Artificial intelligence as a tool for
                diversity analysis of plant genetic
                resources in genebanks"),
                tags$p(class = "texto_inicial",
                       "Plant genetic resources conserved in genebanks have
               the potential to harbor agronomically important traits
               that can be used to overcome current and future
               environmental challenges. However, despite their
               importance, these resources remain largely underutilized
               around the world. During the last years, the development
               of genome sequencing technologies has opened new avenues
               in the discovery of relevant genes for various agronomic
               traits in multiple crop species. Genebanks and the
               scientific community have generated massive amounts of
               phenotypic, genotypic, climatic, and geospatial data
               related to germplasm accessions. Currently, the challenge
               is to efficiently use these high-dimensional datasets to
               discover the true potential of plant accessions for breeding
               purposes. Although artificial intelligence (AI) tools may
               outperform competing methods in extracting hidden insights
               from high-dimensional data, their use in studying genebank
               accessions is still in its infancy. This situation could
               be due to the fact that these tools are relatively new
                and applications in the field of plant genetic resources
               are still incipient. To address this issue, we propose to
               use artificial intelligence tools to comprehensively
               analyze the phenotypic, genomic, and environmental data
               associated with plant germplasm accessions held in genebanks
               to gain a deep understanding of their genetic diversity,
               discover useful germplasm in climate change scenarios and
               promote their use.")
              )
            )
          )
        ),
        shinydashboard::tabItem(
          tabName = 'genotypic1',
          fluidPage(
            shiny::sidebarLayout(
              shiny::sidebarPanel(tags$h3("Germplasm bank management"),
                                  tags$br(),
                                  tags$p('In this demo, we are going to show how
                          to conduct basic genetic diversity
                          analyses of SNP data in a sample of
                          genebank accessions of Lima bean to
                          explore the genetic structure of the
                          sample. In this demo we are going to
                          learn how to calculate distance matrices
                          among individuals and populations,
                          visualize the distance matrices using
                          clustering algorithms (UPGMA and NJ),
                          how to conduct a principal component
                          analyses and a discriminant analysis
                          of principal components to assign
                          individuals into a k number of populations.'),
                                  tags$br(),
                                  shiny::fileInput("filevcf", "Choose a VCF file",
                                                   multiple = F,
                                                   accept = ".vcf",
                                                   buttonLabel = "Uploading..."),
                                  shiny::fileInput("filetxt", "Choose a TXT file",
                                                   multiple = F,
                                                   accept = ".txt",
                                                   buttonLabel = "Uploading...")
              ),
              mainPanel(
                tabsetPanel(
                  type = 'tabs',
                  shiny::tabPanel('Traditional Methods',
                                  tags$h2("Results traditional methods"),
                                  tags$br(),
                                  tags$h4('Let´s calculate a Euclidean distance matrix
                                  between individuals on the basis of the
                                  observed allele frequencies within individuals.
                                  This distance is not a genetic
                                  distance but a geometric distance since it
                                  does not assume any evolutionary model.'),
                                  tags$br(),
                                  tags$h4("Circular dendrogram"),
                                  shiny::plotOutput("plot_data", height = "500px",
                                                    dblclick = "double_click",
                                                    brush = brushOpts(
                                                      id = "brush_plot",
                                                      resetOnNew = TRUE)),
                                  tags$h4("Horizontal dendrogram"),
                                  shiny::plotOutput("plot_data2", height = "500px",
                                                    dblclick = "double_click2",
                                                    brush = brushOpts(
                                                      id = "brush_plot",
                                                      resetOnNew = TRUE))
                  ),
                  shiny::tabPanel('Distance Matrix',
                                  tags$br(),
                                  tags$p("Genetic distances measure how genetically
                                 similar are individuals or populations and
                                 are based on a specific evolutionary model.
                                 Genetic distances can be seen as summary
                                 statistics because they take the whole data
                                 set (for example, data from SNP loci) and
                                 summarize the genetic differentiation between
                                 samples (individuals or populations) in a
                                 value. One of the most used genetic distances
                                 is Nei's genetic distance which measures
                                 the genetic distance among populations on
                                 the basis of their genetic identity, namely
                                 the proportion of alleles that are shared
                                 between populations. For our data set, we
                                 have defined gene pools as populations in the
                                 genind object called LimaBeanData2.
                                 To calculate Nei's genetic distances among
                                 populations we need to add the strata to
                                 the genind object (in this case on the basis
                                 of the gene pool) to define the populations.
                                 Finally, we will build with this distance
                                 matrix a UPGMA tree and a neighbor-joining (
                                 NJ) topology, applying 1000 bootstrap permutations
                                 to get statistical support of the groups in the topologies."),
                                  tags$br(),
                                  tags$h4("UPGMA topology"),
                                  shiny::plotOutput("distance_matrix1", height = "500px"),
                                  tags$h4("NJ topology"),
                                  shiny::plotOutput("distance_matrix2", height = "500px")
                  ),
                  shiny::tabPanel('PCA',
                                  tags$br(),
                                  tags$p('PCA (Principal component analysis) is a multivariate
                                  technique that summarizes
                                 the information provided by the genetic markers
                                 (for example, SNPs) into a few set of components.
                                 We can apply PCA to genlight objects with the
                                 function glPca of the package adgenet'),
                                  tags$br(),
                                  tags$h4("PCA plot"),
                                  shiny::plotOutput("plotPCA", height = "500px")
                  ),
                  shiny::tabPanel('DPCA',
                                  tags$br(),
                                  tags$p('To identify genetic clusters, we can apply
                                 another multivariate approach known as DAPC.
                                 This approach is convenient when we are more
                                 interested in describing the diversity
                                 among groups of individuals than within groups.
                                 This approach is focused on finding discriminant
                                 functions that better describe the differences
                                 among groups, while minimizing the differences
                                 within groups. However, to find the discriminant
                                 functions, DAPC needs the groups to be known a
                                 priori and in many cases we just do not know
                                 how many groups are present in our sample.
                                 To address this issue, the adegenet package
                                 implements the function find.clusters (to find
                                 the number of clusters by running first a
                                 PCA and then using the k-means algorithm)
                                 and the function dapc to establish how are
                                 the relationships among the clusters.'),
                                  tags$br(),
                                  tags$h4("DPCA plot"),
                                  sliderInput("npca", "Number of PCA", value = 5, min = 1, max = 30),
                                  sliderInput("nda", "Number of discriminant functions", value = 5, min = 1, max = 30),
                                  shiny::plotOutput("plotDPCA", height = "500px"),
                                  tags$h4("Clusters"),
                                  sliderInput("npca1", "Number of PCA", value = 5, min = 1, max = 30),
                                  sliderInput("nclust", "Number of cluster", value = 5, min = 1, max = 30),
                                  shiny::dataTableOutput("table1D"),
                                  tags$h4("DPCA custom plot"),
                                  shiny::plotOutput("plotDPCA2", height = "500px"),
                                  tags$h4("Membership probabilities"),
                                  shiny::plotOutput("plotDPCA3", height = "500px")
                  ),
                  shiny::tabPanel('BDSpl',
                                  tags$h3("Basic diverity statistics per locus"),
                                  tags$br(),
                                  tags$p('In this section we are going to estimate
                                 basic diversity statistics per locus in the
                                 whole sample. We are going to use the program
                                 adegenet to estimate the number of alleles
                                 per locus (NA), observed heterosygosity per
                                 locus (Hobs) and expected (Hexp) heterozygosity
                                 per locus. For SNP makers, NA might not be
                                 very useful because these loci are expected
                                 to be biallelic in populations (according to
                                 the infinite site mutational model) and we
                                 have also filtered the vcf file to include
                                 only biallelic SNPs. Hobs is the proportion
                                 of heterozygous individuals that were observed
                                 in the locus. Hexp is the heterozygosity we
                                 expect to observe in the locus assuming that
                                 the population (in this case the whole sample)
                                 is in Hardy-Weinberg equilibrium (HWE)in that
                                 locus. A significant difference among Hobs and
                                 Hexp means that the locus is not in HWE'),
                                  tags$br(),
                                  tags$h4("Hobs per locus"),
                                  shiny::plotOutput("plot_d1", height = "500px"),
                                  tags$h4("Hexp per locus"),
                                  shiny::plotOutput("plot_d2", height = "500px"),
                                  tags$h4("Hexp as a function of Hobs per locus"),
                                  shiny::plotOutput("plot_d3", height = "500px"),
                                  tags$h4("Testing the difference among Hexp and Hobs per locus")
                  ),
                  shiny::tabPanel('BDSpc',
                                  tags$h3("Basic diversity statistics per cluster (or population)"),
                                  tags$br(),
                                  tags$p('Once we have defined the number of clusters
                                 (or populations) in a sample, it is a good idea
                                 to estimate genetic diversity indexes for each
                                 cluster, for example if we want to find out
                                 which cluster is the most or least diverse.
                                 For doing this, we are going to use the program
                                 hierfstat to estimate observed heterozygosity
                                 (Ho) and mean gene diversities (Hs) within
                                 populations. A significant difference among
                                 Ho and Hs means that the population is not in
                                 HWE. Because different populations may not be
                                 in HWE for different reasons, diversity among
                                 populations may be hardly compared on the basis
                                 of Ho. This is because Hs is preferred to compare
                                 the genetic diversity between populations because
                                 with Hs all populations are assumed to be in HWE.
                                 It is also useful to estimate the fixation index
                                 Fis to measure the deviation to the assumption of
                                 HWE within each population.'),
                                  tags$br(),
                                  tags$h4("Summary Hs"),
                                  tags$p('To summarize Hs values within populations. We can observe
                                 that populations are different in their genetic diversity,
                                 with the highest values for the domesticated admixed,
                                 followed by the Andean wild gene pool AI'),
                                  tags$h4("Summary Ho"),
                                  shiny::verbatimTextOutput("summary1"),
                                  tags$p('To summarize Ho values within populations. We can
                                 observe that within populations, Ho is very los, as
                                 expected for a selfing species as Lima bean.'),
                                  tags$h4("Summary Fis"),
                                  shiny::verbatimTextOutput("summary2"),
                                  tags$p('To summarize Fis values within populations.
                                 Fis is calculated as (Hs-Ho)/Hs. We can see
                                 that Fis values are positive for all populations
                                 and in some of them close to 1, as expected for
                                 a predominantly autogamous species as Lima bean.'),
                                  shiny::verbatimTextOutput("summary3"),
                                  tags$p('A useful plot to compare genetic diversity among genepools'),
                                  shiny::plotOutput('boxplot')
                  ),
                  shiny::tabPanel('Results ML',
                                  tags$h3("Genetic divergence among populations"),
                                  tags$br(),
                                  tags$p('We may also be interested in finding out how genetically
                                 divergent are the populations. For this, we can calculate
                                 the fixation index Fst, which measures the difference
                                 between the mean gene diversity within populations (Hs)
                                 and the total genetic diversity expected in the whole
                                 population (Ht). The higher the difference between Ht
                                 and Hs, the higher the fixation index Fst and therefore
                                 the higher the populations we are comparing. We will
                                 use the program hierfstat to calculate Fst indexes.'),
                                  tags$br(),
                                  tags$p('A useful plot to explore the overall difference between
                                 Hs and Ht. We can see that Hs is lower than Ht, therefore
                                 we can envision that there will be certain degree of
                                 differentiation among populations (in this case we are
                                 comparing genepools in Lima bean) as measured by Fst
                                 (Fst=(Ht-Hs)/Ht).'),
                                  shiny::plotOutput('boxplot1'),
                                  tags$p('to calculate Fis and Fst on one level hierarchy
                                 (the populations -or genepools- grouped within the
                                 total population) by the method of Weir and Cockerham (1984)'),
                                  shiny::verbatimTextOutput("summary4"),
                                  tags$p('To obtain pariwise Fst values among genepools'),
                                  shiny::verbatimTextOutput("summary5"),
                                  tags$p('Cluster Validity by NbCLust (and factoextra)
                                using 30 indices from the scientific literature'),
                                  shiny::plotOutput('plot1'),
                                  tags$p('Cluster Validity by OptCluster (an improvement of ClValid)'),
                                  shiny::plotOutput('plot2'),
                                  tags$p('Vizualization of OptCluster Output by means of UPGMA dendogram'),
                                  shiny::plotOutput('plot3')
                  )
                )
              )
            )
          )
        ),
        shinydashboard::tabItem(
          tabName = 'phenotypic',
          fluidPage(
            sidebarLayout(
              sidebarPanel(tags$h3("Phenotypic analysis in genebanks"),
                           tags$p('This work is in the framework of the project "
                          To use AI tools for multidimensional data analysis
                          to understand and promote biodiversity in plant
                          germplasm banks (Phaseolus case)" led by AGROSAVIA,
                          UNAL and Alianza Bioversity-CIAT from a cooperation
                          agreement. the project objective is to "To use AI
                          tools for multidimensional data analysis to understand
                          and promote biodiversity in plant germplasm banks
                          (Phaseolus case).'),
                           shiny::fileInput("filecsv", "Choose a CVS file",
                                            multiple = F,
                                            accept = c("text/csv",
                                                       "text/comma-separated-values, text/plain",
                                                       ".csv"),
                                            buttonLabel = "Uploading...")
              ),
              mainPanel(
                tags$h3("Model"),
                verbatimTextOutput('modelF'),
                plotOutput("modFit"),
                tags$h3('Confusion matrix'),
                tags$p('Table of confusion of the training data of
                       the common bean races classification.'),
                plotOutput('con_matrix')
              )
            )
          )
        ),
        shinydashboard::tabItem(
          tabName = 'optimization',
          fluidPage(
            sidebarLayout(
              sidebarPanel(tags$h3("Oescriptor optimization"),
                           tags$p('Distribution of minimum depth and number of trees for
                          each of the descriptors evaluated. It is observed in
                          order of importance of the descriptor, the number of
                          nodes associated with each descriptor and the interaction
                          between them.'),
                           shiny::fileInput("filecsv1", "Choose a CVS file",
                                            multiple = F,
                                            accept = c("text/csv",
                                                       "text/comma-separated-values, text/plain",
                                                       ".csv"),
                                            buttonLabel = "Uploading...")
              ),
              mainPanel(
                tags$h3("ML Results"),
                plotOutput('min_depth'),
                tags$hr(),
                plotOutput('multi'),
                tags$hr(),
                plotOutput('inter')
              )
            )
          )
        ),
        shinydashboard::tabItem(
          tabName = 'team',
          fluidPage(
            tags$h1(class = "title", 'Miembros del equipo'),
            tags$br(),
            tags$div(class = "cardsteam",
                     tags$div(class = "card",
                              tags$img(src = "Joaquin.jpg", width = "100px"),
                              tags$div(class = "containerteam",
                                       tags$h4(tags$b('Joaquin Guillermo Ramirez Gil')),
                                       tags$img(src = 'unalblack.png', width = '100px'),
                                       tags$h6("ORCID: https://orcid.org/0000-0002-0162-3598")
                                       )
                              ),
                     tags$div(class = "card",
                              tags$img(src = "Camilo.jpg", width = "100px"),
                              tags$div(class = "containerteam",
                                       tags$h4(tags$b('Juan Camilo Henao-Rojas')),
                                       tags$img(src = 'agrosavia.png', width = '100px'),
                                       tags$h6("ORCID: https://orcid.org/0000-0003-0007-6809")
                              )
                     ),
                     tags$div(class = "card",
                              tags$img(src = "Andres.jpg", width = "100px"),
                              tags$div(class = "containerteam",
                                       tags$h4(tags$b('Andrés Cortés-Vera')),
                                       tags$img(src = 'agrosavia.png', width = '100px'),
                                       tags$h6("ORCID: https://orcid.org/0000-0003-4178-0675")
                              )
                     )
                     
                     ),
            tags$br(),
            tags$br(),
            tags$div(class = "cardsteam",
                     tags$div(class = "card",
                              tags$img(src = "Isabel.jpeg", width = "100px"),
                              tags$div(class = "containerteam",
                                       tags$h4(tags$b('Maria Isabel Chacon-Sanchez')),
                                       tags$img(src = 'unalblack.png', width = '100px'),
                                       tags$h6("ORCID: https://orcid.org/0000-0001-7781-6129")
                              )
                     ),
                     tags$div(class = "card",
                              tags$img(src = "Diego.jpg", width = "100px"),
                              tags$div(class = "containerteam",
                                       tags$h4(tags$b('Diego Felipe Conejo')),
                                       tags$img(src = 'ciat.png', width = '100px'),
                                       tags$h6("ORCID: https://orcid.org/0000-0001-7129-4016")
                              )
                     ),
                     tags$div(class = "card",
                              tags$img(src = "Marlon.jpeg", width = "100px"),
                              tags$div(class = "containerteam",
                                       tags$h4(tags$b('Marlon E. Cobos')),
                                       tags$img(src = 'kansas.png', width = '100px'),
                                       tags$h6("ORCID: https://orcid.org/0000-0002-2611-1767")
                              )
                     )
                     ),
            tags$br(),
            tags$br(),
            tags$div(class = "cardsteam",
                     tags$div(class = "card",
                              tags$img(src = "Luis.jpeg", width = "100px"),
                              tags$div(class = "containerteam",
                                       tags$h4(tags$b('Luis Felipe Lopez')),
                                       tags$img(src = 'agrosavia.png', width = '100px'),
                                       tags$h6("ORCID: https://orcid.org/0000-0002-4967-6955")
                              )
                     ),
                     tags$div(class = "card",
                              tags$img(src = "Paula.jpg", width = "100px"),
                              tags$div(class = "containerteam",
                                       tags$h4(tags$b('Paula Helena Reyes Herrera')),
                                       tags$img(src = 'agrosavia.png', width = '100px'),
                                       tags$h6("ORCID:")
                              )
                     ),
                     tags$div(class = "card",
                              tags$img(src = "Kevin.jpg", width = "100px"),
                              tags$div(class = "containerteam",
                                       tags$h4(tags$b('Kevin Steven Quiroga Benavides')),
                                       tags$img(src = 'unalblack.png', width = '100px'),
                                       tags$h6("ORCID: https://orcid.org/0000-0002-2748-942X")
                              )
                     )
                     ),
            tags$br(),
            tags$div(class = "footer",
                     tags$div(class = 'footer_1', 
                              tags$h5(class="title", "Versión 3.0.1"),
                              tags$p(class = "footer_2", "© Copyright GermVersity. All rights reserved")))
            )
          ),
        shinydashboard::tabItem(
          tabName = 'gapit',
          fluidPage(
            shiny::sidebarLayout(
              shiny::sidebarPanel(tags$h3("GAPIT"),
                                  tags$br(),
                                  tags$p('GAPIT is a Genome Association and Prediction
                                 Integrated Tool freely available for Public
                                 since 2011 (Lipka et al. 2012). It has been
                                 updated frequently to incorporate the state
                                 of art methods for Genome Wide Association
                                 Study (GWAS). Currently, Dr. Jiabo Wang is
                                 leading the new development, GAPIT version
                                 3. In addition to the MLM methods, the new
                                 version implemented two new GWAS methods as
                                 FarmCPU (Liu et al. 2016) and BLINK (Huang
                                 et al. 2019). López-Hernández & Cortés (2019)
                                 explored this tool to obtain the Genome-Environment
                                 Associations (GEA) in crops (Common bean)
                                 using Genotyping by Sequencing (GBS) and
                                 Environmental Indices of Heat Stress.'),
                                  tags$br(),
                                  shiny::fileInput("filexlsxG", "Choose a Excel (xlsx) file",
                                                   multiple = F,
                                                   accept = ".xlsx",
                                                   buttonLabel = "Uploading..."),
                                  shiny::fileInput("filetxtG", "Choose a TXT file",
                                                   multiple = F,
                                                   accept = ".txt",
                                                   buttonLabel = "Uploading..."),
                                  tags$br(),
                                  shiny::actionButton('gob', 'Click to continue..')
              ),
              mainPanel(tags$h2("Results Downloading..."),
                        tags$br(),
                        tags$h3('See the files with the result in the next folder: ~/GAPIT/'),
                        tags$br(),
                        shiny::verbatimTextOutput('ready'),
                        tags$br(),
                        shiny::actionButton('Manhattan', 'See plot..'),
                        tags$br(),
                        tags$h3('Manhattan Circular plot for MLM'),
                        tags$br(),
                        shiny::plotOutput('plot_data1G', height = '800px'),
                        tags$h3('QQ plot for MLM'),
                        tags$br(),
                        shiny::plotOutput('plot_data2G', height = '800px'),
                        tags$br(),
                        tags$h3('Manhattan Circular plot for Farm-CPU'),
                        tags$br(),
                        shiny::plotOutput('plot_data3G', height = '800px'),
                        tags$h3('QQ plot for Farm-CPU'),
                        tags$br(),
                        shiny::plotOutput('plot_data4G', height = '800px'),
                        tags$br(),
                        tags$h3('Manhattan Circular plot for Blink'),
                        tags$br(),
                        shiny::plotOutput('plot_data5G', height = '800px'),
                        tags$h3('QQ plot for Blink'),
                        tags$br(),
                        shiny::plotOutput('plot_data6G', height = '800px')

              )
            )
          )
        ),
        shinydashboard::tabItem(
          tabName = 'spatialdm',
          fluidPage(
            shiny::sidebarLayout(
              shiny::sidebarPanel(tags$h3("Spatial Distribution Modeling (SDM)"),
                                  tags$br(),
                                  tags$p("Exploration of the spatial distribution modeling
                                 using the Maxent Algorithm (Phillips et al. 2006,
                                 Phillips et al. 2008, Phillips et al. 2017)
                                 by means of the R-package ENMeval (Kass et al.
                                 2021, https://doi.org/10.1111/2041-210X.13628)"),
                                  tags$br(),
                                  shiny::fileInput("filexlsx", "Choose a Excel file",
                                                   multiple = F,
                                                   accept = ".xlsx",
                                                   buttonLabel = "Uploading..."),
                                  shiny::fileInput("filebio", "Choose the Layers files",
                                                   multiple = T,
                                                   accept = c(".bil",".hdr"),
                                                   buttonLabel = "Uploading..."),
                                  tags$hr(),
                                  tags$p('Use WorldClim data (2022)'),
                                  shiny::actionButton('WC', 'Click me!')
              ),
              mainPanel(tags$h2("Results"),
                        tags$br(),
                        tags$h4("Responses of the variables (predictors)"),
                        shiny::plotOutput("plot_dataS", height = "500px"),
                        tags$h4("Plotting the prediction modeling using the ENMevaluate results"),
                        shiny::plotOutput("plot_data2S", height = "500px"),
                        tags$h4("Plotting the prediction modeling using the GLM results"),
                        plotOutput("plot_data3S", height = "500px")

              )
            )
          )
        ),
        shinydashboard::tabItem(
          tabName = 'cps',
          fluidPage(
            shiny::sidebarLayout(
              shiny::sidebarPanel(tags$h3("Correcting for population stratification"),
                                  tags$br(),
                                  tags$p('We are going to use the approach of Zhao et al.
                                 (2012) which consists in first, regressing
                                 (at each SNP locus and by applying a generalized
                                 linear model, GLM), the environmental variable
                                 and the genotype against population information.
                                 The residuals of the models are then used to
                                 adjust the values of the environmental variable
                                 and the genotype before performing the Random
                                 Forest analysis.'),
                                  tags$br(),

                                  shiny::fileInput("filelfmmC", "Choose a LFMM file",
                                                   multiple = F,
                                                   accept = ".lfmm",
                                                   buttonLabel = "Uploading..."),
                                  tags$br(),

                                  shiny::fileInput("filetxt1C", "Choose a TXT file (genepool)",
                                                   multiple = F,
                                                   accept = ".txt",
                                                   buttonLabel = "Uploading..."),
                                  tags$br(),

                                  shiny::fileInput("filetxt2C", "Choose a TXT file (environment)",
                                                   multiple = F,
                                                   accept = ".txt",
                                                   buttonLabel = "Uploading..."),
                                  tags$br(),

                                  shiny::fileInput("filetxt3C", "Choose a TXT file (position)",
                                                   multiple = F,
                                                   accept = ".txt",
                                                   buttonLabel = "Uploading..."),

                                  tags$br()
              ),
              mainPanel(
                shiny::tabsetPanel(
                  type = 'tabs',
                  shiny::tabPanel('Distribution',
                                  tags$br(),
                                  tags$h3('Explore the overall distribution of the phenotype'),
                                  plotOutput('plot_dataC'),
                                  tags$br(),
                                  tags$h3('Should all be zero (0)'),
                                  verbatimTextOutput('corrected')
                  ),
                  shiny::tabPanel('Random Forest',
                                  tags$br(),
                                  tags$p('RF algorithms should be optimized in order to
                                 maximize the percentage of variance explained
                                 (PVE, for continuos variables) and to ensure results
                                 are repeatable. We are going to adjuts two parameters:
                                 ntree (number of trees grown per forest) and mtry
                                 (number of predictors to randomly sample at each node).
                                 We should establish the values of ntree and mtry where
                                 PVE reaches a plateau. We will try values of ntree from
                                 100 to 10000 (in increments of 100), and values of mtry
                                 of sqrt(p), 2*sqrt(p), 0.1(p), 0.2(p), p/3, and p (
                                 p=number of loci). '),
                                  tags$br(),
                                  plotOutput('plot_data2C')
                  ),
                  shiny::tabPanel('Correlation',
                                  tags$br(),
                                  plotOutput('hist1')),
                  shiny::tabPanel('Results SNP',
                                  tags$br(),
                                  DT::dataTableOutput('data_table'),
                                  tags$br(),
                                  plotOutput('manhattan'))
                )

              )
            )
          )
        ),
        shinydashboard::tabItem(
          tabName = 'lfmm',
          fluidPage(
            shiny::sidebarLayout(
              shiny::sidebarPanel(tags$h3("Latent factor mixed model analysis"),
                                  tags$br(),
                                  tags$p('Code to run a latent factor mixed model analysis (LFMM)
                                 using data from 259 accessions of wild Lima bean
                                 (Phaseolus lunatus L.), that belong to four wild gene
                                 pools (two Mesoamerican: MI and MII, and two Andean:AI
                                 and AII), that have been genotyped at 10668 SNP loci.
                                 LFMM applies a regression model to carry association
                                 tests among genetic variants and environmental variables.
                                 Correction for confounding effects, such as population
                                 structure, is done by including (unobserved) latent
                                 factors (set with K) which are estimated in parallel
                                 with the response variables (SNP loci) and the environmental
                                 variable of interest.'),
                                  tags$br(),
                                  shiny::fileInput("filevcfL", "Choose a VCF file",
                                                   multiple = F,
                                                   accept = ".vcf",
                                                   buttonLabel = "Uploading..."),
                                  shiny::fileInput("filelfmmL", "Choose a LFMM file",
                                                   multiple = F,
                                                   accept = ".lfmm",
                                                   buttonLabel = "Uploading..."),
                                  shiny::fileInput("filetxtL", "Choose a TXT file (environment)",
                                                   multiple = F,
                                                   accept = ".txt",
                                                   buttonLabel = "Uploading..."),
                                  shiny::fileInput("filetxt1L", "Choose a TXT file (position)",
                                                   multiple = F,
                                                   accept = ".txt",
                                                   buttonLabel = "Uploading..."),
                                  tags$br()
              ),
              mainPanel(
                shiny::tabsetPanel(
                  type = 'tabs',
                  shiny::tabPanel('Descriptive',
                                  tags$h2("Results"),
                                  tags$br(),
                                  tags$h3('Estimation of K (number of populations)'),
                                  tags$br(),
                                  plotOutput('plot_dataL'),
                                  tags$br(),
                                  tags$p('The screeplot indicates that there are around K=6-8 main components in the data.'),
                                  tags$br(),
                                  plotOutput('plot_data1L')),
                  shiny::tabPanel('False Discovery Rate',
                                  tags$br(),
                                  tags$h3("Identify LFMM candidates using False Discovery Rate"),
                                  tags$br(),
                                  tags$p('The ridge_lfmm function returns an object that contains
                                  the following matrices:
                                  Matrix U: matrix of latent variables
                                  Matrix B: matrix of effects of all explanatory variables (loci)
                                  Matrix V: matrix of loadings for all latent variables'),
                                  tags$br(),
                                  tags$h3('Unadjusted p-values.'),
                                  tags$br(),
                                  plotOutput('plot_data2L'),
                                  tags$br(),
                                  tags$h3('GIF-adjusted p-values.'),
                                  tags$br(),
                                  plotOutput('plot_data3L'),
                                  tags$br(),
                                  tags$h3('Unadjusted p-values QQplots.'),
                                  tags$br(),
                                  plotOutput('plot_data4L'),
                                  tags$br(),
                                  tags$h3('GIF-adjusted p-values QQplots.'),
                                  tags$br(),
                                  plotOutput('plot_data5L')),
                  shiny::tabPanel('Combining results',
                                  tags$br(),
                                  tags$h3("Combining results with position matrices"),
                                  tags$br(),
                                  tags$h3('Unadjusted and Calibrated p-values.'),
                                  tags$br(),
                                  DT::dataTableOutput('table1L'),
                                  tags$br(),
                                  tags$h3('Manhattan Plot'),
                                  tags$br(),
                                  plotOutput('plot_data6L'),
                                  tags$br(),
                                  tags$h3('Effect size of the environmental variable of interest on each SNP'),
                                  tags$br(),
                                  plotOutput('plot_data7L'))
                )
              )
            )
          )
        )
      )
    )
  )

))
