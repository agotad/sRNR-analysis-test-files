library(shiny)
library(shinyWidgets)
library(DT)
library(data.table)
library(parallel)
library(shinythemes)
library(APERO)
library(Rsamtools)
library(reshape2)
library(dplyr)
library(parallel)
library(ggplot2)
library(VennDiagram)
library(UpSetR)
library(tidyr)
library(corrplot)
library(ggseqlogo)
library(msa)
library(writexl)
library(bamsignals)
library(GenomicRanges)
library(reticulate)
library(rtracklayer)
library(readxl)
library(readr)
library(openxlsx)
library(purrr)
library(stringr)
library(Gviz)
library(bslib)
library(RColorBrewer)
library(shinyBS)
library(shinyhelper)
library(shinydashboard)
library(bs4Dash)

tags <- htmltools::tags


options(shiny.maxRequestSize = 900 * 1024^2)
mylist <- list.files("./www/")



my_theme <- bs_theme(
  bootswatch = "minty",
  base_font = font_google("Roboto"),
  font_scale = 1.5
)


tags_spinner <- tags$head(
  tags$style(HTML("  
#loading-spinner {
  display: none;
  position: fixed;
  z-index: 9999;
  top: 0;
  left: 0;
  width: 100%;
  height: 100%;
  background-color: rgba(255, 255, 255, 0.8);
}
.spinner-border {
  position: absolute;
  top: 50%;
  left: 50%;
  transform: translate(-50%, -50%);
  width: 4rem;
  height: 4rem;
  border: 0.4em solid rgba(0, 0, 0, 0.1);
  border-top: 0.4em solid #3498db;
  border-radius: 50%;
  animation: spin 1s linear infinite;
}
@keyframes spin {
  to { transform: translate(-50%, -50%) rotate(360deg); }
}")),

  tags$script(HTML("  
Shiny.addCustomMessageHandler('show_spinner', function(message) {
  document.getElementById('loading-spinner').style.display = 'block';
});
Shiny.addCustomMessageHandler('hide_spinner', function(message) {
  document.getElementById('loading-spinner').style.display = 'none';
});
"))
)


ui <- bs4DashPage(
  freshTheme = my_theme,
  title = "APERO sRNA Analysis",

  header = bs4DashNavbar(
    skin = "light",
    status = "white",
    border = TRUE,
    controlbarIcon = icon("cogs")
  ),

  sidebar = bs4DashSidebar(
    skin = "light",
    status = "primary",
    title = "APERO",
    brandColor = "primary",
    bs4SidebarMenu(
      bs4SidebarMenuItem("Input", tabName = "main", icon = icon("dna")),
      bs4SidebarMenuItem("Tables", tabName = "tables", icon = icon("table")),
      bs4SidebarMenuItem("Plots", tabName = "plots", icon = icon("chart-bar")),
      bs4SidebarMenuItem("Genome browser", tabName = "genome", icon = icon("map")),
      bs4SidebarMenuItem("Secondary structures", tabName = "structures", icon = icon("project-diagram"))
    )
  ),

  body = bs4DashBody(
    tags_spinner,
    div(id = "loading-spinner", class = "spinner-overlay",
    div(class = "spinner-border"),
    h4("Processing, please wait...", style = "position: absolute; top: 60%; left: 50%; transform: translate(-50%, -50%); color: #555;")),

    bs4TabItems(
      bs4TabItem(tabName = "main",
        fluidRow(
          bs4Card(
            title = "Input Options",
            width = 12,
            prettyRadioButtons(
              inputId = "input_option",
              label = "Choose input option:",
              choices = c("APERO analysis (First step)", "APERO results analysis (Second step)"),
              selected = "APERO analysis (First step)",
              shape = "curve",
              status = "primary",
              animation = "pulse"
            )
          )
        ),
       

        conditionalPanel(
          condition = "input.input_option == 'APERO results analysis (Second step)'",
          fluidRow(
            bs4Card(
              title = "Upload Files",
              width = 12,
              fileInput("group1_files", "Upload APERO Group 1 CSVs", multiple = TRUE, accept = ".csv"),
              textInput("group1_label", "Label for Group 1"),
              fileInput("group2_files", "Upload APERO Group 2 CSVs", multiple = TRUE, accept = ".csv"),
              textInput("group2_label", "Label for Group 2"),
              fileInput("fasta_file", "Upload FASTA Reference Genome (.fasta)", accept = c(".fasta", ".fa")),
              fileInput("file_ptt", "Upload PTT File (.ptt)", accept = ".ptt"),
              fileInput("fasta_files", "FASTA files for BLAST", multiple = TRUE, accept = c(".fasta", ".fa")),
              prettyRadioButtons(
  inputId = "species_choice",
  label = "Choose the species for analysis:",
  choices = c("Lactococcus lactis", "Lactococcus casei"),
  selected = "Lactococcus lactis",
  shape = "curve",
  status = "success",
  animation = "pulse"
))
          ),

          fluidRow(
            bs4Card(
              title = "Analysis Parameters",
              width = 12,
              numericInput("min_length", "Minimum Sequence Length:", value = 50),
              numericInput("max_length", "Maximum Sequence Length:", value = 500),
              numericInput("max_promoter_distance", "Max Distance to Promoter:", value = 20),
              numericInput("max_terminator_distance", "Max Distance to Terminator:", value = 20),
              actionButton("process_btn", "Process Files", class = "btn-primary")
            )
          )
        ),

        conditionalPanel(
          condition = "input.input_option == 'APERO analysis (First step)'",
          fluidRow(
            bs4Card(
              title = "BAM and Genome Inputs",
              width = 12,
              fileInput("bamFile", "Upload BAM file", accept = ".bam"),
              fileInput("annotationFile", "Upload Genome (PTT format)", accept = ".ptt"),
              fileInput("gff3_path", "Upload GFF3 File", accept = ".gff3")
            )
          ),
          fluidRow(
            bs4Card(
              title = "Detection Parameters",
              width = 12,
              numericInput("wmax", "The maximal accepted width of a start peak (wmax)", 10),
              numericInput("min_dist", "Minimum distance (nt) between two detected peaks to consider them as separate (min_dist)", 10),
              numericInput("enrichment", "Threshold for enrichment over background to detect a peak. Lower values make detection more sensitive (enrichment)", 0.3),
              numericInput("min_read_number", "Minimum number of reads required at a position to consider it a potential start site (min_read_number)", 20),
              numericInput("genome_size", "The total length of the reference genome, used to normalize read density or define limits (genome_size)", 2529478),
              numericInput("readthrough_proportion", "Allowed proportion of reads that extend beyond the detected 3′ end (readthrough_proportion). Lower = stricter.", 0.01),
              numericInput("Fmin", "Minimum F-score threshold (used internally to filter poor predictions) (Fmin). Can be NA. ", NA),
              numericInput("thread_number", "Number of CPU threads to use in parallel computation (thread_number)", 8),
              sliderInput("tRNA_thresh", "tRNA overlap threshold (%)", 0, 100, 50),
              sliderInput("rRNA_thresh", "rRNA overlap threshold (%)", 0, 100, 50),
              sliderInput("gene_thresh", "Gene overlap threshold (%)", 0, 100, 70),
              sliderInput("mRNA_thresh", "mRNA overlap threshold (%)", 0, 100, 50),
              actionButton("runAnalysis", "Run APERO", class = "btn-primary"),
              downloadButton("downloadResults", "Download Results")
            )
          )
        )
      ),

      bs4TabItem(tabName = "tables",
      fluidRow(
          bs4ValueBoxOutput("sRNA_count", width = 4),
          bs4ValueBoxOutput("promoter_count", width = 4),
          bs4ValueBoxOutput("terminator_count", width = 4)
        ),
        fluidRow(
          bs4Card(title = "Info table", width = 12, DTOutput("filesTable"))
        ),
        fluidRow(
          bs4Card(width = 12, downloadButton("download_all", "Download All Results in Excel format"))
        ),
        fluidRow(
          bs4Card(title = "Main Table", width = 12, div(style = "overflow-x: auto; width: 100%;", DTOutput("lentele_table")))
        ),
        fluidRow(
          bs4Card(title = "Promoters", width = 12, div(style = "overflow-x: auto; width: 100%;", DTOutput("promoter_table")))
        ),
        fluidRow(
          bs4Card(title = "Terminators", width = 12, div(style = "overflow-x: auto; width: 100%;",DTOutput("terminators_table")))
        ),
        fluidRow(
          bs4Card(title = "Genomic Context", width = 12, div(style = "overflow-x: auto; width: 100%;",DTOutput("genomic_context")))
        ),
        fluidRow(
          bs4Card(title = "Structure info", width = 12, div(style = "overflow-x: auto; width: 100%;",DTOutput("structure_table")))
        ),
        fluidRow(
          bs4Card(title = "BLAST Results", width = 12, div(style = "overflow-x: auto; width: 100%;",DTOutput("blast_results")))
        )
      ),

      bs4TabItem(tabName = "plots",
        fluidRow(
          bs4Card(title = "sRNA Length Histogram", width = 12, plotOutput("length_histogram"))),
fluidRow(
          bs4Card(title = "GC Content Histogram", width = 12, plotOutput("gc_content_plot"))
        ),
        fluidRow(
          bs4Card(title = "Length vs GC Content", width = 12, plotOutput("lengthvsgc"))),
fluidRow(
          bs4Card(title = "Strand vs GC Content", width = 12, plotOutput("strandvsgc"))
        ),
        fluidRow(
          bs4Card(title = "Venn Plot", width = 12, plotOutput("venn_plot"))),
fluidRow(
          bs4Card(title = "UpSet Plot", width = 12, plotOutput("upset_plot"))
        ),
        fluidRow(
          bs4Card(title = "Length by Group", width = 12, plotOutput("length_by_group"))),
fluidRow(
          bs4Card(title = "sRNA Group Distribution", width = 12, plotOutput("groupcount"))
        ),
        fluidRow(
          bs4Card(title="Free 3'end number distribution", width = 12, plotOutput("free3_distribution"))),
fluidRow(
          bs4Card(title="Different groups", width = 12, plotOutput("upset_group_plot"))),

        fluidRow(
          bs4Card(title="Length in different classes", width = 12, plotOutput("lengthclass"))),
fluidRow(
          bs4Card(title="Different groups", width = 12, selectInput("promoter_region", "Select promoter region:", choices = c("Minus10", "Minus35"), selected = "Minus10"),
                 plotOutput("promoter_logo_plot"))),
      ),

      bs4TabItem(tabName = "genome",
        fluidRow(
          bs4Card(title = "Genome Coordinate Selection", width = 6,
            numericInput("region_start", "Start coordinate", 10000),
            numericInput("region_end", "End coordinate", 100000),
            helpText("Blue genes = '-', red genes = '+'"),
            downloadButton("download_genome_browser", "Download Genome PDF")
          ),
          bs4Card(title = "Genome Browser Plot", width = 6,
            plotOutput("genome_browser_plot", height = "600px")
          )
        )
      ),

      bs4TabItem(tabName = "structures",
        fluidRow(
          bs4Card(title = "sRNA Structure Selection", width = 12,
          helpText(HTML("<b>Note:</b> All secondary structures are generated and saved into ./www/ folder.")),
            selectInput("gene", "Select an image from below", choices = mylist),
            imageOutput("image", height = "600px")
          )
        )
      )
    )
  )
)


server <- function(input, output, session) {
  data_ready <- reactiveVal(FALSE)
  data_ready(FALSE)

  species_data <- reactive({
  if (input$species_choice == "Lactococcus lactis") {
    list(
      promoters = read.csv("lacto_promoters.csv"),
      terminators = read.csv("lacto_terminators.csv")
    )
  } else {
    list(
      promoters = read.csv("paracasei_promoters.csv"),
      terminators = read.csv("paracasei_terminators.csv")
    )
  }
})


  chr_id <- reactiveVal(NULL)
  chr_id2 <- reactiveVal(NULL)

  options(ucscChromosomeNames = FALSE)

 observeEvent(input$runAnalysis, {
    {
      session$sendCustomMessage("show_spinner", list())  # <- show spinner
  on.exit({
    session$sendCustomMessage("hide_spinner", list())  # <- hide spinner
  }, add = TRUE)
      req(input$bamFile, input$annotationFile)

      # Anotacijos failas
      ptt <- read.csv(input$annotationFile$datapath, sep = "\t", skip = 2, header = TRUE, stringsAsFactors = FALSE)

      # 5'detekcija
      res <- APERO_start_detection(
        work_dir = "./",
        bam_name = input$bamFile$datapath,
        ptt_file = ptt,
        wmax = input$wmax,
        min_dist = input$min_dist,
        enrichment = input$enrichment,
        min_read_number = input$min_read_number,
        genome_size = input$genome_size
      )

      # 3' detekcija
      res2 <- APERO_end_detection(
        work_dir = "./",
        start_table = res,
        mTEX_bam = input$bamFile$datapath,
        readthrough_proportion = input$readthrough_proportion,
        Fmin = input$Fmin,
        thread_number = input$thread_number,
        genome_size = input$genome_size,
        ptt_file = ptt
      )
      #padarau kad butu End column, kad nereiktu tolimesneje analizeje jo ideti ir iskart butu galima atsisiusti rezultatus kur yra ir sRNR galo koordinates
      res2$End <- res2$Position + res2$lg
      #padarau du naujus stulpelius del koordinaciu platinimo
      res2$Pos2 <- res2$Position - 100

      res2$End2<- res2$End + 100



observeEvent(input$bamFile, {
  req(input$bamFile)
  bamh <- BamFile(input$bamFile$datapath)
  header <- scanBamHeader(bamh)
  chrom_names <- names(header$targets)
  chr_id(chrom_names[1])  
})




      header <- scanBamHeader(BamFile(input$bamFile$datapath))
chrom_names <- names(header$targets)
chr_id(chrom_names[1])
res2$chr <- chr_id()

res2 <- res2[!is.na(res2$str), ]


if (!is.null(res2$str) && length(res2$str) == nrow(res2)) {
  granges_platinti <- GRanges(
    seqnames = rep(res2$chr[1], nrow(res2)),
    ranges = IRanges(start = res2$Pos2, end = res2$End2),
    strand = Rle(res2$str)
  )
} else {
  showNotification("Mismatch between strand and coordinate rows. GRanges not created.", type = "error")
  return(NULL)
}



if (!is.null(res2$str) && length(res2$str) == nrow(res2)) {
  granges_neplatinti <- GRanges(
    seqnames = rep(res2$chr[1], nrow(res2)),
    ranges = IRanges(start = res2$Position, end = res2$End),
    strand = Rle(res2$str)
  )
} else {
  showNotification("Mismatch between strand and coordinate rows. GRanges not created.", type = "error")
  return(NULL)
}



      bampath <- input$bamFile$datapath
    indexBam(bampath)
      bamFile <- BamFile(bampath)
genes <- granges_platinti

covSigs <- bamCoverage(bampath, genes, verbose=FALSE)

reiksmes <- as.list(covSigs)[[1]]


reiksmes <- reiksmes + 1
pradzios <- start(granges_platinti)

analyze_segment <- function(reiksmes, pradzia, grandine) {
    start_position <- NA
    end_position <- NA
    first_peak_found <- FALSE

   
    for (i in 1:(length(reiksmes) - 2)) {
        tikrinamas_ivertis <- reiksmes[i]
        sekantis_ivertis <- reiksmes[i + 1]
        trecias_ivertis <- reiksmes[i + 2]

        if (!first_peak_found) {
            if ((sekantis_ivertis >= 3 * tikrinamas_ivertis & sekantis_ivertis > 5) ||
                (trecias_ivertis >= 3 * tikrinamas_ivertis & trecias_ivertis > 5)) {
                start_position <- i + 1
                first_peak_found <- TRUE
                break
            }
        }
    }

   
    for (i in length(reiksmes):1) {
        tikrinamas_ivertis <- reiksmes[i]
        sekantis_ivertis <- ifelse(i + 1 <= length(reiksmes), reiksmes[i + 1], NA)
        trecias_ivertis <- ifelse(i + 2 <= length(reiksmes), reiksmes[i + 2], NA)

        if ((!is.na(sekantis_ivertis) && sekantis_ivertis <= tikrinamas_ivertis / 3 & sekantis_ivertis > 5) ||
            (!is.na(trecias_ivertis) && trecias_ivertis <= tikrinamas_ivertis / 3 & trecias_ivertis > 5)) {
            end_position <- i
            break
        }
    }
    if (!is.na(start_position) && !is.na(end_position)) {
        adjusted_start_position <- start_position + pradzia - 1
        adjusted_end_position <- end_position + pradzia - 1
        return(data.frame(start_position = adjusted_start_position, end_position = adjusted_end_position, strand = grandine))
    } else {
        return(data.frame(start_position = NA, end_position = NA, strand = NA))
    }
}


rezultatai <- lapply(seq_along(covSigs), function(i) {
  reiksmes <- as.list(covSigs)[[i]]
  grandine <- as.character(strand(granges_platinti)[i])
  analyze_segment(reiksmes, pradzios[i], grandine)
})
unikalus <- do.call(rbind, rezultatai) %>%
  distinct(start_position, end_position, strand, .keep_all = TRUE) %>%
  filter(!is.na(start_position) & !is.na(end_position))


rownames(unikalus) <- NULL
unikalus_filtered <- unikalus %>%
  filter(start_position < end_position & end_position-start_position <= 500)


observeEvent(input$fasta_file, {
  fasta <- readDNAStringSet(input$fasta_file$datapath)
  if (length(names(fasta)) > 1) {
    showNotification("Warning: Multiple sequences detected in FASTA. Only the first will be used.", type = "warning")
  }
  chr_id(names(fasta)[1])
})


if (!is.null(unikalus_filtered$strand) && length(unikalus_filtered$strand) == nrow(unikalus_filtered)) {
naujas_granges <- GRanges(
  seqnames = rep(chr_id(), nrow(unikalus_filtered)),
  ranges = IRanges(start = unikalus_filtered$start_position, end = unikalus_filtered$end_position),
  strand = unikalus_filtered$strand
)
} else {
  showNotification("Mismatch between strand and coordinate rows. GRanges not created.", type = "error")
  return(NULL)
}


df_koreguotas <- data.frame(seqnames=seqnames(naujas_granges),
  starts=start(naujas_granges)-1,
  ends=end(naujas_granges),
  strands=strand(naujas_granges))


df_koreguotas$length= df_koreguotas$ends - df_koreguotas$starts

df_filtruotas<-filter(df_koreguotas, length <= 500)

df_atrinktas = subset(df_filtruotas, select = -length)



#filtruoju pagal iRNR ir tRNR ir rRNR ir genus

if (!is.null(input$gff3_path) && file.exists(input$gff3_path$datapath)) {
  gff3_data <- import(input$gff3_path$datapath, format = "GFF3")
} else {
  showNotification("Please upload a valid GFF3 file.", type = "error")
  return(NULL)
}





annotations <- gff3_data
if (!is.null(df_atrinktas$strands) && length(df_atrinktas$strands) == nrow(df_atrinktas)) {
gr <- GRanges(
  seqnames = df_atrinktas$seqnames,
  ranges = IRanges(start = df_atrinktas$starts, end = df_atrinktas$ends),
  strand = df_atrinktas$strands
)
} else {
  showNotification("Mismatch between strand and coordinate rows. GRanges not created.", type = "error")
  return(NULL)
}



tRNAs <- subset(annotations, type == "tRNA")
mRNAs <- subset(annotations, type == "CDS")
genes<-subset(annotations, type == "gene")
rRNAs<-subset(annotations, type == "rRNA")


overlap_fun <- function(x, y, threshold) {
  overlaps <- findOverlaps(x, y)
  overlap_widths <- width(pintersect(x[queryHits(overlaps)], y[subjectHits(overlaps)]))
  sRNA_widths <- width(x[queryHits(overlaps)])
  percent_overlap <- (overlap_widths / sRNA_widths) * 100
  return(queryHits(overlaps)[percent_overlap >= threshold])
}


overlapping_indices_tRNAs <- overlap_fun(gr, tRNAs, input$tRNA_thresh)
overlapping_indices_rRNAs <- overlap_fun(gr, rRNAs, input$rRNA_thresh)
overlapping_indices_genes <- overlap_fun(gr, genes, input$gene_thresh)
overlapping_indices_mRNAs <- overlap_fun(gr, mRNAs, input$mRNA_thresh)


overlapping_indices <- unique(c(overlapping_indices_tRNAs, overlapping_indices_rRNAs, overlapping_indices_genes))


non_overlapping_sRNAs <- gr[-overlapping_indices]
df_atrinktas <- as.data.frame(non_overlapping_sRNAs)
colnames(df_atrinktas)[colnames(df_atrinktas) == "start"] <- "starts"
colnames(df_atrinktas)[colnames(df_atrinktas) == "end"] <- "ends"
colnames(df_atrinktas)[colnames(df_atrinktas) == "strand"] <- "strands"

df_atrinktas$width <- df_atrinktas$ends - df_atrinktas$starts


colnames(df_atrinktas)<-c("seqnames", "starts" ,  "ends" ,  "width",  "strands")





      # Atsisiust faila
      output$downloadResults <- downloadHandler(
        filename = "apero_results.csv",
        content = function(file) {
          write.csv(df_atrinktas, file, row.names = FALSE)
        }
      )
    }
  })



#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////












files_and_labels <- reactive({
  list(
    list(files = input$group1_files, label = input$group1_label),
    list(files = input$group2_files, label = input$group2_label)
  )
})
observeEvent(input$process_btn, {


req(input$input_option == "APERO results analysis (Second step)")      
  req(input$group1_files, input$group2_files)      
  req(input$fasta_file)                            
 


  #combininu failus ir abieju grupiu
  files_info <- rbind(
    data.frame(group = rep(input$group1_label, length(input$group1_files$name)),
               name = input$group1_files$name,
               size = input$group1_files$size,
               type = input$group1_files$type,
               sRNA_count = sapply(input$group1_files$datapath, function(path) {
                 if (!is.null(path)) {
                   tryCatch({
                     nrow(read.csv(path))
                   }, error = function(e) { NA })
                 } else {
                   NA
                 }
               })),
    data.frame(group = rep(input$group2_label, length(input$group2_files$name)),
               name = input$group2_files$name,
               size = input$group2_files$size,
               type = input$group2_files$type,
               sRNA_count = sapply(input$group2_files$datapath, function(path) {
                 if (!is.null(path)) {
                   tryCatch({
                     nrow(read.csv(path))
                   }, error = function(e) { NA })
                 } else {
                   NA
                 }
               }))
  )


 
 
 


  output$filesTable <- renderDataTable({
    files_info
  }, options = list(pageLength = 5, autoWidth = TRUE))
})


extract_srna_sequences <- function(files_paths, file_names, fasta_file_path, label) {





    fasta_sequences <- readDNAStringSet(fasta_file_path$datapath)




    result_df <- data.frame(File = character(), Type = character(), Sequence = character(),
                            Position = integer(), End = integer(), ID = character(), Strand = character(), Length = integer(), GC_content = integer(),
                            
                            stringsAsFactors = FALSE)

    for (i in 1:length(files_paths)) {
        file_path <- files_paths[i]



        if (!is.null(file_path) && file.exists(file_path)) {
           
            data <- read.csv(file_path, stringsAsFactors = FALSE)





            for (j in 1:nrow(data)) {
               
                start_pos <- data$starts[j]
                end_pos <- data$ends[j]
                extracted_seq <- subseq(fasta_sequences, start_pos, end_pos)
                chromosome_id <- data$seqnames[j]
                strands <- data$strands[j]
                file_name <- file_names[i]
                file_type <- label
                sequence_length <- width(extracted_seq)
                gc_count <- sum(letterFrequency(extracted_seq, letters = c("G", "C"), as.prob = TRUE))
                gc_content <- 100 * gc_count / sequence_length
                
                ID = NA
                extracted_df <- data.frame(
                  File = file_name,
                  Type = file_type,
                  Start = start_pos,
                  End = end_pos,
                  Chromosome = chromosome_id,
                  ID = NA,
                  Strand = strands,
                  Length = sequence_length,
                  GC_content = gc_content,
                                          
                                            Sequence = as.character(extracted_seq),
                                           stringsAsFactors = FALSE)
                result_df <- rbind(result_df, extracted_df)
            }
        }
    }
    return(result_df)
}
 




 processed_sequences <- eventReactive(input$process_btn, {
  session$sendCustomMessage("show_spinner", list())
    tryCatch({
    req(input$fasta_file)

   
   files_labels <- files_and_labels()
files_labels <- files_and_labels()
group1_files <- files_labels[[1]]$files$datapath
group1_names <- files_labels[[1]]$files$name
group2_files <- files_labels[[2]]$files$datapath
group2_names <- files_labels[[2]]$files$name

extracted_sequences_group1 <- extract_srna_sequences(group1_files, group1_names, input$fasta_file, files_labels[[1]]$label)
extracted_sequences_group2 <- extract_srna_sequences(group2_files, group2_names, input$fasta_file, files_labels[[2]]$label)

extracted_sequences <- rbind(extracted_sequences_group1, extracted_sequences_group2)





  # kad butu galima prafiltruoti pagal ilgius
  min_length <- input$min_length
  max_length <- input$max_length
  extracted_sequences <- extracted_sequences[nchar(extracted_sequences$Sequence) >= min_length & nchar(extracted_sequences$Sequence) <= max_length, ]



sRNAs <- extracted_sequences
sRNAs<-sRNAs[ ,c(1,2,3,4,5,6,7,8,9,10)]
genes <- read.delim(input$file_ptt$datapath, skip=2, header=TRUE) %>%
  separate(Location, into = c("gene_start", "gene_end"), sep = "\\.\\.", convert = TRUE)





#pirmiausia surandu genu overlapus abiejuose stranduose
sRNAs <- sRNAs %>%
  rowwise() %>%
  mutate(
    overlapping_same = list(which(genes$gene_start <= End & genes$gene_end >= Start & genes$Strand == Strand)),
    overlapping_opposite = list(which(genes$gene_start <= End & genes$gene_end >= Start & genes$Strand != Strand)),

    overlap_same_idx = ifelse(length(pluck(overlapping_same, 1, .default = integer(0))) > 0,
                              pluck(overlapping_same, 1, which.max(
                                pmin(genes$gene_end[pluck(overlapping_same, 1)], End) -
                                  pmax(genes$gene_start[pluck(overlapping_same, 1)], Start) + 1
                              ), .default = NA_integer_),
                              NA_integer_),

    overlap_opposite_idx = ifelse(length(pluck(overlapping_opposite, 1, .default = integer(0))) > 0,
                                  pluck(overlapping_opposite, 1, which.max(
                                    pmin(genes$gene_end[pluck(overlapping_opposite, 1)], End) -
                                      pmax(genes$gene_start[pluck(overlapping_opposite, 1)], Start) + 1
                                  ), .default = NA_integer_),
                                  NA_integer_),

   
    Group = case_when(
      !is.na(overlap_same_idx) & (Start < genes$gene_start[overlap_same_idx] | End > genes$gene_end[overlap_same_idx]) ~ "UTR-derived",
      !is.na(overlap_same_idx) ~ "Sense overlap",
      is.na(overlap_same_idx) & !is.na(overlap_opposite_idx) ~ "Antisense",
      TRUE ~ "Intergenic"
    ),

    overlap_same_id = ifelse(!is.na(overlap_same_idx), genes$Synonym[overlap_same_idx], NA_character_),
    overlap_same_product = ifelse(!is.na(overlap_same_idx), genes$Product[overlap_same_idx], NA_character_),
    overlap_opposite_id = ifelse(!is.na(overlap_opposite_idx), genes$Synonym[overlap_opposite_idx], NA_character_),
    overlap_opposite_product = ifelse(!is.na(overlap_opposite_idx), genes$Product[overlap_opposite_idx], NA_character_)
  ) %>%
  ungroup()


assign_group <- function(Start, End, Strand, overlap_same_idx, overlap_opposite_idx, genes) {
  groups <- c()
  if (!is.na(overlap_same_idx)) {
   
    groups <- c(groups, "Sense overlap")
    gene_strand <- genes$Strand[overlap_same_idx]
    gene_start  <- genes$gene_start[overlap_same_idx]
    gene_end    <- genes$gene_end[overlap_same_idx]
    if (gene_strand == "+") {
      if (Start < gene_start) groups <- c(groups, "5' UTR-derived")
      if (End > gene_end)   groups <- c(groups, "3' UTR-derived")
    } else {
      if (End > gene_end)   groups <- c(groups, "5' UTR-derived")
      if (Start < gene_start) groups <- c(groups, "3' UTR-derived")
    }
  }
  if (is.na(overlap_same_idx) & !is.na(overlap_opposite_idx)) {
    groups <- c(groups, "Antisense")
  }
  if (is.na(overlap_same_idx) & is.na(overlap_opposite_idx)) {
    groups <- c(groups, "Intergenic")
  }
  paste(unique(groups), collapse = ", ")
}

# antra- surandu artimiausia gena
sRNAs <- sRNAs %>%
  rowwise() %>%
  mutate(
    # artimiausias upstream (5') gene:
    upstream_idx = {
      idx <- which(genes$gene_end < Start)
      if(length(idx) == 0) NA_integer_ else max(idx)
    },
    upstream_gene = ifelse(!is.na(upstream_idx), genes$Synonym[upstream_idx], NA_character_),
    upstream_product = ifelse(!is.na(upstream_idx), genes$Product[upstream_idx], NA_character_),
   
    # artimiausias downstream (3') gene:
    downstream_idx = {
      idx <- which(genes$gene_start > End)
      if(length(idx) == 0) NA_integer_ else min(idx)
    },
    downstream_gene = ifelse(!is.na(downstream_idx), genes$Synonym[downstream_idx], NA_character_),
    downstream_product = ifelse(!is.na(downstream_idx), genes$Product[downstream_idx], NA_character_),
   
    # genomic context
    # primas + arba -→ grandinė sRNA, jei overlapina tanme pačiame strande (overlap_same_idx). jei nėra overlapo- automatiškai +, antras + arba -: The strand of the overlapping gene on the same strand (using genes$Strand[overlap_same_idx]). jei ne- default "+". trecias + arba -: grandin4 geno kitoje grandin4je  (genes$Strand[overlap_opposite_idx]). jei nera geno- "+".
    genomic_context = paste0(
      ifelse(!is.na(overlap_same_idx), ifelse(Strand == "+", "+", "-"), "+"), "/",
      ifelse(!is.na(overlap_same_idx), ifelse(genes$Strand[overlap_same_idx] == "+", "+", "-"), "+"), "/",
      ifelse(!is.na(overlap_opposite_idx), ifelse(genes$Strand[overlap_opposite_idx] == "+", "+", "-"), "+")
    ),
   
   
    Group = assign_group(Start, End, Strand, overlap_same_idx, overlap_opposite_idx, genes)
  ) %>%
  ungroup()

updated_my_data <- sRNAs

species_code <- if (input$species_choice == "Lactococcus lactis") "LBL" else "LCB"
updated_my_data <- updated_my_data %>%
  mutate(
    bin = sprintf("%04d", Start %/% 100),
    base_id = paste0("s", species_code, bin)
  ) %>%
  group_by(base_id, Start, End, Strand) %>%
  mutate(dup_count = n()) %>%
  ungroup() %>%
  group_by(base_id, Strand) %>%
  mutate(letter = if (n() > 1) letters[row_number()] else "") %>%
  ungroup() %>%
  mutate(ID = paste0(base_id, letter, Strand)) %>%
  select(-bin, -base_id, -dup_count, -letter)

updated_my_data_cleaned <- updated_my_data %>%
  select(File, Type, Start, End, Chromosome, ID, Strand, Group, Length, GC_content, Sequence)



sRNA <- updated_my_data_cleaned





final_table <- updated_my_data %>%
  transmute(
    ID = ID,
    Group = Group,
    Strand = Strand,
    Start = Start,
    End = End,
    Genomic_Context = genomic_context,
    `5' Flanking Gene ID` = upstream_gene,
    `5' Flanking Gene Product` = upstream_product,
    `3' Flanking Gene ID` = downstream_gene,
    `3' Flanking Gene Product` = downstream_product,
    `Sense Gene ID` = overlap_same_id,
    `Sense Gene Product` = overlap_same_product,
    `Gene at 5' or 3' End ID` = overlap_opposite_id,
    `Gene at 5' or 3' End Product` = overlap_opposite_product
  )




genomic_context_lentele <- final_table 




Promoters <- species_data()$promoters
Terminators <- species_data()$terminators

merged_promoters<-Promoters



sRNAs<-sRNA
max_gap <- input$max_promoter_distance

results_list <- list()

for(i in seq_len(nrow(sRNAs))) {
  srna <- sRNAs[i, ]
 
  if(srna$Strand == "+") {
   
    matches <- merged_promoters %>%
      filter((srna$Start - minus10_end) >= 0,
             (srna$Start - minus10_end) <= max_gap)
   
    if(nrow(matches) > 0) {
      m <- matches[1, ]  
      gap <- srna$Start - m$minus10_end
      promoter_presence <- "+"
    } else {
      m <- data.frame(minus10_start = NA, minus10_end = NA, minus10_seq = NA,
                      minus35_start = NA, minus35_end = NA, minus35_seq = NA,
                      promoter_pos = NA)
      gap <- NA
      promoter_presence <- "-"
    }
   
  } else if(srna$Strand == "-") {
    matches <- merged_promoters %>%
      filter((minus10_start - srna$End) >= 0,
             (minus10_start - srna$End) <= max_gap)
   
    if(nrow(matches) > 0) {
      m <- matches[1, ]
      gap <- m$minus10_start - srna$End
      promoter_presence <- "+"
    } else {
      m <- data.frame(minus10_start = NA, minus10_end = NA, minus10_seq = NA,
                      minus35_start = NA, minus35_end = NA, minus35_seq = NA,
                      promoter_pos = NA)
      gap <- NA
      promoter_presence <- "-"
    }
   
  } else {
    m <- data.frame(minus10_start = NA, minus10_end = NA, minus10_seq = NA,
                    minus35_start = NA, minus35_end = NA, minus35_seq = NA,
                    promoter_pos = NA)
    gap <- NA
    promoter_presence <- "-"
  }
 
  out <- data.frame(
    ID             = srna$ID,
    Strand         = srna$Strand,
    Group          = srna$Group,
    Start          = srna$Start,
    End            = srna$End,
    Promoter       = promoter_presence,  
    Distance       = gap,
    `Minus10 start`  = m$minus10_start,
    `Minus10 end`   = m$minus10_end,
    `Minus10 seq`    = m$minus10_seq,
    `Minus35 start`  = m$minus35_start,
    `Minus35 end`    = m$minus35_end,
    `Minus35 seq`    = m$minus35_seq,
    stringsAsFactors = FALSE
  )
 
  results_list[[i]] <- out
}

final_df <- do.call(rbind, results_list)
su_promotoriais<-final_df

sRNAs_with_promoter <- final_df[final_df$Promoter == "+", ]




terminators<-Terminators



find_sRNA_terminators <- function(sRNAs, terminators, max_distance = input$max_terminator_distance) {
  sRNAs <- sRNAs %>%
    mutate(Lacto_Terminator = "-")  
 
  for (i in 1:nrow(sRNAs)) {
    current_sRNA <- sRNAs[i, ]
   
    relevant_terms <- terminators %>% filter(Chain == current_sRNA$Strand)
   
    if (current_sRNA$Strand == "+") {
      close_terms <- relevant_terms %>%
        filter(Start >= current_sRNA$End & Start <= current_sRNA$End + max_distance)
    } else {
      close_terms <- relevant_terms %>%
        filter(Start <= current_sRNA$Start & Start >= current_sRNA$Start - max_distance)
    }
   
    if (nrow(close_terms) > 0) {
      sRNAs$Lacto_Terminator[i] <- "+"
    }
  }
 
  return(sRNAs)
}

# darau analize
result_df <- find_sRNA_terminators(sRNAs, terminators, max_distance = 20)

su_terminatoriais<-result_df




sRNAs_df <- sRNA
lacto_terminators_df <- terminators


matched_sRNAs <- data.frame()

for (i in 1:nrow(sRNAs_df)) {
  sRNA <- sRNAs_df[i, ]
 
  relevant_terms <- lacto_terminators_df %>%
    filter(Chain == sRNA$Strand)
 
  if (sRNA$Strand == "+") {
    close_terms <- relevant_terms %>%
      filter(Start >= sRNA$End & Start <= sRNA$End + 20)
  } else {
    close_terms <- relevant_terms %>%
      filter(Start <= sRNA$Start & Start >= sRNA$Start - 20)
  }
 
  if (nrow(close_terms) > 0) {
    for (j in 1:nrow(close_terms)) {
      term <- close_terms[j, ]
     
      matched_sRNAs <- rbind(matched_sRNAs, data.frame(
        ID = sRNA$ID,
        Strand = sRNA$Strand,
        `sRNA Start` = sRNA$Start,
        `sRNA End` = sRNA$End,
        `sRNA Group` = sRNA$Group,
        Terminator = "+",
        `Terminator Start` = term$Start,
        `Terminator End` = term$End,
        `Terminator Sequence` = term$Sequence
      ))
    }
  }
}


sRNAs_with_terminators<-matched_sRNAs

sRNAs$Promoter <- su_promotoriais$Promoter

sRNAs$Terminator <- su_terminatoriais$Lacto_Terminator

# ///////////////////// blastas

xlsx_data <- sRNAs


    xlsx_data$GC_content <- as.numeric(str_replace(xlsx_data$GC_content, ",", "."))

 
    results <- data.frame(
      Name = xlsx_data$ID,
      Group = xlsx_data$Group,
      Strand = xlsx_data$Strand,
      Start = xlsx_data$Start,
      End = xlsx_data$End,
      stringsAsFactors = FALSE
    )

   
    query_file <- tempfile(fileext = ".fasta")
    sequences <- DNAStringSet(xlsx_data$Sequence)
    names(sequences) <- xlsx_data$ID
    writeXStringSet(sequences, filepath = query_file)

    results_blast <- results  

for (fasta in input$fasta_files$datapath) {
 
  system(paste("makeblastdb -in", shQuote(fasta), "-dbtype nucl"))

 
  blast_out <- tempfile(fileext = ".tsv")
  cmd <- paste("blastn -query", shQuote(query_file),
               "-db", shQuote(fasta),
               "-outfmt '6 qseqid pident'",
               "-max_target_seqs 1 -evalue 1e-5 -out", shQuote(blast_out))
  system(cmd)


  if (file.info(blast_out)$size > 0) {
    blast_result <- fread(blast_out, header = FALSE)
    fasta_header <- names(readDNAStringSet(fasta, nrec = 1))[[1]]
    clean_name <- gsub("[^A-Za-z0-9_]+", "_", fasta_header)
    colnames(blast_result) <- c("ID", clean_name)
  } else {
    blast_result <- data.frame(ID = xlsx_data$ID, temp_col = NA)
    colnames(blast_result)[2] <- basename(fasta)
  }


  blast_result <- blast_result[!duplicated(blast_result$ID), ]


  results_blast <- merge(results_blast, blast_result, by.x = "Name", by.y = "ID", all.x = TRUE)
}



colnames(matched_sRNAs) <- c("ID", "Strand", "Start", "End", "Group", "Terminator", "Terminator start", "Terminator end", "Terminator Sequence")

  return(list(
    sRNAs = sRNAs,
    sRNAs_with_promoter = sRNAs_with_promoter,
    sRNAs_with_terminators = matched_sRNAs,
    sRNAs_genomic_context = genomic_context_lentele,
    sRNAs_blast = results_blast
  ))
}, error = function(e) {
    
    print(e$message)
    NULL  
  })
 })


output$lentele_table <- renderDT({
  datatable(processed_sequences()$sRNAs, options = list(scrollX = TRUE))
})

output$promoter_table <- renderDT({
  datatable(
    processed_sequences()$sRNAs_with_promoter %>%
      distinct(Start, End, .keep_all = TRUE),
    options = list(scrollX = TRUE)
  )
})

output$terminators_table <- renderDT({
  datatable(
    processed_sequences()$sRNAs_with_terminators %>%
      distinct(Start, End, .keep_all = TRUE),
    options = list(scrollX = TRUE)
  )
})

output$genomic_context <- renderDT({
  datatable(
    processed_sequences()$sRNAs_genomic_context %>%
      distinct(Start, End, .keep_all = TRUE),
    options = list(scrollX = TRUE)
  )
})

output$blast_results <- renderDT({
  datatable(
    processed_sequences()$sRNAs_blast %>%
      distinct(Start, End, .keep_all = TRUE),
    options = list(scrollX = TRUE)
  ) %>%
    formatStyle(
      columns = 6:ncol(processed_sequences()$sRNAs_blast),
      background = styleColorBar(range(0, 100, na.rm = TRUE), '#059605'),
      backgroundSize = '100% 90%',
      backgroundRepeat = 'no-repeat',
      backgroundPosition = 'center'
    )
})





 
  output$download_all <- downloadHandler(
  filename = "All_Tables.xlsx",
  content = function(file){
    wb <- createWorkbook()

    addWorksheet(wb, "Main information")
    writeData(wb, "Main information",
      processed_sequences()$sRNAs %>%
        distinct(Start, End, .keep_all = TRUE)
    )

    addWorksheet(wb, "Promoter information")
    writeData(wb, "Promoter information",
      processed_sequences()$sRNAs_with_promoter %>%
        distinct(Start, End, .keep_all = TRUE)
    )

    addWorksheet(wb, "Terminator information")
    writeData(wb, "Terminator information",
      processed_sequences()$sRNAs_with_terminators %>%
        distinct(Start, End, .keep_all = TRUE)
    )

    addWorksheet(wb, "Genomic context information")
    writeData(wb, "Genomic context information",
      processed_sequences()$sRNAs_genomic_context %>%
        distinct(Start, End, .keep_all = TRUE)
    )

    addWorksheet(wb, "Structure information")
    writeData(wb, "Structure information",
      structures_data() %>%
        distinct(Start, End, .keep_all = TRUE)
    )

    blast_data <- processed_sequences()$sRNAs_blast %>%
      distinct(Start, End, .keep_all = TRUE)
    addWorksheet(wb, "BLAST Results")
    writeData(wb, "BLAST Results", blast_data)

    if (ncol(blast_data) > 5) {
      blast_cols <- 6:ncol(blast_data)

      rotated_style <- createStyle(textRotation = 45, halign = "center", valign = "bottom")

      addStyle(
        wb, "BLAST Results",
        style = rotated_style,
        cols = blast_cols,
        rows = 1,
        gridExpand = TRUE
      )

      setColWidths(wb, "BLAST Results", cols = blast_cols, widths = 15)

      conditionalFormatting(
        wb, "BLAST Results",
        cols = blast_cols, rows = 2:(nrow(blast_data) + 1),
        style = c("#ffffff", "#0e870e"),
        type = "colorScale"
      )
    }

    saveWorkbook(wb, file, overwrite = TRUE)
  }
)



observeEvent(processed_sequences(), {
  main_info <- processed_sequences()$sRNAs
  updateNumericInput(session, "region_start", value = min(main_info$Start, na.rm = TRUE))
  updateNumericInput(session, "region_end", value = min(main_info$End, na.rm = TRUE))
})


observeEvent(processed_sequences(), {
  main_info <- processed_sequences()$sRNAs
  updateSelectInput(session, "sRNA_choice", choices = main_info$ID)
})


observeEvent(processed_sequences(), {
  main_info <- processed_sequences()$sRNAs
  groups <- unique(main_info$Group)
  updateSelectizeInput(session, "group_choice", choices = groups, selected = groups)
})



(ucscChromosomeNames=FALSE)
process_data <- reactive({
    req(processed_sequences())
    options(ucscChromosomeNames=FALSE)

    sRNAs <- processed_sequences()$sRNAs
    promoters <- processed_sequences()$sRNAs_with_promoter
    terminators <- processed_sequences()$sRNAs_with_terminators

chrom <- chr_id()
if (is.null(chrom)) {
  warning("chr_id() returned NULL. Defaulting chromosome to 'unknown'.")
  chrom <- "unknown"
}

sRNAs$Chromosome <- rep(chrom, nrow(sRNAs))
promoters$Chromosome <- rep(chrom, nrow(promoters))
terminators$Chromosome <- rep(chrom, nrow(terminators))

   

    genes <- NULL
    if (!is.null(input$file_ptt)) {
      genes <- read.delim(input$file_ptt$datapath, skip = 2, header = TRUE)


    }

    list(
      sRNAs = sRNAs,
      promoters = promoters,
      terminators = terminators,
      genes = genes
    )
})



genome_browser_plot_obj <- reactive({
  data <- process_data()
  sRNAs <- data$sRNAs
  promoters <- data$promoters
  terminators <- data$terminators
  genes <- data$genes

  correct_chr <- chr_id()

  sRNA_gr <- GRanges(
      seqnames = sRNAs$Chromosome,
      ranges = IRanges(start = sRNAs$Start, end = sRNAs$End),
      strand = sRNAs$Strand,
      Group = sRNAs$Group
    )

    promoter_gr <- GRanges(
      seqnames = promoters$Chromosome,
      ranges = IRanges(
        start = pmin(promoters$`Minus35.start`, promoters$`Minus10.end`),
        end = pmax(promoters$`Minus35.start`, promoters$`Minus10.end`)
      ),
      strand = promoters$Strand
    )

    terminator_gr <- GRanges(
      seqnames = terminators$Chromosome,
      ranges = IRanges(start = terminators$`Terminator start`, end = terminators$`Terminator end`),
      strand = terminators$Strand
    )

    gtrack <- GenomeAxisTrack()

    sRNA_track <- AnnotationTrack(
      sRNA_gr, name = "sRNAs", shape = "arrow", stacking = "dense",
      group = sRNA_gr$Group, feature = sRNA_gr$Group,
      fill = "lightblue"
    )

    promoter_track <- AnnotationTrack(
      promoter_gr, name = "Promoters", shape = "box", stacking = "dense", fill = "green"
    )

    terminator_track <- AnnotationTrack(
      terminator_gr, name = "Terminators", shape = "box", stacking = "dense", fill = "red"
    )

    tracks <- list(gtrack)

   
    gene_track <- NULL
    if (!is.null(genes)) {
     
      gene_gr <- GRanges(
        seqnames = rep(unique(sRNAs$Chromosome)[1], nrow(genes)),
        ranges = IRanges(start = as.integer(sub("\\..*", "", genes$Location)),
                         end = as.integer(sub(".*\\.\\.", "", genes$Location))),
        strand = genes$Strand,
        gene = genes$Synonym,
        product = genes$Product
      )

      if (!is.null(gene_gr)) {
  gene_colors <- ifelse(as.character(strand(gene_gr)) == "+", "dodgerblue", "firebrick")

  gene_track <- GeneRegionTrack(
    gene_gr,
    genome = "unknown",
    chromosome = as.character(unique(seqnames(gene_gr))),
    name = "Genes",
    transcriptAnnotation = "gene",
    showId = TRUE,
    col = "black",
    fill = gene_colors
  )
}

    }

 
    if (!is.null(gene_track)) {
      tracks <- c(tracks, list(gene_track))
    }

    tracks <- c(tracks, list(sRNA_track, promoter_track, terminator_track))

    tracks
  })



output$genome_browser_plot <- renderPlot({
  req(genome_browser_plot_obj())

  validate(
    need(input$region_start < input$region_end, "Start coordinate must be less than end coordinate."),
    need(input$region_start >= 0, "Start coordinate must be positive."),
    need(input$region_end >= 0, "End coordinate must be positive.")
  )

 

  plotTracks(
    genome_browser_plot_obj(),
    from = input$region_start,
    to = input$region_end,
    background.title = "white",
    col.title = "black",
    cex.title = 1.2,
    cex.axis = 0.9
  )
})



output$download_genome_browser <- downloadHandler(
  filename = function() {
    paste0("Genome_Browser_", Sys.Date(), ".pdf")
  },
  content = function(file) {
    pdf(file, width = 12, height = 8)
    plotTracks(
      genome_browser_plot_obj(),
      from = input$region_start,
      to = input$region_end,
      background.title = "white",
      col.title = "black",
      cex.title = 1.2,
      cex.axis = 0.9
    )
    dev.off()
  }
)



output$venn_plot <- renderPlot({
  req(processed_sequences()$sRNAs)

  
  sRNAs <- processed_sequences()$sRNAs

  
  set_list <- split(sRNAs$Sequence, sRNAs$Type)
  set_list <- lapply(set_list, unique)

  
  n_sets <- length(set_list)
  if (n_sets < 2) {
    showNotification("At least two groups are required for a Venn diagram.", type = "error")
    return(NULL)
  }

  
  set_colors <- if (n_sets >= 3) {
    RColorBrewer::brewer.pal(min(n_sets, 8), "Set3")
  } else {
    c("#66c2a5", "#fc8d62")[1:n_sets]
  }

  grid::grid.newpage()

  if (n_sets == 2) {
    venn_plot <- VennDiagram::draw.pairwise.venn(
      area1 = length(set_list[[1]]),
      area2 = length(set_list[[2]]),
      cross.area = length(intersect(set_list[[1]], set_list[[2]])),
      category = names(set_list),
      fill = set_colors,
      cex = 1,
      cat.cex = 1,
      cat.col = "black",
      cat.fontface = "bold"
    )
    grid::grid.draw(venn_plot)
  } else if (n_sets <= 5) {
    venn_plot <- VennDiagram::venn.diagram(
      x = set_list,
      category.names = names(set_list),
      filename = NULL,
      fill = set_colors,
      cex = 1,
      cat.cex = 1,
      cat.col = "black",
      cat.fontface = "bold",
      main = "Overlap of sRNAs Between Files",
      main.cex = 2,
      margin = 0.05,
      alpha = 0.7
    )
    grid::grid.draw(venn_plot)
  } else {
    showNotification("Too many sets for Venn plot. Use UpSet plot instead.", type = "warning")
  }
})




output$upset_plot <- renderPlot({
  req(processed_sequences()$sRNAs)
  sRNAs <- processed_sequences()$sRNAs

  id_by_file <- sRNAs %>%
    dplyr::select(Type, Sequence) %>%
    distinct() %>%
    mutate(present = 1) %>%
    tidyr::pivot_wider(names_from = Type, values_from = present, values_fill = list(present = 0))

  rownames(id_by_file) <- id_by_file$Sequence
  id_by_file <- id_by_file[, -1, drop = FALSE]

 
  UpSetR::upset(
  as.data.frame(id_by_file),
  sets = colnames(id_by_file),
  order.by = "freq",
  keep.order = TRUE,
  sets.bar.color = "#7b68ee",
  main.bar.color = "#4169e1",
  matrix.color = "#6a5acd",
  text.scale = c(2, 2, 1.5, 1.5, 2, 2),
  mb.ratio = c(0.6, 0.4)
)
})


output$promoter_logo_plot <- renderPlot({
  req(processed_sequences())

  promoters <- processed_sequences()$sRNAs_with_promoter

  if (input$promoter_region == "Minus10") {
    seqs <- na.omit(promoters$Minus10.seq)
    region_label <- "–10 promoter region"
  } else {
    seqs <- na.omit(promoters$Minus35.seq)
    region_label <- "–35 promoter region"
  }

  seqs <- seqs[nchar(seqs) > 0 & grepl("^[ACGTacgt]+$", seqs)]

  if (length(seqs) == 0) {
    showNotification("No valid sequences available for plotting.", type = "warning")
    return(NULL)
  }

  ggseqlogo::ggseqlogo(seqs, method = "prob") +
    ggtitle(paste("Sequence Logo for", region_label)) +
    labs(subtitle = paste("Number of sequences:", length(seqs))) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    )
})

output$sRNA_count <- renderValueBox({
  req(processed_sequences()$sRNAs)
  bs4ValueBox(
    value = nrow(processed_sequences()$sRNAs),
    subtitle = "Total sRNAs Detected",
    icon = icon("dna"),
    color = "indigo"
  )
})

output$promoter_count <- renderValueBox({
  req(processed_sequences()$sRNAs_with_promoter)
  bs4ValueBox(
    value = nrow(processed_sequences()$sRNAs_with_promoter),
    subtitle = "With Promoter Detected",
    icon = icon("arrow-up"),
    color = "indigo"
  )
})

output$terminator_count <- renderValueBox({
  req(processed_sequences()$sRNAs_with_terminators)
  bs4ValueBox(
    value = nrow(processed_sequences()$sRNAs_with_terminators),
    subtitle = "With Terminator Detected",
    icon = icon("arrow-down"),
    color = "indigo"
  )
})



output$length_histogram <- renderPlot({
  req(processed_sequences()$sRNAs)
  data <- processed_sequences()$sRNAs

  ggplot(data, aes(x = Length)) +
    geom_histogram(binwidth = 10, fill = "#574c85", color = "black") +
    labs(
      title = "Distribution of sRNA Sequence Lengths",
      subtitle = paste("Across", length(unique(data$File)), "file(s)"),
      x = "Sequence Length (nt)",
      y = "Count"
    ) +
    facet_wrap(~ File, scales = "free_x") +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
})


output$group_plot <- renderPlot({
  req(processed_sequences()$sRNAs)
  data <- processed_sequences()$sRNAs

  ggplot(data, aes(x = Group, fill = Group)) +
    geom_bar(color = "black") +
    labs(
      title = "Distribution of sRNA Functional Groups",
      subtitle = paste("Across", length(unique(data$File)), "file(s)"),
      x = "Group Classification",
      y = "sRNA Count"
    ) +
    facet_wrap(~ File, scales = "free_x") +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
})

output$group_pie_chart <- renderPlot({
  req(processed_sequences()$sRNAs)
  df <- processed_sequences()$sRNAs
  df %>% count(Group) %>%
    ggplot(aes(x = "", y = n, fill = Group)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y") +
    theme_void() +
    labs(title = "sRNA Group Distribution")
})

output$strand_distribution_plot <- renderPlot({
  req(processed_sequences()$sRNAs)
  df <- processed_sequences()$sRNAs
  ggplot(df, aes(x = Strand, fill = Strand)) +
    geom_bar(color = "black") +
    labs(title = "Strand Distribution of sRNAs", x = "Strand", y = "Count") +
    theme_minimal()
})


output$gc_content_plot <- renderPlot({
  req(processed_sequences()$sRNAs)
  data <- processed_sequences()$sRNAs

  ggplot(data, aes(x = GC_content)) +
    geom_histogram(binwidth = 1, fill = "#8676d0", color = "black") +
    labs(
      title = "GC Content Distribution of sRNAs",
      subtitle = "Calculated per sequence and grouped by file",
      x = "GC Content (%)",
      y = "Number of sRNAs"
    ) +
    facet_wrap(~ File, scales = "free_x") +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
})


output$lengthclass <- renderPlot({
  req(processed_sequences()$sRNAs)
  data <- processed_sequences()$sRNAs

  ggplot(data, aes(x = Length, fill = Group)) +
    geom_histogram(binwidth = 10, color = "black") +
    labs(
      title = "Distribution of sRNA Sequence Lengths by Group",
      subtitle = "Each group shows binned lengths of sRNAs",
      x = "sRNA Length (nucleotides)",
      y = "Number of sRNAs"
    ) +
    facet_wrap(~ Group, scales = "free_y") +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
})



output$comparison <- renderPlot({
  req(processed_sequences()$sRNAs)
  data <- processed_sequences()$sRNAs

  ggplot(data, aes(x = Type, y = Length)) +
    geom_boxplot(fill = "#b1a2f0", color = "black", outlier.shape = 16, outlier.size = 2) +
    labs(
      title = "Comparison of sRNA Lengths by Experimental Condition",
      subtitle = "Boxplot showing the distribution of sRNA lengths across groups",
      x = "Experimental Condition (Group Label)",
      y = "sRNA Length (nucleotides)"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
})





output$lengthvsgc <- renderPlot({
  req(processed_sequences()$sRNAs)
  data <- processed_sequences()$sRNAs

  ggplot(data, aes(x = Length, y = GC_content, color = Type)) +
    geom_point(alpha = 0.6, size = 2) +
    labs(
      title = "Relationship Between sRNA Length and GC Content",
      subtitle = "Colored by experimental condition",
      x = "sRNA Length (nucleotides)",
      y = "GC Content (%)",
      color = "Condition"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "bottom"
    ) +
    guides(color = guide_legend(override.aes = list(size = 4)))
})


output$strandvsgc <- renderPlot({
  req(processed_sequences()$sRNAs)
  data <- processed_sequences()$sRNAs

  ggplot(data, aes(x = Strand, y = GC_content, fill = Strand)) +
    geom_violin(trim = FALSE, alpha = 0.7, color = "black") +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    labs(
      title = "GC Content Distribution by sRNA Strand",
      subtitle = "Violin and boxplots show distribution and median",
      x = "Strand (+ / -)",
      y = "GC Content (%)"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "none"
    )
})



output$groupcount <- renderPlot({
  req(processed_sequences()$sRNAs)
  data <- processed_sequences()$sRNAs

  ggplot(data, aes(x = Group)) +
    geom_bar(fill = "#B22222", color = "black", width = 0.7) +
    labs(
      title = "Distribution of sRNAs Across Different Groups",
      subtitle = "Counts of sRNAs categorized by group type",
      x = "Group Type",
      y = "Number of sRNAs"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
})

output$upset_group_plot <- renderPlot({
  req(processed_sequences()$sRNAs)
  sRNAs <- processed_sequences()$sRNAs

  sRNAs_grouped <- sRNAs %>%
    select(ID, Group, File) %>%
    distinct()

  id_by_group_file <- sRNAs_grouped %>%
    unite("GroupFile", Group, File, sep = "_") %>%
    mutate(present = 1) %>%
    pivot_wider(names_from = GroupFile, values_from = present, values_fill = list(present = 0))

  rownames(id_by_group_file) <- id_by_group_file$ID
  id_by_group_file <- id_by_group_file[, -1]

 
  UpSetR::upset(
  as.data.frame(id_by_group_file),
  sets = colnames(id_by_group_file),
  order.by = "freq",
  keep.order = TRUE,
  main.bar.color = "#4C72B0",
  sets.bar.color = "#55A868",
  number.angles = 45,
  point.size = 3,
  line.size = 1,
  text.scale = c(1.4, 1.4, 1.2, 1.2, 1.5, 1.2)
)

})







output$length_by_group <- renderPlot({
  req(processed_sequences()$sRNAs)

  data <- processed_sequences()$sRNAs

  ggplot(data, aes(x = Group, y = Length, fill = Group)) +
    geom_boxplot(outlier.color = "red", outlier.size = 1.5) +
    labs(
      title = "Distribution of sRNA Lengths Across Groups",
      subtitle = "Each box represents the distribution of sRNA lengths in a group",
      x = "sRNA Group",
      y = "Length (nucleotides)"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 12, color = "gray40"),
      legend.position = "none"
    )
})



output$free3_distribution <- renderPlot({
  req(structures_data())

  structures <- structures_data()

  
  structures$`free nt at the 3 end`  <- as.numeric(as.character(structures$`free nt at the 3 end`))

  ggplot(structures, aes(x = `free nt at the 3 end` )) +
    geom_histogram(binwidth = 1, fill = "#117711", color = "black") +
    labs(
      title = "Distribution of Free 3' End Nucleotides",
      subtitle = "Number of free (unpaired) nucleotides at the 3' end of sRNA structures",
      x = "Free 3' Nucleotides",
      y = "Count (Number of sRNAs)"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 12, color = "gray40"),
      axis.text.x = element_text(angle = 0),
      panel.grid.major.x = element_blank()
    )
})




 reactiveMyList <- reactive({
    
    list.files("./www/", pattern = "\\.jpg$", full.names = FALSE)
  })


observe({
  new_choices <- reactiveMyList()
 
  if (!is.null(input$gene) && input$gene %in% new_choices) {
    updateSelectInput(session, "gene", choices = new_choices, selected = input$gene)
  } else if (length(new_choices) > 0) {
    updateSelectInput(session, "gene", choices = new_choices, selected = new_choices[1])
  } else {
    updateSelectInput(session, "gene", choices = character(0), selected = NULL)
  }
})





  observeEvent(input$process_btn, {
    session$sendCustomMessage("show_spinner", list())

    req(input$fasta_file, input$file_ptt)

   custom_dir <- "www"
  if (!dir.exists(custom_dir)) {
    dir.create(custom_dir)
  }
 

  print("Starting to process sequences...")


    addResourcePath("customData", custom_dir)

    csv_data <- processed_sequences()$sRNAs
final_df <- csv_data
    fasta_data <- readDNAStringSet(input$fasta_file$datapath)

    ptt_data <- read.delim(input$file_ptt$datapath, skip = 2, header = TRUE, stringsAsFactors = FALSE)

write_multifasta <- function(df, filepath) {
  con <- file(filepath, "w")  
  on.exit(close(con))  
 
  for (i in seq_len(nrow(df))) {
    header <- paste0(">", df$ID[i])
    sequence <- df$Sequence[i]
   
   
    writeLines(header, con)
    writeLines(sequence, con)
  }
}


custom_fasta_folder <- "sRNA_FASTA_files"
if (!dir.exists(custom_fasta_folder)) {
  dir.create(custom_fasta_folder)
}





#naudoju
fasta_file_path <- "sequences.fasta"
write_multifasta(final_df, fasta_file_path)


sequences <- readDNAStringSet(fasta_file_path)


seq_ids <- names(sequences)


names(sequences) <- seq_ids


intermediate_dir <- "Generated_additional_files"
if (!dir.exists(intermediate_dir)) {
  dir.create(intermediate_dir)
}
file.copy("relplot.pl", "Generated_additional_files/relplot.pl", overwrite = TRUE)

Sys.chmod("Generated_additional_files/relplot.pl", mode = "0755")
process_sequence_1 <- function(sequence, seq_id) {
  fasta_filename <- file.path(custom_fasta_folder, paste0(seq_id, ".fasta"))
  writeXStringSet(sequence, fasta_filename)
 
  fasta_out <- file.path(intermediate_dir, paste0(seq_id, ".out"))
  ss_ps_filename <- file.path(intermediate_dir, paste0(seq_id, "_ss.ps"))
  dp_ps_filename <- file.path(intermediate_dir, paste0(seq_id, "_dp.ps"))
  rss_ps_filename <- file.path(intermediate_dir, paste0(seq_id, "_rss.ps"))
 
  command <- sprintf("cd %s && RNAfold -p < %s > %s", intermediate_dir, shQuote(file.path("..", fasta_filename)), shQuote(basename(fasta_out)))
system(command)

 
  relplot_command <- paste("Generated_additional_files/relplot.pl",
                         shQuote(ss_ps_filename),
                         shQuote(dp_ps_filename),
                         ">",
                         shQuote(rss_ps_filename))
system(relplot_command)


}

for (seq_id in seq_ids) {
  sequence <- sequences[seq_id]  
  process_sequence_1(sequence, seq_id)
}

process_sequence_2 <- function(seq_id) {
  input_ps <- file.path(intermediate_dir, paste0(seq_id, "_rss.ps"))
  output_jpg <- file.path(custom_dir, paste0(seq_id, ".jpg"))

  command <- sprintf("gs -sDEVICE=jpeg -dNOPAUSE -dQUIET -dBATCH -dSAFER -dFitPage -r300 -dJPEGQ=95 -sOutputFile='%s' '%s'",
                     output_jpg, input_ps)
  system(command)
}

for (seq_id in seq_ids) {
  process_sequence_2(seq_id)
}



expected_jpgs <- paste0(seq_ids, ".jpg")
output_dir <- custom_dir
timeout <- 15  
start_time <- Sys.time()

repeat {
  existing_jpgs <- list.files(output_dir, pattern = "\\.jpg$", full.names = FALSE)
  if (all(expected_jpgs %in% existing_jpgs)) {
    break
  }
  if (as.numeric(Sys.time() - start_time, units = "secs") > timeout) {
    warning("Timeout waiting for JPG files.")
    break
  }
  Sys.sleep(0.5)
}

updateSelectInput(session, "gene", choices = expected_jpgs, selected = expected_jpgs[1])



on.exit({
  data_ready(TRUE)
    session$sendCustomMessage("hide_spinner", list())
  }, add = TRUE)


    })

structures_data <- eventReactive(input$process_btn, {
parse_rnafold_output <- function(file) {
  lines <- readLines(file)
  if (length(lines) < 3) {
   
    return(data.frame(
      ID = NA,
      Sequence_5to3 = NA,
      Structure = NA,
      MFE = NA,
      free_3prime = NA,
      stringsAsFactors = FALSE
    ))
  }
 
 
  id_line <- lines[1]
  id_val <- sub("^>", "", id_line)  
 
 
  seq_val <- lines[2]
 
 
  mfe_line_idx <- grep("\\(\\-?\\d+\\.\\d+\\)$", lines)
 
  if (length(mfe_line_idx) == 0) {
   
    return(data.frame(
      ID = id_val,
      Sequence_5to3 = seq_val,
      Structure = NA,
      MFE = NA,
      free_3prime = NA,
      stringsAsFactors = FALSE
    ))
  }
 
 
  mfe_line <- lines[mfe_line_idx[1]]
 
 
  pattern <- "^(.*) \\((-?\\d+\\.\\d+)\\)$"
 
  structure_val <- sub(pattern, "\\1", mfe_line)
  mfe_val       <- sub(pattern, "\\2", mfe_line)
 
 
  count_trailing_dots <- function(x) {
    rev_chars <- rev(strsplit(x, "")[[1]])
    count <- 0
    for (ch in rev_chars) {
      if (ch == ".") count <- count + 1 else break
    }
    return(count)
  }
  free_3prime <- count_trailing_dots(structure_val)
 
 
  data.frame(
    ID = id_val,
    Sequence_5to3 = seq_val,
    Structure = structure_val,
    MFE = as.numeric(mfe_val),
    free_3prime = free_3prime,
    stringsAsFactors = FALSE
  )
}


out_files <- list.files(path = "./Generated_additional_files/", pattern = "\\.out$", full.names = TRUE)


parsed_list <- lapply(out_files, parse_rnafold_output)

df_structures <- bind_rows(parsed_list)




                 
sekos<-processed_sequences()$sRNAs

processed_filtered <- subset(processed_sequences()$sRNAs, ID %in% df_structures$ID)
  processed_filtered <- processed_filtered[!duplicated(processed_filtered$ID), ]
  df_structures <- df_structures[!duplicated(df_structures$ID), ]

  if (!"ID" %in% colnames(processed_filtered) || !"ID" %in% colnames(df_structures)) {
    showNotification("Missing ID column in structure merge", type = "error")
    return(NULL)
  }

df_final <- merge(
  x = processed_filtered,
  y = df_structures,
  by = "ID",
  all = FALSE  
)

df_final<- df_final[ ,c(1,4,5,7,8,14,15,16,17)]


colnames(df_final) <- c("ID", "Start", "End",  "Strand", "Group", "Sequence (5' → 3')", "Structure", "Minimum free energy", "free nt at the 3 end")




  df_final

})

   

    output$image <- renderImage({
      req(input$gene)
       
        filename <- normalizePath(file.path("./www/", input$gene))
       
       
        list(src = filename, style = "max-width: 200%; max-height: 1200px; height: auto;")
    }, deleteFile = FALSE)







output$structure_table <- renderDT({
  datatable(
    structures_data() %>%
      distinct(Start, End, .keep_all = TRUE),
    options = list(scrollX = TRUE)
  )
})



output$dataReady <- reactive({
  data_ready()
})
outputOptions(output, "dataReady", suspendWhenHidden = FALSE)


}

shinyApp(ui, server)