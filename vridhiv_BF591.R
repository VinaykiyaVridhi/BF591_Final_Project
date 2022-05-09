library(shiny)
library(ggplot2)
library(colourpicker)
library(rsconnect)
library(xtable)
library(tidyverse)
library('RColorBrewer')
library(DT)
library(shinyWidgets)
library(ggbeeswarm)
options(shiny.maxRequestSize=40*1024^2) 

js <- '.nav-tabs-custom .nav-tabs li.active {border-top-color:#ac2d1d;}"'

ui <- fluidPage(
    titlePanel("BF591 Final Project"),
    mainPanel(
        tabsetPanel(
            
            
            #TAB1
            tabPanel("Sample Data",
                     h5("This tab displays the metadata characterstics"),
                     sidebarPanel(fileInput("sample_data", "Load data", accept = ".csv"),
                                  radioButtons("hist", "Select parameter to make a graph", 
                                               c("mRNA_seq_reads", "PMI", "RIN","age_death")),
                                  colourInput("histcolor", "Histogram color"),
                                  actionButton("go", "Go")
                     ),
                     mainPanel(
                         tabsetPanel(
                             tabPanel("Summary",
                                      p("Summary of the data uploaded"),
                                      tableOutput("summary")
                             ),
                             tabPanel("Table",
                                      p("Table showing the data uploaded"),
                                      DT::dataTableOutput("table")
                             ),
                             tabPanel("Graph",
                                      p("Histogram showing the data uploaded"),
                                      plotOutput("histogram")
                             )
                         )
                     ),
            ), 
            
            
            
            #TAB2
            tabPanel("Counts",
                     h5("This tab displays the characterstics of normalised counts of huntington's disease experiments"),
                     sidebarPanel(fileInput("counts_data","Load counts data",accept = ".csv"),
                                  sliderInput("per_v", min = 0.1, max = 3.5, 
                                              label = "Choose X to Include genes with a percentile of variance",value = 1.2, step = 0.2),
                                  sliderInput("per_z", min = 5, max = 60, 
                                              label = "Choose X to Include genes with samples that are non-zero",value = 30, step = 5),
                                  selectInput("xaxis", "Select PC for x axis", 
                                              choices = c('PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10'), selected = 1),
                                  selectInput("yaxis", "Select PC for y axis", 
                                              choices = c('PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10'), selected = 3),
                                  actionButton("go", "Go")
                     ),
                     mainPanel(
                         tabsetPanel(
                             tabPanel("Data Table",
                                      p("This is a summary of the uploaded data"),
                                      tableOutput("filter_counts")
                             ),
                             tabPanel("Scatter Plot",
                                      h5("Diagnostic scatter plots"),
                                      p("A graph for median count vs variance"),
                                      plotOutput("plot1")
                             ),
                             tabPanel("PCA PLot",
                                      p("Scatter plot of PCA projections"),
                                      plotOutput("PCA_plot")
                             ),
                             tabPanel("Heat Map",
                                      p("Heatmap of the filtered normalised counts"),
                                      plotOutput("heat_plot")
                             )
                         )
                     ),
            ),
            
            
            
            #TAB3
            tabPanel("Differential Expression",
                     h5("This tab displays the characterstics of differential expression values of Huntington's disease experiment"),
                     sidebarPanel(fileInput("DE_data", "Load differential expression data", accept = ".csv"),
                                  radioButtons("x_para", "Choose the parameter for the x-axis", 
                                               c("baseMean","HD.mean","Control.mean","log2FoldChange","lfcSE","stat","pvalue","padj"), selected = "log2FoldChange"),
                                  radioButtons("y_para", "Choose the parameter for the y-axis", 
                                               c("baseMean","HD.mean","Control.mean","log2FoldChange","lfcSE","stat","pvalue","padj"), selected = "padj"),
                                  colourInput("color1", "Base point color"),
                                  colourInput("color2", "Highlight point color"),
                                  sliderInput("slider", min = -25, max = 0,
                                              label = "Select the slider of the p adjusted",value = -15, step = 5),
                                  actionButton("go", "Go")
                     ),
                     mainPanel(
                         tabsetPanel(
                             tabPanel("Table",
                                      DT::dataTableOutput('DE_table')
                             ),
                             
                             tabPanel("Volcano Plot",
                                      p('This is a volcano plot'),
                                      plotOutput('vol_plot')
                             ), 
                         ),
                     ), 
            ),
            
            
            #TAB 4
            tabPanel("Gene Expression",
                     h5("This tab displays the characterstics gene expression of huntington's disease experiment"),
                     sidebarPanel(fileInput("sample_data2", "Load sample data", accept = ".csv"),
                                  fileInput("counts_data2", "Load normalized counts data", accept = ".csv"),
                                  selectInput("categorical_value", label = ("Choose category to plot the data"), 
                                              choices = c('Control vs HD_Samples' = 1, 'Age at Death' = 2)),
                                  searchInput("input_id",
                                              label = NULL,
                                              value = "ENSG00000000003.10",
                                              placeholder = "ENSG00000000003.10",
                                              btnSearch = NULL,
                                              btnReset = NULL,
                                              resetValue = "",
                                              width = NULL
                                  ),
                                  radioButtons("plot_type","Choose a Plot type to view results", c('Barplot'= 1,'BeeSwarm' = 2, 'Boxplot'= 3,'Violin plot'= 4)),
                                  actionButton("go", "Go")
                     ),
                     mainPanel(
                         tabsetPanel(
                             tabPanel("Samples",
                                      DT::dataTableOutput('sample4')
                             ),
                             tabPanel("Normalized Counts Table",
                                      DT::dataTableOutput("table4")
                             ),
                             tabPanel("Plot",
                                      plotOutput("plot4")
                             ),
                         ),
                     ), 
            ),
        ) 
    )
)




server <- function(input, output, session)
{
    
    #SAMPLES TAB
    load_data <- reactive({read.csv(input$sample_data$datapath , row.names = 1)})
    
    sample_data_summary<- function(dataf)
    {
        categories <- colnames(dataf)
        #To find the datatype of data in each column
        type <- c(typeof(dataf$sample_id),typeof(dataf$PMI),typeof(dataf$age_death),
                  typeof(dataf$RIN),typeof(dataf$mRNA_seq_reads)) 
        
        #To calculate the mean of each numerical column, round to three digits
        mean <- c("unique",round(mean(dataf$PMI %>% na.omit),3),round(mean(dataf$age_death),3),
                  round(mean(dataf$RIN),3),round(mean(dataf$mRNA_seq_reads),3))
        
        #To calculate the standard deviation
        std_dev <- ifelse(sapply(dataf, is.numeric) == TRUE,round(apply(dataf,2,sd,na.rm=TRUE),3),"unique")
        summary <- data.frame(categories, type, mean, std_dev)
        return(summary)
    }
    
    #Graphical Representation of sample data
    histogram_plot <- function(df,category,bins)
    {
        df$category<- as.numeric(df[[category]])
        p <- ggplot(df, aes(category)) +
            geom_histogram(fill = "blue",bins = bins)
        return(p)
    }
    
    #To create summary table
    output$summary <- renderTable({req(input$sample_data)
        df<- load_data() 
        p<- sample_data_summary(df)
        return(p)
    }, height = 350)
    
    #Sortable data table
    output$table <- DT::renderDataTable({req(input$sample_data)
        DT::datatable(load_data())})
    
    output$histogram <- renderPlot({req(input$sample_data)
        df <- load_data()
        p <- hist(as.numeric(df[[input$hist]]), xlab= paste(input$hist), 
                  main = paste("Histogram showing ", input$hist), col = input$histcolor,
                  breaks=10) 
        return(p)},
        height = 350)
    
    
    #COUNTS TAB
    load_data_counts <- reactive({read.csv(input$counts_data$datapath,row.names = 1)})
    
    filter_counts <- function(counts,per_v,per_z)
    {
        df_var <- mutate(counts, variance = (apply(counts,1,sd)/apply(counts,1,mean)) > per_v)
        df_notzero <- mutate(df_var, zeros = rowSums(df_var == 0) < per_z) 
        return(df_notzero)
    }
    
    
    summary_counts <- function (counts,filtered)
    {
        filtered <- filtered %>% 
            filter(variance ==TRUE) %>% 
            select(-variance) %>% 
            filter(zeros ==TRUE) %>% 
            select(-zeros)
        options(scipen = 50)
        samples <- ncol(counts)
        genes <- nrow(counts)
        excluded_genes <- genes - nrow(filtered)
        new_genes_list <- nrow(filtered)
        df<- data.frame("Number of samples at the start"= c(samples,''),"Number of genes at the start"= c(genes,''),
                        "Genes excluded while filtering"= c(excluded_genes, excluded_genes/genes),
                        "Genes filtered"= c(new_genes_list, new_genes_list/genes))
        rownames(df) <- c('Numbers of Samples',"Percentage of Total Genes")
        return(df %>%  mutate_if(is.numeric, round, digits = 3))
    }
    
    #To create a Scatter plot
    scatter_cv<- function(filtered)
    {
        df <- mutate(filtered, median = apply(filtered,1,median)) %>%  
            mutate(variance_cv = (apply(filtered,1,sd)/apply(filtered,1,mean))) %>% 
            mutate(Filtered = ifelse(variance == FALSE | zeros ==FALSE ,'Filtered Out','Not Filtered Out')) 
        p <- ggplot(df,aes(median,variance_cv)) + geom_point(aes(color = factor(Filtered)),alpha = 0.5) + 
            scale_x_continuous(trans='log2')+
            scale_color_manual(values = c('Filtered Out' = "blue", "Not Filtered" = "grey"))+
            labs(title="Median counts vs Variance", x ="Log scale median counts", y = "Variance")
        return(p)
    } 
    
    #To create a Heatmap
    heatmap_plot<- function(f)
    {
        #Creating heatmap and filterng out NA
        f <- f %>% 
            filter(variance ==TRUE) %>% 
            select(-variance) %>% 
            filter(zeros ==TRUE)%>% select(-zeros)  
        f[f==0] <- NA
        f<- f  %>%
            log10() %>% 
            as.matrix() %>%  
            na.omit() %>% heatmap()
        return(p)
    }
    
    #To create PCA plot
    pca_plot<- function(counts,x,y)
    {
        #Transposing the data
        pca_results <- prcomp(scale(t(counts)), center=FALSE, scale=FALSE)
        matrix <- do.call(rbind.data.frame, pca_results)
        #Mutating rows to columns
        matrix_new <- mutate (matrix, col1 = rownames(matrix)) 
        Std_dev <- (pca_results$sdev) 
        Var <- Std_dev^2 
        total <- sum(Var) 
        var_div <- Var/total  
        new_tibble <- tibble(var_div = var_div)
        pc <- mutate(new_tibble, principal_components = seq.int(nrow(new_tibble)))
        pc$principal_components <- sub("^","PC", pc$principal_components) 
        pc_cum <- mutate(pc,cumulative = cumsum(pc$var_div))
        x_plot <- matrix_new[,x]
        y_plot <- matrix_new[,y]
        #Changing variance into percentage
        x_var <- round(pc_cum$var_div[strtoi(str_sub(x, start= -1))] *100,3) 
        y_var <- round(pc_cum$var_div[strtoi(str_sub(y, start= -1))] *100,3) 
        p <- ggplot()+geom_point(aes(x_plot,y_plot)) + 
            labs(title = paste(x,"vs", y), x = paste(x, ' (Variance Explained ',x_var,"%)"),
                 y = paste0(y,' (Variance Explainied ',y_var,"%)"))
        return(p)
    }
    
    output$filter_counts <- renderTable({req(input$counts_data)
        df <- load_data_counts()
        filter <- filter_counts(df,input$per_v,input$per_z)
        p <-summary_counts(df,filter)
        return(p)},rownames = TRUE, height = 450)
    output$plot1 <- renderPlot({req(input$counts_data)
        df <- load_data_counts()
        filter <- filter_counts(df,input$per_v,input$per_z)
        p <- scatter_cv(filter) 
        return(p)},height = 450)
    output$heat_plot <- renderPlot({req(input$counts_data)
        df <- load_data_counts()
        filter <- filter_counts(df,input$per_v,input$per_z)
        p<- heatmap_plot(filter)
        return(p)}, height = 450)
    output$PCA_plot <- renderPlot({req(input$counts_data)
        df <- load_data_counts()
        p <- pca_plot(df,input$x,input$y)
        return(p)}, height = 450)
    
    
    #DESeq TAB
    load_data_difexp <- reactive({read.csv(input$DE_data$datapath,row.names = 1)})
    
    draw_table <- function(dataf, slider) 
    {
        dataf <- subset(dataf,padj<1*10^slider) 
        dataf <- mutate(dataf, pvalue = format(pvalue, scientific = TRUE)) %>% mutate(padj = format(padj, scientific = TRUE))
        return(dataf)
    }
    
    #Volcano plot
    vol_plot <- function(dataf, xname, yname, slider, color1, color2) {
      p <- ggplot(data = dataf, 
                  aes(x =!!sym(xname), y=-log10(!!sym(yname)))) + 
        geom_point(aes(color = padj< 1*10^(slider))) +
        theme(legend.position = "bottom") +
        scale_color_manual(values = c('TRUE' = color1, 'FALSE' = color2)) +
        labs(x= xname, y=str_glue("-log10({yname})"), color=str_glue("{yname} < 10^{slider}"))
      return(p)
    }
    
    output$DE_table <- DT::renderDataTable({req(input$DE_data)
        DT::datatable(load_data_difexp())})
    

    output$vol_plot <- renderPlot({vol_plot(dataf = load_data_difexp(),
                                                 slider = input$slider,
                                                 xname = input$x_para,
                                                 yname = input$y_para,
                                                 color1 = input$color1,
                                                 color2 = input$color2)})
    
    
    #GeneExpression TAB
    counts_cyodata <- reactive({read.csv(input$counts_data2$datapath,row.names = 1)})
    
    sample_cyodata <- reactive({read.csv(input$sample_data2$datapath,row.names = 1)})
    
    gene<- reactive({input$input_id})
    
    age_div<- function(sample,norm,gene_id){
        dataf<-subset(norm,gene==gene_id)
        a_1 <- filter(sample,age_death>=36) %>% filter(age_death<55)
        a1_col <- a_1$sample_id
        a_2 <-filter(sample,age_death>=55) %>% filter(age_death<66)
        a2_col <- a_2$sample_id
        a_3 <-filter(sample,age_death>=66) %>% filter(age_death<75)
        a3_col <- a_3$sample_id
        a_4 <-filter(sample,age_death>=75) 
        a4_col <- a_4$sample_id
        a1 <- dataf %>% select(all_of(a1_col)) %>% unlist() %>% list() %>% lapply(function(x) {x[x!=0]}) 
        a2 <- dataf %>% select(all_of(a2_col)) %>% unlist() %>% list() %>% lapply(function(x) {x[x!=0]})
        a3 <- dataf %>% select(all_of(a3_col)) %>% unlist() %>% list() %>% lapply(function(x) {x[x!=0]})
        a4 <- dataf %>% select(all_of(a4_col)) %>% unlist() %>% list() %>% lapply(function(x) {x[x!=0]})
        t_1 <- data.frame(a1,letter = "Ages 36-55")
        colnames(t_1) = c("counts","group")
        t_2 <- data.frame(a2,letter = "Ages 55-66")
        colnames(t_2) = c("counts","group")
        t_3 <- data.frame(a3,letter = "Ages 66-75")
        colnames(t_3) = c("counts","group")
        t_4 <- data.frame(a4,letter = "Ages >75")
        colnames(t_4) = c("counts","group")
        ages <- rbind(t_1,t_2,t_3,t_4)
        return(ages)
    }
    
    control_vs_sample<- function(sample,norm,gene_id){
        dataf<-subset(norm,gene==gene_id)
        ch <-mutate(sample,letter=substr(sample_id,1,1)) 
        c <- filter(ch,letter=='C') 
        c_col <- c$sample_id
        h<- filter(ch,letter =="H")
        h_col <- h$sample_id
        c_cnts <- dataf %>% select(all_of(c_col)) %>% unlist() %>% list() %>% lapply(function(x) {x[x!=0]})
        h_cnts <- dataf %>% select(all_of(h_col)) %>% unlist() %>% list() %>% lapply(function(x) {x[x!=0]})
        t_h<-data.frame(h_cnts,letter = "H")
        t_c<- data.frame(c_cnts,letter = "C")
        colnames(t_h) = c("counts","group")
        colnames(t_c) = c("counts","group")
        control_vs_sample_df <- rbind(t_h,t_c)
        return(control_vs_sample_df)
    }
    
    bar_plot <-function(dataframe){
        p<- ggplot(dataframe, aes(counts, fill = group)) + 
            geom_histogram(alpha = 0.5, bins = 20, aes(y = ..density..), position = 'identity')
        return(p)
    }
    
    bee_plot<-function(dataframe){
        p<-ggplot(dataframe, aes(x=group, y=counts, color = group)) + geom_beeswarm() +
            xlab("Controls vs Affected Samples")
        return(p)
    }
    
    box_plot <- function(dataframe){
        p <- ggplot(dataframe, aes(x=group, y=counts, fill=group)) + 
            geom_boxplot(alpha=0.3) +
            theme(legend.position="none")
        return(p)
    }
    violin_plot <- function(dataframe){
        p <- ggplot(dataframe, aes(x=group, y=counts, fill=group)) + geom_violin() +
            xlab("Controls vs Affected Samples")
        return(p)
    }   
    
    output$table4 <- DT::renderDataTable({req(input$counts_data2)
        DT::datatable(counts_cyodata())})
    
    output$sample4 <- DT::renderDataTable({req(input$sample_data2)
        DT::datatable(sample_cyodata())})
    
    output$count3<- renderTable({req(input$counts_data2)
        df_counts<- load_counts() %>% rownames_to_column("gene")
        df_samples<- load_sample()
        gene_name <- input$input_id
        p<- age_div(df_samples,df_counts,gene_name)
        return(p)})
    
    output$plot4 <- renderPlot({req(input$counts_data2)
        df_counts<- counts_cyodata() %>% rownames_to_column("gene")
        df_samples<- sample_cyodata()
        gene_name <- gene()
        if(input$categorical_value ==1){
            data<- control_vs_sample(df_samples,df_counts,gene_name)
        }
        else{
            data<- age_div(df_samples,df_counts,gene_name)
        }
        if(input$plot_type==1){
            p<-bar_plot(data)
        }
        else if (input$plot_type==2){
            p<-bee_plot(data)
        }
        else if (input$plot_type==3){
            p<- box_plot(data)
        }
        else{
            p<- violin_plot(data)
        }
        return(p)})
    
    
}

shinyApp(ui = ui, server = server)