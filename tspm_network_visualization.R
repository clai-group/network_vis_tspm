if(!require(pacman)) install.packages("pacman")
require(dplyr)
require(tidyr)
#devtools::install_github("clai-group/mlho")
require(mlho)
require(DT)
require(igraph)
require(ggraph)
require(RColorBrewer)
pacman::p_load(data.table, devtools, backports, Hmisc, tidyr,dplyr,ggplot2,plyr,scales,readr,RcppParallel,
               httr, DT, lubridate, tidyverse,reshape2,foreach,doParallel,caret,gbm,lubridate,praznik,epitools,tcltk)

### Load dbmart -----------------
load() ### include at least three columns, patient_num, start_date and phenx (ICD codes here)
dbmart <- dplyr::select(dbmart,patient_num,start_date,phenx) %>% distinct()

### TSPM+ -----------------

db <- tSPMPlus::transformDbMartToNumeric(dbmart)
phenxlookup <- db$phenxLookUp
patlookup <- db$patientLookUp

phenxOfInterest = c(as.numeric(db$phenxLookUp$num_Phenx)) ### for the dementia/AD work we only are interested in those codes. 
## it should be other codes of interest and not general at this point
temporalBucket =  c(0,1,3)
minDuration = 0 #techical parameter, ignore for now ##TODO if not working with 0 use 1
bitShift = 0 #techical parameter, ignore for now
lengthOfPhenx = 7 #techical parameter
storeSequencesDuringCreation = FALSE  #if true, old way -> writing out "plain" sparse sequences in patient based files, FALSE-> do in memory sparsity

mem_buffer <- 1 #in GB. just a buffer to make sure the computer wont crash
cores_buffer <- 1
buffer<- 100000000*mem_buffer
dbmart_adapt <- tSPMPlus::splitdbMartInChunks(db$dbMart, includeCorSeq = TRUE, buffer = buffer)
numOfChunks = length(dbmart_adapt$chunks)
dbmart_num <- dbmart_adapt$chunks[[1]]

sparsity = 0.0001
numOfThreads = detectCores()-cores_buffer

corseq <- tSPMPlus::getCandidateSequencesForPOI(dbmart_num,
                                                minDuration,
                                                bitShift,
                                                lengthOfPhenx,
                                                temporalBucket,
                                                phenxOfInterest,
                                                storeSequencesDuringCreation,
                                                numOfThreads = numOfThreads,
                                                sparsityValue = sparsity)


corseq <- dplyr::distinct(corseq, .keep_all = TRUE)


corseq <- corseq %>%
  dplyr::group_by(patient_num,sequence,endPhenx,durationBucket) %>%
  dplyr::summarise(count=length((patient_num)))

corseq$value.var <- 1

corseq$startPhen <-substr(corseq$sequence, start = 1, stop = (nchar(corseq$sequence)-lengthOfPhenx))

corseq$startPhen_dur <- paste0(corseq$startPhen,"-",corseq$durationBucket)

dat <- data.table(corseq[c("patient_num","endPhenx","value.var","startPhen_dur")])
wide.dat.start <- dcast.data.table(dat, patient_num ~ startPhen_dur, value.var="value.var", fun=length)

end <- phenxOfInterest

corseq$value.var <- 1

corseq$startPhen <-substr(corseq$sequence, start = 1, stop = (nchar(corseq$sequence)-lengthOfPhenx))

corseq$startPhen_dur <- paste0(corseq$startPhen,"-",corseq$durationBucket)

dat <- data.table(corseq[c("patient_num","endPhenx","value.var","startPhen_dur")])
wide.dat.start <- dcast.data.table(dat, patient_num ~ startPhen_dur, value.var="value.var", fun=length)

end <- phenxOfInterest

for (i in seq(1:numOfChunks)) {
  gc()
  cores<-detectCores()
  cl <- parallel::makeCluster(cores-cores_buffer)
  doParallel::registerDoParallel(cl)
  tryCatch({
    corrs <- foreach(j = 1:length(end),#
                     .combine='rbind',
                     .multicombine=TRUE,
                     .packages = c("data.table")) %dopar% {
                       tryCatch({
                         dat.j <- data.table(corseq[corseq$endPhenx %in% end[j],c("patient_num","endPhenx","value.var")])
                         lab.j <- dcast.data.table(dat.j, patient_num ~ endPhenx, value.var="value.var", fun=length)
                         colnames(lab.j)[2] <- "label"
                         wide.j <- merge(wide.dat.start,lab.j,by="patient_num",all.x = T)
                         wide.j[is.na(wide.j)] <- 0
                         wide.j$patient_num <- NULL
                         
                         
                         out <- apply(wide.j[, -("label")], 2, cor.test, wide.j$label, method="spearman")
                         p <- data.frame(sapply(out, "[[", "p.value"))
                         p$startPhen_dur <- rownames(p)
                         rownames(p) <- NULL
                         colnames(p)[1] <- "p.value"
                         p$p.adjust <- p.adjust(p$p.value, method = "holm", n = nrow(p))# consider "holm" or "hochberg" or "bonferroni"
                         
                         rho <- data.frame(sapply(out, "[[", "estimate"))
                         rho$startPhen_dur <- rownames(rho)
                         rownames(rho) <- NULL
                         colnames(rho)[1] <- "rho"
                         rho$startPhen_dur <- substr(rho$startPhen_dur,1,nchar(rho$startPhen_dur)-4)
                         cor.j <- merge(rho,p,by="startPhen_dur")
                         cor.j$rho.abs <- abs(cor.j$rho)
                         cor.j$endPhenx <- end[j]
                         
                         
                         rm(rho,p);gc()
                         
                         
                         cor.j
                         
                       },
                       error = function(foll) {cat("ERROR :",conditionMessage(foll), "\n")})
                     }
    
    
    corrs <- merge(corrs,db$phenxLookUp,by.x = "endPhenx",by.y ="num_Phenx" )
    
    corrs$startPhen <- sub("\\-.*", "", corrs$startPhen_dur)
    corrs$sequence <- ifelse (corrs$endPhenx <10,paste0(corrs$startPhen,"000000",corrs$endPhenx),NA)
    corrs$sequence <- ifelse(corrs$endPhenx <100 & corrs$endPhenx >9 ,paste0(corrs$startPhen,"00000",corrs$endPhenx),corrs$sequence)
    corrs$sequence <- ifelse(corrs$endPhenx <1000 & corrs$endPhenx >99 ,paste0(corrs$startPhen,"0000",corrs$endPhenx),corrs$sequence)
    
    corrs$sequence <- as.numeric(corrs$sequence)

  },
  error = function(foll) {cat("ERROR in chunk ",i, ": ", conditionMessage(foll), "\n")})
  stopCluster(cl)
  
}

corrs_final <- corrs %>%
  dplyr::select(-phenx) %>%
  dplyr::group_by(endPhenx,startPhen_dur,sequence,startPhen) %>%
  dplyr::summarise_all(mean, na.rm=TRUE)


### Initial Vizualisations
## 1.) TOP N Sequence in unique patients

# Get frequency of sequences from corseq df
corseq_freq <- corseq %>%
  group_by(sequence) %>%
  summarise(frequency = n())
corseq_freq$sequence <- as.character(corseq_freq$sequence)

# Plot a simple line plot to visualize available sequences' frequency 
# Data frames with a high number of unique sequences might not visualize well
ggplot(corseq_freq, aes(x = seq_along(frequency), y = frequency)) +
  # Add a smooth line
  geom_line(color = "blue") +
  # Add labels to axes
  labs(x = "Sequence", y = "Frequency") +
  # Add a theme for better appearance
  theme_minimal()

### COMPARISON OF TWO DATASETS PLOT ###
# Generate random integers for an artifical frequency. Comment out if comparison is not wanted
# or add code to load in a second data frame containing corss_final/corseq_freq
corseq_freq$artificial_freq <- abs(round(rnorm(nrow(corseq_freq), mean = 20, sd = 62))) #distribution and parameters are randomly chosen

# For two line plots
# Requires corseq_freq to have a new column called "artificial_freq" containing values from the data frame to compare from 
ggplot(corseq_freq, aes(x = seq_along(sequence))) + 
  geom_line(aes(y = frequency, color = "Observed Frequency")) +
  geom_line(aes(y = artificial_freq, color = "Artificial Frequency")) +
  labs(title = "Line plots for Frequency of Sequences", y = "Frequency", x = "Sequence") +
  theme_minimal() +
  scale_color_manual(name = "Frequency Type", values = c("Observed Frequency" = "blue", "Artificial Frequency" = "red"))
### END COMPARISON OF TWO DATASETS PLOT ###

## 2.) TOP N Sequence for correlation

data_graph <- corrs_final[corrs_final$p.adjust<=0.05, ] %>% 
  dplyr::select (startPhen, endPhenx, rho) %>% 
  dplyr::group_by(startPhen, endPhenx) %>%
  dplyr::summarize(mean_val = mean(rho, na.rm = TRUE), freq=n()) %>% 
  distinct()

data_graph <- merge(data_graph, phenxlookup, by.x ="startPhen", by.y = "num_Phenx")
colnames(data_graph)[colnames(data_graph) == 'phenx'] <- 'start'
data_graph <- merge(data_graph, phenxlookup, by.x ="endPhenx", by.y = "num_Phenx")
colnames(data_graph)[colnames(data_graph) == 'phenx'] <- 'end'

# Smooth line plot
# Order data frame first
data_graph_rho <- data_graph[order(-data_graph$mean_val), ]
ggplot(data_graph_rho, aes(x = seq_along(mean_val), y = mean_val)) +
  # Add a smooth line
  geom_line(color = "blue") +
  # Add labels to axes
  labs(x = "Sequence", y = "Correlation Coefficient") +
  # Add a theme for better appearance
  ggtitle("TOP N Sequence (Correlation)") +
  theme_minimal()

### COMPARISON OF TWO DATASETS PLOT ###
# For two line plots
# data_graph_rho must contain a column "artificial_mean_val" that contains the mean values of a second data frame
data_graph_rho$artificial_mean_val <- runif(nrow(data_graph), min = -0.2, max = 1)

# For two line plots
# data_graph_rho must contain a column "artificial_mean_val" that contains the mean values of a second data frame
ggplot(data_graph_rho, aes(x = seq_along(mean_val))) + 
  geom_line(aes(y = mean_val, color = "Observed Correlation Coefficient")) +
  geom_line(aes(y = artificial_mean_val, color = "Artificial Correlation Coefficient")) +
  labs(title = "Line plots for Correlation Coefficients") +
  theme_minimal() +
  scale_color_manual(name = "Data Set", values = c("Observed Correlation Coefficient" = "blue", "Artificial Correlation Coefficient" = "red"))
### END COMPARISON OF TWO DATASETS PLOT ###

## 3.) TOP N in frequency of ICD Codes

# Create a data frame connect that includes the start and end phenX of each sequence
connect <- data_graph %>% dplyr::select (start, end, mean_val, freq)
colnames(connect) <- c('from','to','value', 'freq') 

# Create nodes data frame containing the frequency of each phenX
c(as.character(connect$from), as.character(connect$to)) %>%
  as.tibble() %>%
  group_by(value) %>%
  dplyr::summarize(n=n()) -> nodes

# Create a barplot from data frame nodes that depicts the frequency of each phenX (ICD Code)
# Currently does not label Y axis due to too many phenX being present
nodes <- nodes[order(nodes$n, decreasing = TRUE), ]
barplot(nodes$n,
        xlab = "ICD Codes", ylab = "Frequency",
        main = "Barplot of ICD Codes",
        las = 2,
        col = "skyblue",
        names.arg = rep("", nrow(nodes)))

### COMPARISON OF TWO DATASETS PLOT ###
# Create an artifical variable with the same length of icd codes. Make sure to have the same length when comparing two data sets.
nodes$n2 <- sample(1:65, size = nrow(nodes), replace = TRUE)
nodes_long <- pivot_longer(nodes, cols = c(n2, n), names_to = "group", values_to = "frequency")

ggplot(nodes_long, aes(x = value, y = frequency, fill = group)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "Bar plots for ICDs of two data sets", x = "ICD Code", y = "Frequency") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
### END COMPARISON OF TWO DATASETS PLOT ###

### Dendrogram: requires a df with at least two columns phenotype and CONCEPT_CD -----------
## Q: edges color and edges width (frequency of each connection)
## ref: https://r-graph-gallery.com/310-custom-hierarchical-edge-bundling.html
corseq <- as.data.frame(corseq)
corseq$patient_num <- as.character(corseq$patient_num)
corseq$sequence <- as.character(corseq$sequence)
freq <- dplyr::select(corseq, startPhen, endPhenx, patient_num) %>%
  distinct() %>%
  dplyr::group_by(startPhen, endPhenx) %>%
  dplyr::summarise(freq = n_distinct(patient_num))

data_graph <- merge(data_graph, freq, by.x = c("startPhen", "endPhenx"))

codes <- unique(c(unique(data_graph$start), unique(data_graph$end)))

phen_label <- df %>% dplyr::select (phenotype, CONCEPT_CD) %>% 
  distinct() %>% 
  group_by(CONCEPT_CD) %>%
  dplyr::mutate(phenotype = ifelse(n() > 1, paste(phenotype, collapse = "&"), phenotype)) %>%
  distinct()
unique(phen_label$phenotype) # 4

phen_label <- subset(phen_label, phen_label$CONCEPT_CD %in% codes)

d1 <- data.frame(from="origin", to=unique(phen_label$phenotype))
d2 <- phen_label
colnames(d2)<- c("from", "to")
hierarchy <- rbind(d1, d2)

phe <- phen_label  %>% 
  group_by(phenotype) %>%  dplyr::summarize(n=n())

colnames(nodes) <- c("name", "value")
colnames(phe) <- c("name", "value")
orign <- data.frame(name ="origin", value="68")

vertices <- rbind(phe, nodes, orign)
vertices$group  <-  hierarchy$from[match( vertices$name, hierarchy$to ) ]
vertices$value = as.numeric(vertices$value)
vertices <- data.frame(vertices)

mygraph <- graph_from_data_frame(hierarchy, vertices=vertices)

from  <-  match(connect$from, vertices$name)
to  <-  match(connect$to, vertices$name)


p <- ggraph(mygraph, layout = 'dendrogram', circular = TRUE) +
  geom_conn_bundle(data = get_con(from = from, to = to), aes(colour= connect$value[after_stat(index)], 
                                                             width = connect$freq[after_stat(index)] )) +
  scale_edge_color_continuous(low="white", high="red")
#scale_edge_colour_distiller(palette = "RdPu")
p

p+ 
  geom_node_point(aes(filter = leaf, x = x*1.05, y=y*1.05, colour=group, size=vertices$value)) +
  scale_colour_manual(values= rep( brewer.pal(9,"Paired") , 30)) +
  scale_size_continuous(range = c(0.1,10) )  +
  geom_node_text(aes(label = name),  colour = 'black', size=3, nudge_x = p$data$x * .15, nudge_y = p$data$y * .15) +
  theme_void()









### Network --------------
# Q: same phenotype didn't group together, legend, edge color
# ref: https://r.igraph.org/articles/igraph.html#structural-properties-of-graphs

data_graph <- data_graph %>% dplyr::select (start, end, mean_val, freq)

c( as.character(data_graph$start), as.character(data_graph$end)) %>%
  as.tibble() %>%
  group_by(value) %>%
  dplyr::summarize(n=n()) -> nodes
colnames(nodes) <- c("name", "n")
nodes <- merge(nodes, phen_label, by.x="name", by.y="CONCEPT_CD")

graph <- graph_from_data_frame(data_graph, directed = FALSE, vertices=nodes)

edge_values <- E(graph)$mean_val
edge_values_norm <- (edge_values - min(edge_values)) / (max(edge_values) - min(edge_values)) 
edge_color_gradient <- colorRamp(c("white", "red")) 
edge_colors <- edge_color_gradient(edge_values_norm)


coul  <- brewer.pal(4, "Set1") 
#vertice_colors <- coul[as.numeric(as.factor(V(graph)$phenotype))]

#vertex_order <- order(V(graph)$phenotype) 
plot(graph,
     layout = layout_in_circle(graph),
     #layout = layout_with_fr(graph, weights = abs(E(graph)$mean_val)),
     #layout = layout_in_circle(graph)[vertex_order, ], 
     vertex.size = sqrt(nodes$n),
     # vertex.color = vertice_colors,
     vertex.label.color = "black",     
     vertex.label.cex = 0.8, 
     edge.width = sqrt(data_graph$freq),
     edge.color = edge_colors
)
# legend("bottomright", legend=levels(as.factor(V(graph)$phenotype)),
#        col = coul , 
#        bty = "n", 
#        pch=20 , pt.cex = 2 , cex = 1, 
#        text.col=coul , horiz = FALSE, inset = c(-0.35, 0.1))

