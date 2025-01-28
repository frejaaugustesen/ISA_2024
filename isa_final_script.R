# description -------------------------------------------------------------

# ISA project for fall semester 2024, 10 ECTS
# comparison of data quality from data provided by snATAC-seq methods from 10X Genomics
# using their monomodal snATAC-seq method (mono) and their multimodal snATAC-seq + RNA-seq method (multi) 
# using data from 20 different donors, each providing one sample per method, equaling
# 40 paired samples, to detemine differences in quality from the two methods.

# the samples were obtained from ND, pre-T2D and T2D patients, pancreas - more specific
# from the islet of langerhans in the pancreas.

# see source article of data by wang et al. 2022 - "Integrating genetics with single-cell multiomic 
# measurements across disease states identifies mechanisms of beta cell dysfunction in type 2 diabetes"
# doi: https://doi.org/10.1038/s41588-023-01397-9

# setup ------------------------------------------------------------------
.libPaths(here::here("library"))
# usethis::use_package("ArchR", min_version = TRUE)

# package load
library(ArchR)
library(parallel)
library(rstatix)
library(qs)
library(ggrastr)
library(ggpubr)
library(coin)

addArchRThreads(getArchRThreads()/2) # half of what is avaliable on the computer

addArchRGenome("hg38") # humane reference genome used for alignment of fragment files

set.seed(100) # to make sure of the same outcome every time

# load (ArrowFiles) --------------------------------------------------------------------

## inputFiles ----

#vector with paths to all fragment files to be used for inputFiles in createArrowFiles function
raw_data <- list.files(path = "/work/isa/data-raw/fragment_files", # found using getwd()
                       pattern = NULL, 
                       all.files = FALSE, # only visible files are written
                       full.names = TRUE, # gives path to all files
                       recursive = FALSE, 
                       ignore.case = FALSE, 
                       include.dirs = FALSE, 
                       no.. = FALSE)

## sampleNames ----

# using the stringr package to create a vector of names with only the sample ID's
raw_data_names <- str_extract(raw_data, "[^/]+(?=_frags.sort.bed.gz)")

## Arrowfile ----
raw_ArrowFile <- createArrowFiles(
  inputFiles = raw_data,
  sampleNames = raw_data_names,
  minTSS = 4, # not to be set to high - can be increased later
  minFrags = 1000,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)


# doublet scores ---------------------------------------------------------------
doubScores <- addDoubletScores(
  input = raw_ArrowFile,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
  LSIMethod = 1,
  threads = 1
)

# ArchRProject -----------------------------------------------------------------

## create project ----
projISA_v0 <- ArchRProject(
  ArrowFiles = raw_ArrowFile,
  outputDirectory = "ISA_project",
  copyArrows = TRUE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)

## save project ----
saveArchRProject(ArchRProj = projISA_v0, outputDirectory = "Save-ISA_project", load = FALSE)

## load project ----
projISA_v0 <- loadArchRProject(path = here::here("Save-ISA_project/"))

# metafilt color ---------------------------------------------------------------

# Load necessary library
library(colorspace) # for destaturation
library(pals) # for color pallete
library(dplyr) # dataframe wrangling
library(readxl) # excel file processing

# Load meta data
meta <- read_excel ("/work/isa/data-raw/meta.xlsx")

# Generate 20 distinct HEX colors for the donors
set.seed(100)  # For reproducibility
donor_id_colors <- pals::glasbey(n = 20)  # 20 distinct colors
names(donor_id_colors) <- unique(meta$donor_id) # add donor_ids to each color

# Create a data frame to runique()# Create a data frame to represent the samples, donors and method
meta_filt <- meta %>%
  dplyr::mutate(method = case_when(multiome == "Yes" ~ "multi", # when multiome == "Yes", then write multi
                                   .default = "mono")) %>%  # .dafult when the multiome == "Yes" returns FALSE or NA, then write mono
  dplyr::select(id, donor_id, method) # keep only these columns

# Add colors to meta_filt dataframe based on the donor_ids
meta_filt$color <- donor_id_colors[match(meta_filt$donor_id, names(donor_id_colors))]

# ligten the color if the method is "multi" otherwise keep the original color
meta_filt <- meta_filt %>% dplyr::mutate(color = case_when(method == "multi" ~ colorspace::lighten(color, 0.5),
                                                           .default = color))
meta_filt # look at it again
View(meta_filt)

# extract colors - this vector can later be used to color each sample
sample_colors <- meta_filt %>%
  dplyr::pull(color, id)

# plot the colors
meta_filt %>%
  ggplot(aes(x = donor_id, y = rep(1, 40), fill = id)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = sample_colors) +
  theme(legend.position = "none")

# plot colors and split by method
meta_filt %>%
  ggplot(aes(x = donor_id, y = rep(1, 40), fill = id)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = sample_colors) +
  facet_wrap(~ method) +
  theme(legend.position = "none")

# simple color -----------------------------------------------------------------
# Define sample names for each treatment group
mono_samples <- c("MM_108", "MM_59", "MM_93", "MM_60", "MM_95", "MM_55", "MM_61", "MM_87",
                  "MM_88", "MM_89", "MM_77", "MM_78", "MM_109", "MM_51", "MM_96", "MM_97",
                  "MM_122", "MM_98", "MM_110", "MM_123")

multi_samples <- c("MM_747", "MM_746", "MM_743", "MM_724", "MM_745", "MM_734", "MM_744",
                   "MM_722", "MM_723", "MM_733", "MM_742", "MM_748", "MM_717", "MM_735",
                   "MM_718", "MM_719", "MM_736", "MM_503", "MM_741", "MM_721")

# Create a named vector of colors
sample_colors_simple <- c(
  setNames(rep("mediumpurple1", length(mono_samples)), mono_samples),
  setNames(rep("darkorange", length(multi_samples)), multi_samples))

# subsetting ArchR project ----------------------------------------------
# load object
projISA_test <- loadArchRProject(path = "/work/isa/Save-ISA_project_v1/")
# save test
saveArchRProject(ArchRProj = projISA_test, outputDirectory = "ISA_project_test", load = TRUE)

# dividing archr project into two based on mono and multi:
# Make a new isa_dataframe_merged dataframe
isa_dataframe <- as.data.frame(getCellColData(projISA_test))
isa_dataframe$row <- rownames(isa_dataframe)
isa_dataframe_merged_2 <- merge(isa_dataframe, meta_filt, by.x = "Sample", by.y = "id", all.x = TRUE)
rownames(isa_dataframe_merged_2) <- isa_dataframe_merged_2$row

# Barcodes in mono
cells_mono <- isa_dataframe_merged_2 %>% dplyr::filter(method == "mono") %>% rownames()
# Subset archr
archr_mono <- subsetArchRProject(
  ArchRProj = projISA_test,
  cells = cells_mono,
  outputDirectory = "archr_mono",
  dropCells = TRUE,
  force = TRUE,
  threads = 16)
# Add method
archr_mono <- addCellColData(
  archr_mono,
  data = rep("mono", length(archr_mono$cellNames)),
  cells = archr_mono$cellNames,
  name = "method")
# Barcodes in multi
cells_multi <- isa_dataframe_merged_2 %>% dplyr::filter(method == "multi") %>% rownames()
# Subset archr
archr_multi <- subsetArchRProject(
  ArchRProj = projISA_test,
  cells = cells_multi,
  outputDirectory = "archr_multi",
  dropCells = TRUE,
  force = TRUE,
  threads = 16)
# Add method
archr_multi <- addCellColData(
  archr_multi,
  data = rep("multi", length(archr_multi$cellNames)),
  cells = archr_multi$cellNames,
  name = "method")

# TSS enrichment profiles -------------------------------------------------
TSSenrichment_profile <- plotTSSEnrichment(ArchRProj = projISA_v0,
                                           groupBy = "Sample")
# changing colors on profile
TSSenrichment_profile_simplecolor <- TSSenrichment_profile +
  ggplot2::scale_colour_manual(values = sample_colors_simple) +
  theme(legend.position = "none")

TSSenrichment_profile_simplecolor +
  labs(title = "TSS enrichment profile for all samples") +
  annotation_custom(
    grob = grid::rectGrob(gp = grid::gpar(fill = "mediumpurple1", col = "black")),
    xmin = 1400, xmax = 1500, ymin = 11.25, ymax = 11.75
  ) +
  annotation_custom(
    grob = grid::rectGrob(gp = grid::gpar(fill = "darkorange", col = "black")),
    xmin = 1400, xmax = 1500, ymin = 10.5, ymax = 11
  ) +
  annotate("text", x = 1600, y = 11.5, label = "mono", hjust = 0) +
  annotate("text", x = 1600, y = 10.75, label = "multi", hjust = 0) +
  annotate("text", x = 1600, y = 12.5, label = "Method", hjust = 0.5, fontface = "bold")

## individual profiles ----
tss_sep_sample_simple <- TSSenrichment_profile +
  ggplot2::facet_wrap(~group,
                      scales = "free") +
  ggplot2::theme(legend.position = "none") +
  ggplot2::scale_colour_manual(values = sample_colors_simple)

## profile per method -----------------------------------------------------------
### creating profile from subset ArchR projects
TSSenrichment_profile_mono <- plotTSSEnrichment(ArchRProj = archr_mono,
                                                groupBy = "Sample")

TSSenrichment_profile_multi <- plotTSSEnrichment(ArchRProj = archr_multi,
                                                 groupBy = "Sample")

### plot
TSSenrichment_profile_mono_simple <- TSSenrichment_profile_mono +
  ggplot2::scale_colour_manual(values = sample_colors_simple) +
  theme(legend.position = "none") +
  labs(title = "TSS enrichment profile for mono")

TSSenrichment_profile_multi_simple <- TSSenrichment_profile_multi +
  ggplot2::scale_colour_manual(values = sample_colors_simple) +
  theme(legend.position = "none") +
  labs(title = "TSS enrichment profile for multi")

# fragment size distribution ----------------------------------------------
fragment_size <- ArchR::plotFragmentSizes(ArchRProj = projISA_v0,
                                          threads = 1,
                                          groupBy = "Sample")
# plotting with colors
fragment_size_simplecolor <- fragment_size +
  ggplot2::scale_colour_manual(values = sample_colors_simple,
                               labels = c("Mono", "Multi")) +
  theme(legend.position = "none")

fragment_size_simplecolor +
  labs(title = "Fragment size distribution for all samples") +
  annotation_custom(
    grob = grid::rectGrob(gp = grid::gpar(fill = "darkorange", col = "black")),
    xmin = 635, xmax = 655, ymin = 1.07, ymax = 1.12
  ) +
  annotation_custom(
    grob = grid::rectGrob(gp = grid::gpar(fill = "mediumpurple", col = "black")),
    xmin = 635, xmax = 655, ymin = 1.17, ymax = 1.22
  ) +
  annotate("text", x = 660, y = 1.2, label = "mono", hjust = 0) +
  annotate("text", x = 660, y = 1.1, label = "multi", hjust = 0) +
  annotate("text", x = 660, y = 1.3, label = "Method", hjust = 0.5, fontface = "bold")

## individual distribution plots ----
fragment_size_simple <- fragment_size +
  ggplot2::facet_wrap(~group,
                      scales = "free") +
  ggplot2::theme(legend.position = "none") +
  ggplot2::scale_colour_manual(values = sample_colors_simple)

## distribution plot per method ----
fragment_size_mono <- ArchR::plotFragmentSizes(ArchRProj = archr_mono,
                                               threads = 1,
                                               groupBy = "Sample")

fragment_size_multi <- ArchR::plotFragmentSizes(ArchRProj = archr_multi,
                                                threads = 1,
                                                groupBy = "Sample")

## plot
fragment_size_mono_simple <- fragment_size_mono +
  ggplot2::scale_colour_manual(values = sample_colors_simple) +
  theme(legend.position = "none") +
  labs(title = "Fragment size distribution for mono") +
  ggplot2::ylim(0,1.2)

fragment_size_multi_simple <- fragment_size_multi +
  ggplot2::scale_colour_manual(values = sample_colors_simple) +
  theme(legend.position = "none") +
  labs(title = "Fragment size distribution for multi") +
  ggplot2::ylim(0,1.2)

# filtering cells ---------------------------------------------------------

#load project
projISA_filter <- loadArchRProject(path = "/work/isa/Save-ISA_project_filter/")
# save test
saveArchRProject(ArchRProj = projISA_filter, outputDirectory = "ISA_project_filter", load = TRUE)

## filtering by threshold ----
projISA_filter <- projISA_filter[projISA_filter$nFrags >= 2500 &
                                   projISA_filter$TSSEnrichment >= 7 &
                                   projISA_filter$BlacklistRatio <= 0.05]

ArchR::getCellColData(projISA_filter)

saveArchRProject(ArchRProj = projISA_filter, outputDirectory = "Save-ISA_project_filter", load = FALSE)

projISA_filter <- loadArchRProject(path = "/work/isa/Save-ISA_project_filter/")

# reviewing data and comparing
isa_dataframe_filtered <- as.data.frame(getCellColData(projISA_filter))

# creating merged dataframe with all data and multi/mono and ID
isa_dataframe_merged_filtered <- merge(isa_dataframe_filtered, meta_filt, by.x = "Sample", by.y = "id", all.x = TRUE)

# counting cells
count_before <- isa_dataframe_merged %>%
  dplyr::count(method)
count_before

count_after <- isa_dataframe_merged_filtered %>%
  dplyr::count(method)
count_after

## cell count difference ----
percentage_difference <- count_before %>%
  dplyr::rename(before = n) %>%
  dplyr::inner_join(count_after %>% dplyr::rename(after = n), by = "method") %>%
  mutate(percentage_diff = ((after - before) / before) * 100)

# cell count per sample before and after filtering:

# creating dataframes before and after
data_count_filtered <- isa_dataframe_merged_filtered %>%
  dplyr::count(donor_id, Sample, method)

data_count_unfiltered <- isa_dataframe_merged %>%
  dplyr::count(donor_id, Sample, method)

# renaming columns
data_count_filtered <- data_count_filtered %>%
  dplyr::rename(filtered = n)

data_count_unfiltered <- data_count_unfiltered %>%
  dplyr::rename(unfiltered = n)

## making complete dataframe 

data_count_unfiltered$filtered <- data_count_filtered$filtered

cell_count_all <- data_count_unfiltered

cell_count_all$difference <- c(((cell_count_all$filtered - cell_count_all$unfiltered)/cell_count_all$unfiltered)*100)

## plot ----

# Creating the bar plot - by method
cell_count_all %>%
  mutate(Sample = factor(Sample, levels = unique(Sample[order(method, donor_id)]))) %>%
  ggplot(aes(x = donor_id, y = abs(difference), fill = method)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~method, scales = "free_x") +
  scale_fill_manual(values = c("mono" = "mediumpurple1", "multi" = "darkorange")) +
  labs(
    x = "Donor",
    y = "'Bad' cells (%)",
    fill = "Method",
    title = "Percentage of poor quality nuclei per sample"
  ) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray92"),
        panel.grid.minor = element_line(color = "gray92"),
        panel.border = element_rect(color = "gray92", fill = NA) )

# Creating the bar plot - by donor
cell_count_all %>%
  mutate(Sample = factor(Sample, levels = unique(Sample[order(donor_id, method)]))) %>%
  ggplot(aes(x = donor_id, y = abs(difference), fill = method)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("mono" = "mediumpurple1", "multi" = "darkorange")) +
  labs(
    x = "Donor",
    y = "Percentage of poor quality nuclei",
    fill = "Method",
    title = "Percentage of poor quality nuclei per sample"
  ) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray92"),
        panel.grid.minor = element_line(color = "gray92"),
        panel.border = element_rect(color = "gray92", fill = NA) )

# Creating the bar plot - per method
cell_count_all %>%
  mutate(method = factor(method, levels = c("mono", "multi"))) %>%
  ggplot(aes(x = method, y = abs(difference), fill = method)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("mono" = "mediumpurple1", "multi" = "darkorange")) +
  labs(
    x = "Donor",
    y = "'Bad' cells (%)",
    fill = "Method",
    title = "Percentage of poor quality nuclei per method"
  ) +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray92"),
        panel.grid.minor = element_line(color = "gray92"),
        panel.border = element_rect(color = "gray92", fill = NA) ) +
  annotate("text", x = 1.5, y = 70, label = "T-test, p-value = 0.0297", size = 4, color = "black")

# Creating the bar plot - per method (mean)
cell_count_all %>%
  mutate(difference_abs = abs(difference)) %>% # Convert difference to absolute values
  group_by(method) %>% # Group by method
  summarize(mean_difference = mean(difference_abs, na.rm = TRUE)) %>% # Calculate mean for each method
  ggplot(aes(x = method, y = mean_difference, fill = method)) +
  geom_bar(stat = "identity", position = "dodge") + # Create a bar plot
  scale_fill_manual(values = c("mono" = "mediumpurple1", "multi" = "darkorange")) +
  labs(
    x = "Method",
    y = "percentage of poor quality nuceli (mean)",
    fill = "Method",
    title = "Percentage of poor quality nuclei per method"
  ) +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray92"),
        panel.grid.minor = element_line(color = "gray92"),
        panel.border = element_rect(color = "gray92", fill = NA)) +
  annotate("text", x = 1.5, y = 40, label = "T-test, p-value = 0.0297", size = 4, color = "black")

# peak calling -----------------------------------------------------------------------
# using subsetted ArchRProject
## Add reproducible peakset
# Find peaks
pathToMacs2 <- "/work/isa/conda/myenv/bin/macs2"
# Generate pseudobulk replicates for peak calling
# We set it up so that all cells in each replicates are aggregated together (pseudobulk)
# so that we can do tradition peak calling with macs2
archr_mono <- addGroupCoverages(
  ArchRProj = archr_mono,
  groupBy = "method",
  minCells = min(table(archr_mono@cellColData$Sample)), # minimum number of cells for each replicate
  maxCells = max(table(archr_mono@cellColData$Sample)), # maximum number of cells for each replicate
  minReplicates = 2, # minimum number of replicates
  maxReplicates = 20, # maximum number of replciates, 20 because we have 20 replicates.
  force = TRUE,
  threads = 16
)
archr_multi <- addGroupCoverages(
  ArchRProj = archr_multi,
  groupBy = "method",
  minCells = min(table(archr_multi@cellColData$Sample)), # minimum number of cells for each replicate
  maxCells = max(table(archr_multi@cellColData$Sample)), # maximum number of cells for each replicate
  minReplicates = 2, # minimum number of replicates
  maxReplicates = 20, # maximum number of replciates, 20 because we have 20 replicates.
  force = TRUE,
  threads = 16
)
# Find peaks
archr_mono <- addReproduciblePeakSet(
  ArchRProj = archr_mono,
  groupBy = "method",
  pathToMacs2 = pathToMacs2,
  threads = 16)
archr_multi <- addReproduciblePeakSet(
  ArchRProj = archr_multi,
  groupBy = "method",
  pathToMacs2 = pathToMacs2,
  threads = 16)
# Add peakmatrix
# Compute counts for each peak per cell in the provided ArchRProject
archr_mono <- ArchR::addPeakMatrix(
  ArchRProj = archr_mono,
  threads = 16,
  logFile = createLogFile("addPeakMatrix")
)
archr_multi <- ArchR::addPeakMatrix(
  ArchRProj = archr_multi,
  threads = 16,
  logFile = createLogFile("addPeakMatrix")
)
saveArchRProject(ArchRProj = archr_mono, outputDirectory = "archr_mono", load = TRUE)
saveArchRProject(ArchRProj = archr_multi, outputDirectory = "archr_multi", load = TRUE)

# paired t test, preliminary tests --------------------------------------------------

##  TSS enrichment score ----

#summary statistics
isa_dataframe_merged %>%
  get_summary_stats(TSSEnrichment,type = "mean_sd")

#summary statistics by method
isa_dataframe_merged %>%
  group_by(method) %>%
  get_summary_stats(TSSEnrichment,type = "mean_sd")


# assumptions and preliminary tests

TSS_long <- isa_dataframe_merged %>%
  dplyr::select(donor_id, method, TSSEnrichment) %>% 
  group_by(donor_id, method) %>% 
  dplyr::summarise(mean_TSS = mean(TSSEnrichment)) %>%
  dplyr::ungroup()


TSS_wide <- TSS_long %>%
  tidyr::pivot_wider(id_cols = donor_id,
                     names_from = method,
                     values_from = mean_TSS) 

#visualization
ggpaired(TSS_long, x = "method", y = "mean_TSS",
         order = c("mono", "multi"),
         ylab = "TSSEnrichment", xlab = "method")

# Summary statistics

TSS_long %>%
  group_by(method) %>%
  get_summary_stats(mean_TSS, type = "mean_sd")

# Visualization
box_comparison_TSS <- ggpaired(TSS_long, x = "method", y = "mean_TSS",
                               order = c("mono", "multi"),
                               ylab = "mean_TSS", xlab = "method")

# comparison
TSS_wide$difference <- c(TSS_wide$mono - TSS_wide$multi)

# identifying outliers
TSS_wide %>%
  identify_outliers(difference)

# Shapiro-Wilk normality test for the differences
shapiro_test(TSS_wide$difference)

# historgram
TSS_wide %>%
  ggplot(aes(x = difference)) +
  geom_density()

# QQ plot for the difference
ggqqplot(TSS_wide, "difference")

## fragment count ----

#summary statistics
isa_dataframe_merged %>%
  get_summary_stats(nFrags,type = "mean_sd")

#summar statistics by method
isa_dataframe_merged %>%
  group_by(method) %>%
  get_summary_stats(nFrags,type = "mean_sd")

# assumptions and preliminary tests

nFrags_long <- isa_dataframe_merged %>%
  dplyr::select(donor_id, method, nFrags) %>% # vælger kolonner jeg skal arbejde med
  group_by(donor_id, method) %>% # grupperer efter metode og id (for hver donor per metode)
  dplyr::summarise(mean_nFrags = mean(nFrags)) %>%
  dplyr::ungroup()

nFrags_wide <- nFrags_long %>%
  tidyr::pivot_wider(id_cols = donor_id,
                     names_from = method,
                     values_from = mean_nFrags) # opdeler dataframe i flere koloner så det er 20 rækker og mono+multi for hver

#visualization
ggpaired(nFrags_long, x = "method", y = "mean_nFrags",
         order = c("mono", "multi"),
         ylab = "nFrags", xlab = "method")

# Summary statistics

nFrags_long %>%
  group_by(method) %>%
  get_summary_stats(mean_nFrags, type = "mean_sd")

# Visualization
box_comparison_nFrags <- ggpaired(nFrags_long, x = "method", y = "mean_nFrags",
                                  order = c("mono", "multi"),
                                  ylab = "mean_nFrags", xlab = "method")
box_comparison_nFrags

# comparison
nFrags_wide$difference <- c(nFrags_wide$mono - nFrags_wide$multi)

# identifying outliers
nFrags_wide %>%
  identify_outliers(difference)

# Shapiro-Wilk normality test for the differences
shapiro_test(nFrags_wide$difference)

# laver historgram
nFrags_wide %>%
  ggplot(aes(x = difference)) +
  geom_density()

# QQ plot for the difference
ggqqplot(nFrags_wide, "difference")

## normalized fragment count ----

#summary statistics
isa_dataframe_merged %>%
  get_summary_stats(normFragCount,type = "mean_sd")

#summar statistics by method
isa_dataframe_merged %>%
  group_by(method) %>%
  get_summary_stats(normFragCount,type = "mean_sd")

# assumptions and preliminary tests

normFragCount_long <- isa_dataframe_merged %>%
  dplyr::select(donor_id, method, normFragCount) %>% 
  group_by(donor_id, method) %>% 
  dplyr::summarise(mean_normFragCount = mean(normFragCount)) %>%
  dplyr::ungroup()

normFragCount_wide <- normFragCount_long %>%
  tidyr::pivot_wider(id_cols = donor_id,
                     names_from = method,
                     values_from = mean_normFragCount) 

#visualization
ggpaired(normFragCount_long, x = "method", y = "mean_normFragCount",
         order = c("mono", "multi"),
         ylab = "normFragCount", xlab = "method")

# Summary statistics

normFragCount_long %>%
  group_by(method) %>%
  get_summary_stats(mean_normFragCount, type = "mean_sd")

# Visualization
box_comparison_normFragCount <- ggpaired(normFragCount_long, x = "method", y = "mean_normFragCount",
                                         order = c("mono", "multi"),
                                         ylab = "mean_normFragCount", xlab = "method")

# comparison
normFragCount_wide$difference <- c(normFragCount_wide$mono - normFragCount_wide$multi)

# identifying outliers
normFragCount_wide %>%
  identify_outliers(difference)

# Shapiro-Wilk normality test for the differences
shapiro_test(normFragCount_wide$difference)

# historgram
normFragCount_wide %>%
  ggplot(aes(x = difference)) +
  geom_density()

# QQ plot for the difference
ggqqplot(normFragCount_wide, "difference")


## blacklist ratio ----

#summary statistics
isa_dataframe_merged %>%
  get_summary_stats(BlacklistRatio,type = "mean_sd")

#summar statistics by method
isa_dataframe_merged %>%
  group_by(method) %>%
  get_summary_stats(BlacklistRatio,type = "mean_sd")

# assumptions and preliminary tests

BlacklistRatio_long <- isa_dataframe_merged %>%
  dplyr::select(donor_id, method, BlacklistRatio) %>% 
  group_by(donor_id, method) %>% 
  dplyr::summarise(mean_BlacklistRatio = mean(BlacklistRatio)) %>%
  dplyr::ungroup()

BlacklistRatio_wide <- BlacklistRatio_long %>%
  tidyr::pivot_wider(id_cols = donor_id,
                     names_from = method,
                     values_from = mean_BlacklistRatio) 

#visualization
ggpaired(BlacklistRatio_long, x = "method", y = "mean_BlacklistRatio",
         order = c("mono", "multi"),
         ylab = "BlacklistRatio", xlab = "method")

# Summary statistics

BlacklistRatio_long %>%
  group_by(method) %>%
  get_summary_stats(mean_BlacklistRatio, type = "mean_sd")

# Visualization
box_comparison_BlacklistRatio <- ggpaired(BlacklistRatio_long, x = "method", y = "mean_BlacklistRatio",
                                          order = c("mono", "multi"),
                                          ylab = "mean_BlacklistRatio", xlab = "method")
box_comparison_BlacklistRatio

# comparison
BlacklistRatio_wide$difference <- c(BlacklistRatio_wide$mono - BlacklistRatio_wide$multi)

# identifying outliers
BlacklistRatio_wide %>%
  identify_outliers(difference)

# Shapiro-Wilk normality test for the differences
shapiro_test(BlacklistRatio_wide$difference)

# historgram
BlacklistRatio_wide %>%
  ggplot(aes(x = difference)) +
  geom_density()

# QQ plot for the difference
ggqqplot(BlacklistRatio_wide, "difference")

## Doublet score ----

#summary statistics
isa_dataframe_merged %>%
  get_summary_stats(DoubletScore,type = "mean_sd")

#summar statistics by method
isa_dataframe_merged %>%
  group_by(method) %>%
  get_summary_stats(DoubletScore,type = "mean_sd")

# assumptions and preliminary tests

DoubletScore_long <- isa_dataframe_merged %>%
  dplyr::select(donor_id, method, DoubletScore) %>% 
  group_by(donor_id, method) %>% 
  dplyr::summarise(mean_DoubletScore = mean(DoubletScore)) %>%
  dplyr::ungroup()

DoubletScore_wide <- DoubletScore_long %>%
  tidyr::pivot_wider(id_cols = donor_id,
                     names_from = method,
                     values_from = mean_DoubletScore) 

#visualization
ggpaired(DoubletScore_long, x = "method", y = "mean_DoubletScore",
         order = c("mono", "multi"),
         ylab = "DoubletScore", xlab = "method")

# Summary statistics

DoubletScore_long %>%
  group_by(method) %>%
  get_summary_stats(mean_DoubletScore, type = "mean_sd")

# Visualization
box_comparison_DoubletScore <- ggpaired(DoubletScore_long, x = "method", y = "mean_DoubletScore",
                                        order = c("mono", "multi"),
                                        ylab = "mean_DoubletScore", xlab = "method")

# comparison
DoubletScore_wide$difference <- c(DoubletScore_wide$mono - DoubletScore_wide$multi)

# identifying outliers
DoubletScore_wide %>%
  identify_outliers(difference)

# Shapiro-Wilk normality test for the differences
shapiro_test(DoubletScore_wide$difference)

# historgram
DoubletScore_wide %>%
  ggplot(aes(x = difference)) +
  geom_density()

# QQ plot for the difference
ggqqplot(DoubletScore_wide, "difference")

##  FRiP ----

#summary statistics
archr_combined_df_2 %>%
  get_summary_stats(FRIP_percent,type = "mean_sd")

#summar statistics by method
archr_combined_df_2 %>%
  group_by(method) %>%
  get_summary_stats(FRIP_percent,type = "mean_sd")

# assumptions and preliminary tests

FRIP_long <- archr_combined_df_2 %>%
  dplyr::select(donor_id, method, FRIP_percent) %>% 
  group_by(donor_id, method) %>% 
  dplyr::summarise(mean_FRIP = mean(FRIP_percent)) %>% 
  dplyr::ungroup()

FRIP_wide <- FRIP_long %>%
  tidyr::pivot_wider(id_cols = donor_id,
                     names_from = method,
                     values_from = mean_FRIP) 

#visualization
ggpaired(FRIP_long, x = "method", y = "mean_FRIP",
         order = c("mono", "multi"),
         ylab = "FRIP", xlab = "method")

# Summary statistics
FRIP_long %>%
  group_by(method) %>%
  get_summary_stats(mean_FRIP, type = "mean_sd")

# Visualization
box_comparison_FRIP <- ggpaired(FRIP_long, x = "method", y = "mean_FRIP",
                                order = c("mono", "multi"),
                                ylab = "mean_FRIP", xlab = "method")

# comparison
FRIP_wide$difference <- c(FRIP_wide$mono - FRIP_wide$multi)

# identifying outliers
FRIP_wide %>%
  identify_outliers(difference)

# Shapiro-Wilk normality test for the differences
shapiro_test(FRIP_wide$difference)

# historgram
FRIP_wide %>%
  ggplot(aes(x = difference)) +
  geom_density()

# QQ plot for the difference
ggqqplot(FRIP_wide, "difference")

## Poor Quality nuclei  ----

#summary statistics
cell_count_all %>%
  get_summary_stats(difference,type = "mean_sd")

#summar statistics by method
cell_count_all %>%
  group_by(method) %>%
  get_summary_stats(difference,type = "mean_sd")


# assumptions and preliminary tests
Quality_long <- cell_count_all %>%
  dplyr::select(donor_id, method, difference) %>% 
  group_by(donor_id, method) %>% 
  dplyr::summarise(mean_Quality = mean(difference)) %>% 
  dplyr::ungroup()


Quality_wide <- Quality_long %>%
  tidyr::pivot_wider(id_cols = donor_id,
                     names_from = method,
                     values_from = mean_Quality) 

#visualization
ggpaired(Quality_long, x = "method", y = "mean_Quality",
         order = c("mono", "multi"),
         ylab = "Difference", xlab = "method")

# Summary statistics
Quality_long %>%
  group_by(method) %>%
  get_summary_stats(mean_Quality, type = "mean_sd")

# Visualization
box_comparison_Quality <- ggpaired(Quality_long, x = "method", y = "mean_Quality",
                                   order = c("mono", "multi"),
                                   ylab = "mean_Quality", xlab = "method")

# comparison
Quality_wide$difference <- c(Quality_wide$mono - Quality_wide$multi)

# identifying outliers
Quality_wide %>%
  identify_outliers(difference)

# Shapiro-Wilk normality test for the differences
shapiro_test(Quality_wide$difference)

# historgram
Quality_wide %>%
  ggplot(aes(x = difference)) +
  geom_density()

# QQ plot for the difference
ggqqplot(Quality_wide, "difference")

# paired t-test --------------------------------------------------------------

## fragment count -----

# significant difference 
stat.test_nFrags <- nFrags_long  %>%
  t_test(mean_nFrags ~ method, paired = TRUE) %>%
  add_significance()

# effect size
nFrags_long  %>%
  dplyr::ungroup() %>%
  cohens_d(mean_nFrags ~ method, paired = TRUE)

# report
report_stat.test_nFrags <- stat.test_nFrags %>% add_xy_position(x = "method")

box_comparison_nFrags +
  stat_pvalue_manual(report_stat.test_nFrags, tip.length = 0) +
  labs(subtitle = get_test_label(report_stat.test_nFrags, detailed= TRUE))

## normalized fragment count -----

# significant difference 
stat.test_normFragCount <- normFragCount_long  %>%
  t_test(mean_normFragCount ~ method, paired = TRUE) %>%
  add_significance()

# effect size
normFragCount_long  %>%
  dplyr::ungroup() %>%
  cohens_d(mean_normFragCount ~ method, paired = TRUE)

# report
report_stat.test_normFragCount <- stat.test_normFragCount %>% add_xy_position(x = "method")

box_comparison_normFragCount +
  stat_pvalue_manual(report_stat.test_normFragCount, tip.length = 0) +
  labs(subtitle = get_test_label(report_stat.test_normFragCount, detailed= TRUE))


## blacklist ratio ----

# significant difference
stat.test_BlacklistRatio <- BlacklistRatio_long  %>%
  dplyr::ungroup() %>%
  t_test(mean_BlacklistRatio ~ method, paired = TRUE) %>%
  add_significance()

# effect size 
BlacklistRatio_long  %>%
  cohens_d(mean_BlacklistRatio ~ method, paired = TRUE)

# report
report_stat.test_BlacklistRatio <- stat.test_BlacklistRatio %>% add_xy_position(x = "method")

box_comparison_BlacklistRatio +
  stat_pvalue_manual(report_stat.test_BlacklistRatio, tip.length = 0) +
  labs(subtitle = get_test_label(report_stat.test_BlacklistRatio, detailed= TRUE))

## FRiP -----

# significant difference 
stat.test_FRIP <- FRIP_long  %>%
  t_test(mean_FRIP ~ method, paired = TRUE) %>%
  add_significance()

# effect size
FRIP_long  %>%
  dplyr::ungroup() %>%
  cohens_d(mean_FRIP ~ method, paired = TRUE)

# report
report_stat.test_FRIP <- stat.test_FRIP %>% add_xy_position(x = "method")

box_comparison_FRIP +
  stat_pvalue_manual(report_stat.test_FRIP, tip.length = 0) +
  labs(subtitle = get_test_label(report_stat.test_FRIP, detailed= TRUE))

## poor quality nuclei -----

# significant difference - kan ikke få denne til at virke
stat.test_Quality <- Quality_long  %>%
  t_test(mean_Quality ~ method, paired = TRUE) %>%
  add_significance()

# effect size
Quality_long  %>%
  dplyr::ungroup() %>%
  cohens_d(mean_Quality ~ method, paired = TRUE)

# report
report_stat.test_Quality <- stat.test_Quality %>% add_xy_position(x = "method")

box_comparison_Quality +
  stat_pvalue_manual(report_stat.test_Quality, tip.length = 0) +
  labs(subtitle = get_test_label(report_stat.test_Quality, detailed= TRUE))

# wilcoxon test, preliminary tests --------------------------------------------------------

## TSS enrichment ----

# assumptions and preliminary tests
TSSEnrichment_long <- isa_dataframe_merged %>%
  dplyr::select(donor_id, method, TSSEnrichment) %>% 
  group_by(donor_id, method) %>% 
  dplyr::summarise(mean_TSSEnrichment = mean(TSSEnrichment)) %>%  
  dplyr::ungroup()

TSSEnrichment_wide <- TSSEnrichment_long %>%
  tidyr::pivot_wider(id_cols = donor_id,
                     names_from = method,
                     values_from = mean_TSSEnrichment) 

# summary statistics
TSSEnrichment_long %>%
  group_by(method) %>%
  get_summary_stats(mean_TSSEnrichment, type = "median_iqr")

# visualization
box_comparison_TSS <- ggpaired(TSSEnrichment_long, x = "method", y = "mean_TSSEnrichment",
                               order = c("mono", "multi"),
                               ylab = "mean_TSSEnrichment", xlab = "method")

# assumptions and preliminary tests
TSSEnrichment_wide$differences <- c(TSSEnrichment_wide$mono - TSSEnrichment_wide$multi)

wilcox_assumption_TSSEnrichment <- gghistogram(TSSEnrichment_wide, x = "differences", y = "..density..",
                                               fill = "steelblue",bins = 5, add_density = TRUE)

## doublet Score ----

# assumptions and preliminary tests
DoubletScore_long <- isa_dataframe_merged %>%
  dplyr::select(donor_id, method, DoubletScore) %>% 
  group_by(donor_id, method) %>% 
  dplyr::summarise(mean_DoubletScore = mean(DoubletScore)) %>%  
  dplyr::ungroup()

DoubletScore_wide <- DoubletScore_long %>%
  tidyr::pivot_wider(id_cols = donor_id,
                     names_from = method,
                     values_from = mean_DoubletScore) 

# summary statistics
DoubletScore_long %>%
  group_by(method) %>%
  get_summary_stats(mean_DoubletScore, type = "median_iqr")

# visualization
box_comparison_DoubletScore <- ggpaired(DoubletScore_long, x = "method", y = "mean_DoubletScore",
                                        order = c("mono", "multi"),
                                        ylab = "mean_DoubletScore", xlab = "method")

# assumptions and preliminary tests
DoubletScore_wide$differences <- c(DoubletScore_wide$mono - DoubletScore_wide$multi)

wilcox_assumption_DoubletScore <- gghistogram(DoubletScore_wide, x = "differences", y = "..density..",
                                              fill = "steelblue",bins = 5, add_density = TRUE)

# wilcoxon test ----------------------------------------------------------------

## TSS enrichment ----
wilcox_stat.test_TSSEnrichment <- TSSEnrichment_long  %>%
  dplyr::ungroup() %>%
  rstatix::wilcox_test(mean_TSSEnrichment ~ method, paired = TRUE) %>%
  add_significance()

wilcox.test(TSSEnrichment_wide$mono, TSSEnrichment_wide$multi,
            paired = TRUE)

# effect size
TSSEnrichment_long  %>%
  dplyr::ungroup() %>%
  rstatix::wilcox_effsize(mean_TSSEnrichment ~ method, paired = TRUE)

# report
wilcox_report_TSSEnrichment <- wilcox_stat.test_TSSEnrichment %>% add_xy_position(x = "method")

box_comparison_TSS +
  stat_pvalue_manual(wilcox_report_TSSEnrichment, tip.length = 0) +
  labs(subtitle = get_test_label(wilcox_report_TSSEnrichment, detailed= TRUE))

## doublet Score ----
wilcox_stat.test_DoubletScore <- DoubletScore_long  %>%
  dplyr::ungroup() %>%
  rstatix::wilcox_test(mean_DoubletScore ~ method, paired = TRUE) %>%
  add_significance()

# effect size
DoubletScore_long  %>%
  wilcox_effsize(mean_DoubletScore ~ method, paired = TRUE)

# report
wilcox_report_DoubletScore <- wilcox_stat.test_DoubletScore %>% add_xy_position(x = "method")
box_comparison_DoubletScore +
  stat_pvalue_manual(wilcox_report_DoubletScore, tip.length = 0) +
  labs(subtitle = get_test_label(wilcox_report_DoubletScore, detailed= TRUE))

# plotting results ------------------------------------------------------------

## normalized fragment count ----

# 1/3 norm fragment count - per sample, sorted by method
isa_dataframe_merged %>%
  mutate(Sample = factor(Sample, levels = unique(Sample[order(method, donor_id)]))) %>%
  ggplot(aes(x = donor_id, y = normFragCount, fill = method)) +
  geom_boxplot(outlier.size = 0.2) +
  facet_wrap(~method, scales = "free_x") +
  scale_fill_manual(values = c("mono" = "mediumpurple1", "multi" = "darkorange")) +
  expand_limits(y = 0) +
  labs(title = "Fragment count normalized to secuencing depth per sample",
       x = "Donor",
       y = "Fragment count normalized to sequencing depth") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "none",
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray92"),
        panel.grid.minor = element_line(color = "gray92"),
        panel.border = element_rect(color = "gray92", fill = NA) )

# 2/3 fragment count - per sample, sorted by donor ID
isa_dataframe_merged %>%
  mutate(Sample = factor(Sample, levels = unique(Sample[order(donor_id, method)]))) %>%
  ggplot(aes(x = donor_id, y = normFragCount, fill = method)) +
  geom_boxplot(outlier.size = 0.2) +
  scale_fill_manual(values = c("mono" = "mediumpurple1", "multi" = "darkorange")) +
  expand_limits(y = 0) +
  labs(title = "Fragment count normalized to secuencing depth per sample",
       x = "Donor",
       y = "Fragment count normalized to sequencing depth") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray92"),
        panel.grid.minor = element_line(color = "gray92"),
        panel.border = element_rect(color = "gray92", fill = NA) )

# 3/3 fragment count - per method
isa_dataframe_merged %>%
  mutate(method = factor(method, levels = c("mono", "multi"))) %>%
  ggplot(aes(x = method, y = normFragCount, fill = method)) +
  geom_boxplot(outlier.size = 0.2) +
  scale_fill_manual(values = c("mono" = "mediumpurple1", "multi" = "darkorange")) +
  expand_limits(y = 0) +
  labs(title = "Fragment count normalized to secuencing depth per method",
       x = "Method",
       y = "Fragment count normalized to sequencing depth") +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray92"),
        panel.grid.minor = element_line(color = "gray92"),
        panel.border = element_rect(color = "gray92", fill = NA) ) +
  annotate("text", x = 1.5, y = 0.0011, label = "T-test, p-value = 0.000483", size = 4, color = "black")


## TSS enrichment ----

# 1/3 TSS enrichment - per sample, sorted by method
isa_dataframe_merged %>%
  mutate(Sample = factor(Sample, levels = unique(Sample[order(method, donor_id)]))) %>%
  ggplot(aes(x = donor_id, y = TSSEnrichment, fill = method)) +
  geom_boxplot(outlier.size = 0.2) +
  facet_wrap(~method, scales = "free_x") +
  scale_fill_manual(values = c("mono" = "mediumpurple1", "multi" = "darkorange")) +
  geom_hline(yintercept = 7, color = "red") +
  expand_limits(y = 0) +
  labs(title = "TSS enrichment per sample",
       x = "Donor",
       y = "TSS enrichment") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "none",
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray92"),
        panel.grid.minor = element_line(color = "gray92"),
        panel.border = element_rect(color = "gray92", fill = NA) )

# 2/3 TSS enrichment - per sample, sorted by donor ID
isa_dataframe_merged %>%
  mutate(Sample = factor(Sample, levels = unique(Sample[order(donor_id, method)]))) %>%
  ggplot(aes(x = donor_id, y = TSSEnrichment, fill = method)) +
  geom_boxplot(outlier.size = 0.2) +
  scale_fill_manual(values = c("mono" = "mediumpurple1", "multi" = "darkorange")) +
  geom_hline(yintercept = 7, color = "red") +
  expand_limits(y = 0) +
  labs(title = "TSS enrichment per sample",
       x = "Donor",
       y = "TSS enrichment") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray92"),
        panel.grid.minor = element_line(color = "gray92"),
        panel.border = element_rect(color = "gray92", fill = NA) )

# 3/3 TSS enrichment -  per method
isa_dataframe_merged %>%
  mutate(method = factor(method, levels = c("mono", "multi"))) %>%
  ggplot(aes(x = method, y = TSSEnrichment, fill = method)) +
  geom_boxplot(outlier.size = 0.2) +
  scale_fill_manual(values = c("mono" = "mediumpurple1", "multi" = "darkorange")) +
  geom_hline(yintercept = 7, color = "red") +
  expand_limits(y = 0) +
  labs(title = "TSS enrichment per method",
       x = "Method",
       y = "TSS enrichment") +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray92"),
        panel.grid.minor = element_line(color = "gray92"),
        panel.border = element_rect(color = "gray92", fill = NA)) +
  annotate("text", x = 1.5, y = 35, label = "Wilcox, p-value = 0.00169", size = 4, color = "black")


## blacklist ratio ----

# 1/3 BlacklistRatio - per sample, sorted by method
isa_dataframe_merged %>%
  mutate(Sample = factor(Sample, levels = unique(Sample[order(method, donor_id)]))) %>%
  ggplot(aes(x = donor_id, y = BlacklistRatio, fill = method)) +
  geom_boxplot(outlier.size = 0.2) +
  facet_wrap(~method, scales = "free_x") +
  scale_fill_manual(values = c("mono" = "mediumpurple1", "multi" = "darkorange")) +
  scale_y_log10() +
  expand_limits(y = 0) +
  labs(title = "Blacklist ratio per sample",
       x = "Donor",
       y = "Blacklist ratio (log10 scale)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "none",
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray92"),
        panel.grid.minor = element_line(color = "gray92"),
        panel.border = element_rect(color = "gray92", fill = NA) )

# 2/3 BlacklistRatio - per sample, sorted by donor ID
isa_dataframe_merged %>%
  mutate(Sample = factor(Sample, levels = unique(Sample[order(donor_id, method)]))) %>%
  ggplot(aes(x = donor_id, y = BlacklistRatio, fill = method)) +
  geom_boxplot(outlier.size = 0.2) +
  scale_fill_manual(values = c("mono" = "mediumpurple1", "multi" = "darkorange")) +
  scale_y_log10() +
  expand_limits(y = 0) +
  labs(title = "Blacklist ratio per sample",
       x = "Donor",
       y = "Blacklist ratio (log10 scale)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray92"),
        panel.grid.minor = element_line(color = "gray92"),
        panel.border = element_rect(color = "gray92", fill = NA) )

# 3/3 BlacklistRatio -  per method
isa_dataframe_merged %>%
  mutate(method = factor(method, levels = c("mono", "multi"))) %>%
  ggplot(aes(x = method, y = BlacklistRatio, fill = method)) +
  geom_boxplot(outlier.size = 0.2) +
  scale_fill_manual(values = c("mono" = "mediumpurple1", "multi" = "darkorange")) +
  scale_y_log10() +
  expand_limits(y = 0) +
  labs(title = "Blacklist ratio per method",
       x = "Method",
       y = "Blacklist ratio (log10 scale)") +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray92"),
        panel.grid.minor = element_line(color = "gray92"),
        panel.border = element_rect(color = "gray92", fill = NA)) +
  annotate("text", x = 1.5, y = 0.09, label = "T-test, p-value = 0.295", size = 4, color = "black")

## FRiP ----

archr_mono_df <- as.data.frame(getCellColData(archr_mono))
archr_multi_df <- as.data.frame(getCellColData(archr_multi))

archr_mono_df$FRIP_percent <- c(archr_mono_df$FRIP*100)
archr_multi_df$FRIP_percent <- c(archr_multi_df$FRIP*100)

# mono
archr_mono_df %>%
  mutate(Sample = factor(Sample, levels = unique(Sample[order(method, Sample)]))) %>%
  ggplot(aes(x = Sample, y = FRIP_percent, fill = method)) +
  geom_boxplot(outlier.size = 0.2) +
  scale_fill_manual(values = c("mono" = "mediumpurple1", "multi" = "darkorange")) +
  expand_limits(y = 0) +
  geom_hline(yintercept = 20, color = "red") +
  labs(title = "FRiP percentage for mono per sample",
       x = "Sample",
       y = "FRiP percentage") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "none",
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray92"),
        panel.grid.minor = element_line(color = "gray92"),
        panel.border = element_rect(color = "gray92", fill = NA) )

# multi
archr_multi_df %>%
  mutate(Sample = factor(Sample, levels = unique(Sample[order(method, Sample)]))) %>%
  ggplot(aes(x = Sample, y = FRIP_percent, fill = method)) +
  geom_boxplot(outlier.size = 0.2) +
  scale_fill_manual(values = c("mono" = "mediumpurple1", "multi" = "darkorange")) +
  expand_limits(y = 0) +
  geom_hline(yintercept = 20, color = "red") +
  labs(title = "FRiP percentage for multi per sample",
       x = "Sample",
       y = "FRiP percentage") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "none",
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray92"),
        panel.grid.minor = element_line(color = "gray92"),
        panel.border = element_rect(color = "gray92", fill = NA) )

# plotting them togehter

archr_combined_df <- rbind(archr_mono_df, archr_multi_df)

# adding donor id to dataframe
archr_combined_df_2 <- merge(archr_combined_df, meta_filt, by.x = "Sample", by.y = "id", all.x = TRUE)
archr_combined_df_2$method.y <- NULL

archr_combined_df_2 <- archr_combined_df_2 %>%
  dplyr::rename(method = method.x)

# 1/3 FRiP - per sample, by method
archr_combined_df_2 %>%
  mutate(Sample = factor(Sample, levels = unique(Sample[order(method, donor_id)]))) %>%
  ggplot(aes(x = donor_id, y = FRIP_percent, fill = method)) +
  geom_boxplot(outlier.size = 0.2) +
  facet_wrap(~method, scales = "free_x") +
  scale_fill_manual(values = c("mono" = "mediumpurple1", "multi" = "darkorange")) +
  geom_hline(yintercept = 20, color = "red") +
  expand_limits(y = 0) +
  labs(title = "FRiP in percentage per sample",
       x = "Donor",
       y = "FRiP (percentage)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "none",
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray92"),
        panel.grid.minor = element_line(color = "gray92"),
        panel.border = element_rect(color = "gray92", fill = NA) )

# 2/3 FRiP - per sample, sorted by donor ID
archr_combined_df_2 %>%
  mutate(Sample = factor(Sample, levels = unique(Sample[order(donor_id, method)]))) %>%
  ggplot(aes(x = donor_id, y = FRIP_percent, fill = method)) +
  geom_boxplot(outlier.size = 0.2) +
  scale_fill_manual(values = c("mono" = "mediumpurple1", "multi" = "darkorange")) +
  geom_hline(yintercept = 20, color = "red") +
  expand_limits(y = 0) +
  labs(title = "FRiP in percentage per sample",
       x = "Donor",
       y = "FRiP (percentage)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray92"),
        panel.grid.minor = element_line(color = "gray92"),
        panel.border = element_rect(color = "gray92", fill = NA) )

# 3/3 FRiP - per method
archr_combined_df_2 %>%
  mutate(method = factor(method, levels = c("mono", "multi"))) %>%
  ggplot(aes(x = method, y = FRIP_percent, fill = method)) +
  geom_boxplot(outlier.size = 0.2) +
  scale_fill_manual(values = c("mono" = "mediumpurple1", "multi" = "darkorange")) +
  geom_hline(yintercept = 20, color = "red") +
  expand_limits(y = 0) +
  labs(title = "FRiP in percentage per method",
       x = "Method",
       y = "FRiP (percentage)") +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray92"),
        panel.grid.minor = element_line(color = "gray92"),
        panel.border = element_rect(color = "gray92", fill = NA)) +
  annotate("text", x = 1.5, y = 65, label = "T-test, p-value = 0.144", size = 4, color = "black")


## poor quality nuclei ----

# 1/4 poor quality nuclei - per sample, sorted by method
cell_count_all %>%
  mutate(Sample = factor(Sample, levels = unique(Sample[order(method, donor_id)]))) %>%
  ggplot(aes(x = donor_id, y = abs(difference), fill = method)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~method, scales = "free_x") +
  scale_fill_manual(values = c("mono" = "mediumpurple1", "multi" = "darkorange")) +
  labs(
    x = "Donor",
    y = "'Bad' cells (%)",
    fill = "Method",
    title = "Percentage of poor quality nuclei per sample"
  ) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray92"),
        panel.grid.minor = element_line(color = "gray92"),
        panel.border = element_rect(color = "gray92", fill = NA) )

# 2/4 poor quality nuclei - per sample, by donor ID
cell_count_all %>%
  mutate(Sample = factor(Sample, levels = unique(Sample[order(donor_id, method)]))) %>%
  ggplot(aes(x = donor_id, y = abs(difference), fill = method)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("mono" = "mediumpurple1", "multi" = "darkorange")) +
  labs(
    x = "Donor",
    y = "Percentage of poor quality nuclei",
    fill = "Method",
    title = "Percentage of poor quality nuclei per sample"
  ) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray92"),
        panel.grid.minor = element_line(color = "gray92"),
        panel.border = element_rect(color = "gray92", fill = NA) )

# 3/4 poor quality nuclei - per method
cell_count_all %>%
  mutate(method = factor(method, levels = c("mono", "multi"))) %>%
  ggplot(aes(x = method, y = abs(difference), fill = method)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("mono" = "mediumpurple1", "multi" = "darkorange")) +
  labs(
    x = "Donor",
    y = "'Bad' cells (%)",
    fill = "Method",
    title = "Percentage of poor quality nuclei per method"
  ) +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray92"),
        panel.grid.minor = element_line(color = "gray92"),
        panel.border = element_rect(color = "gray92", fill = NA) ) +
  annotate("text", x = 1.5, y = 70, label = "T-test, p-value = 0.0297", size = 4, color = "black")

# 4/4 poor quality nuclei - per method (mean)
cell_count_all %>%
  mutate(difference_abs = abs(difference)) %>% # Convert difference to absolute values
  group_by(method) %>% # Group by method
  summarize(mean_difference = mean(difference_abs, na.rm = TRUE)) %>% # Calculate mean for each method
  ggplot(aes(x = method, y = mean_difference, fill = method)) +
  geom_bar(stat = "identity", position = "dodge") + # Create a bar plot
  scale_fill_manual(values = c("mono" = "mediumpurple1", "multi" = "darkorange")) +
  labs(
    x = "Method",
    y = "percentage of poor quality nuceli (mean)",
    fill = "Method",
    title = "Percentage of poor quality nuclei per method"
  ) +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray92"),
        panel.grid.minor = element_line(color = "gray92"),
        panel.border = element_rect(color = "gray92", fill = NA)) +
  annotate("text", x = 1.5, y = 40, label = "T-test, p-value = 0.0297", size = 4, color = "black")
