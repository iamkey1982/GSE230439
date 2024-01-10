#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~seurat_multiple_normalization~~~~~~~~~~~~#
#~~~~~Normalization+CCA整合~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#!/usr/bin/Rscript
### Set VERSION
VERSION = "1.0.1"
### Computational Parameters
my.seed <- 202106L

# single-sample
analysis_path <<- c("1_Cellranger", "2_data_qc_statistics", "3_Cell_clustering_analysis", "4_marker_gene_analysis", "5_celltype_annotation", "6_cell_proportion_among_samples", "7_cell_cycle_scoring", "8_Go_KEGG_FUN")
CLUSTER_RESOLUTION_RANGE <- seq(0.1, 1, 0.1)


suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(futile.logger))
suppressPackageStartupMessages(library(future))
suppressPackageStartupMessages(library(future.apply))
suppressPackageStartupMessages(library(yaml))
suppressPackageStartupMessages(library(doParallel))

#### Read parameters ####
#--- save the parameter into R object
AllOptions <- function(){
  parser <- OptionParser()
    ### 1. Project path
	parser <- add_option(parser, c("--project.path"), type="character", metavar="character",
                       help="Path to Project folder", default="/home/project/Single_cells_project/scrna_project")
					   
	### 2. pipeline path
	parser <- add_option(parser, c("--pipeline.path"), type="character", metavar="character",
                       help="Path to Project folder",
                       default="/home/project/Single_cells_project/scrna_project/pipeline/seurat_script/seurat_multiple_normalization")

	### 3. Computational Parameters
	parser <- add_option(parser, c("-n", "--numberofcores"), type="integer", default=1,
						   help="Parallel Run Number of cores [default %default]",
						   metavar="number")
	parser <- add_option(parser, c("--MaxMemMega"), type="integer", default=20000,
						   help="Parallel Max Memory size megabytes [default %default]",
						   metavar="number")

	### 4. Project configfile
	parser <- add_option(parser, c("-c", "--configfile"), type="character", default="config.yaml",
                       help="Config file for the run, input files and executing plan settings  [default %default]",
                       metavar="character")
					   
	### 5. Project logfile
	parser <- add_option(parser, c("-l", "--logname"), type="character", default="log",
                       help="add a unique name for log files [default %default]",
                       metavar="character")
	return(parser)
}


parser <- AllOptions()
#debug(parse_args)
args <- parse_args(parser)

#### Formatting Parameters ####
#convert "NULL"/"FALSE"/"TRUE" (in character) into NULL/FALSE/TRUE
for (i in names(args$options)){
  if (toupper(args$options[i]) == "NULL") { args$options[i] <- NULL
  } else if (toupper(args$options[i]) == "FALSE") { args$options[i] <- FALSE
  } else if (toupper(args$options[i]) == "TRUE") { args$options[i] <- TRUE
  }
}

#### Get Paramaters ####
### 1.Project path
scrna_project_path <- args$project.path
### 2.pipeline path
seurat_pipeline_path <-  args$pipeline.path
### 3.Computational Parameters
WORKER_NUM <- if (!is.null(args$numberofcores)) as.numeric(args$numberofcores)
MAXMEMMEGA <- if (!is.null(args$MaxMemMega)) as.numeric(args$MaxMemMega)
registerDoParallel(cores=WORKER_NUM)
### 3.1 Set up future for parallelization
plan("multicore", workers = WORKER_NUM)
# plan(strategy = "multicore", workers = WORKER_NUM)
options(future.globals.maxSize = MAXMEMMEGA * 1024^2)
### 4.Project configfile
if(!file.exists(args$configfile)){
  logger.error("Please check the config file!")
  print("Please check the config file!")
  parse_args(parser, args = c("--help"))
}
yaml.list <- vector("list")
yaml.list <- yaml.load_file(args$configfile) 
project_code <- gsub("^\\s+|\\s+$", "", yaml.list$project_code)
project_specie <- tolower(gsub("^\\s+|\\s+$", "", yaml.list$project_specie))
project_dir <- file.path(scrna_project_path, project_code)
### 5.logfile
f = function(x){if(x == ""){ return("")}else{ return("-")}}
LOGNAME <- paste0(args$logname, f(args$logname))
if(!file.exists("logs")){
  dir.create("logs")
}
cur_date <- as.character(Sys.Date())
errorLog  =  file.path("logs", sprintf("%sERROR-%s.log", LOGNAME, cur_date))
warnLog   =  file.path("logs", sprintf("%sWARN-%s.log", LOGNAME, cur_date))
infoLog   =  file.path("logs",  sprintf("%sINFO-%s.log", LOGNAME, cur_date))
traceLog   =  file.path("logs",  sprintf("%sTRAC-%s.log", LOGNAME, cur_date))

invisible(flog.logger("error", ERROR, appender.file(errorLog)))
invisible(flog.logger("warn", WARN, appender.file(warnLog)))
invisible(flog.logger("info", INFO, appender.file(infoLog)))
invisible(flog.logger("trace", TRACE, appender.file(traceLog)))
invisible(flog.appender(appender.console(), name = "ROOT"))

logger.info <- function(msg, ...) {
  flog.info(msg, ..., name = "ROOT")
  flog.info(msg, ..., name = "info")
}

logger.warn <- function(msg, ...) {
  flog.warn(msg, ..., name = "ROOT")
  flog.warn(msg, ..., name = "info")
  flog.warn(msg, ..., name = "warn")
}

logger.error <- function(msg, ...) {
  flog.error(msg, ..., name = "ROOT")
  flog.error(msg, ..., name = "info")
  flog.error(msg, ..., name = "warn")
  flog.error(msg, ..., name = "error")
  flog.error(msg, ..., name = "trace")
}


#### Get the local pipeline path  ####
if(is.null(seurat_pipeline_path)) stop("--pipeline.path parameter must be set!")
#### Check non-optional parameters ####
if (is.null(scrna_project_path)) stop("project.path parameter can't be empty!")
if (is.null(project_specie)) stop("project_specie parameter can't be empty!")

### Clean
rm(AllOptions,parser,args)

source(file.path(seurat_pipeline_path, 'sc-funcs.R'), encoding = 'UTF-8')
logger.info("----------------------Project: %s --------------------------", project_code)


##--------------keep variables----------------------
## project_work_dir <- file.path(project_dir, paste(project_code, as.character(format(Sys.Date(), "%Y%m%d")), sep = "_"))
project_work_dir <- file.path(project_dir, paste(project_code, 'rna_seq','report_normalization', sep = "_"))
project_work_data_dir <- file.path(project_work_dir, 'Rdata')
makedirs(c(project_work_dir, project_work_data_dir))
# create analysis report 
makedirs(file.path(project_work_dir, analysis_path))

## sample list summary
sample_data_src = vector("list")
samplelist <- list.files(file.path(project_dir, 'data', 'rna_seq'))
for (sample in samplelist){
  sample_analysis_path <- file.path(project_work_dir, sample, analysis_path)
  sample_data_src[sample] = file.path(project_dir, 'data', 'rna_seq', sample)
}

## project work directory
logger.info("project work directory: %s", project_work_dir)

begin <- Sys.time()
#################################################
# 1_Cellranger
#################################################
print("start 1_Cellranger analysis...")
for(i in seq_along(sample_data_src)){
  sample <- names(sample_data_src)[i]
  from_dir <- sample_data_src[[sample]]
  sample_cellranger_path <- file.path(project_work_dir, analysis_path[1])
  file.copy(from_dir,sample_cellranger_path, recursive=T)
}

#################################################
# 2_data_qc_statistics 
#################################################
print("start 2_data_qc_statistics analysis...")
source(file.path(seurat_pipeline_path,'02_data_qc_statistics.R'), encoding = 'UTF-8')
if (!file.exists(file.path(project_work_data_dir, 'anchors.Rdata'))){
  data_qc_statistics()
}


#################################################
# 3_Cell_clustering_analysis
#################################################
print("start 3_Cell_clustering_analysis...")
source(file.path(seurat_pipeline_path,'03_cell_clustering_analysis.R'), encoding = 'UTF-8')
if (!file.exists(file.path(project_work_data_dir, 'Cluster_all.Rdata'))){
  cell_clustering_analysis()
}


#################################################
# 4_marker_gene_analysis
#################################################
print("start 4_marker_gene_analysis...")
source(file.path(seurat_pipeline_path,'04_marker_gene_analysis.R'), encoding = 'UTF-8')
marker_gene_analysis()


####################################
# 5_celltype_annotation
####################################
print("start 05_celltype_annotation analysis...")
source(file.path(seurat_pipeline_path,'05_celltype_annotation.R'), encoding = 'UTF-8')
celltype_annotation()


####################################
# 6_cell_proportion_among_samples
####################################
print("start 6_cell_proportion_among_samples...")
source(file.path(seurat_pipeline_path,'06_cell_proportion_among_samples.R'), encoding = 'UTF-8')
cell_proportion_among_samples()


####################################
# 7_cell_cycle_scoring
####################################
print("start 7_cell_cycle_scoring...")
source(file.path(seurat_pipeline_path,'07_cell_cycle_scoring.R'), encoding = 'UTF-8')
cell_cycle_scoring()

####################################
# 8_Go_KEGG_FUN.R
####################################
print("start 08_Go_KEGG_FUN...")
source(file.path(seurat_pipeline_path,'08_Go_KEGG_FUN.R'), encoding = 'UTF-8')
Go_KEGG_FUN()


end <- Sys.time()
print(sprintf('分析累计用时：%s mins', difftime(end, begin, units = "mins")))











