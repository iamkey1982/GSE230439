#!/usr/bin/Rscript

###Set VERSION
VERSION = "1.0.1"

args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 1 && (args[1] == "-v" | args[1] == "--version")){
  message("scRNA VDJ pipeline\nVersion: \n\t", VERSION, "\n")
  quit(status = 0)
}

seurat_script <<- "/home/project/Single_cells_project/scrna_project/pipeline/vdj_script/seurat_bcr_single"
scrna_project_path <<- "/home/project/Single_cells_project/scrna_project"
analysis_path <<- c("1_Cellranger", "2_Clonotype_analysis", "3_Scrna_bcr_analysis")


args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 1 && (args[1] == "-v" | args[1] == "--version")){
  message("scRNA seurat BCR analysis pipeline\nVersion: \n\t", VERSION, "\n")
  quit(status = 0)
}

suppressPackageStartupMessages(library(optparse))      ## Options
suppressPackageStartupMessages(library(futile.logger)) ## logger
suppressPackageStartupMessages(library(future.apply))
suppressPackageStartupMessages(library(yaml))
suppressPackageStartupMessages(library(doParallel))

#--- save the parameter into R object
AllOptions <- function(){
  parser <- OptionParser()
  parser <- add_option(parser, c("-n", "--numberofcores"), type="integer", default=1,
                       help="Parallel Run Number of cores [default %default]",
                       metavar="number")
  
  parser <- add_option(parser, c("--MaxMemMega"), type="integer", default=20000,
                       help="Parallel Max Memory size megabytes [default %default]",
                       metavar="number")
  
  parser <- add_option(parser, c("-c", "--configfile"), type="character", default="config.yaml",
                       help="Config file for the run, input files and executing plan settings  [default %default]",
                       metavar="character")
  
  parser <- add_option(parser, c("-l", "--logname"), type="character", default="log",
                       help="add a unique name for log files [default %default]",
                       metavar="character")
  return(parser)
}

parser <- AllOptions()
#debug(parse_args)
pa <- parse_args(parser)

# 程序运行CPU核数
WORKER_NUM            = pa$numberofcores
# 并行最大内存大小(M)
MAXMEMMEGA            = pa$MaxMemMega
f = function(x){if(x == ""){ return("")}else{ return("-")}}
LOGNAME               = paste0(pa$logname, f(pa$logname))

registerDoParallel(cores=WORKER_NUM)
# Set up future for parallelization
plan("multiprocess", workers = WORKER_NUM)
options(future.globals.maxSize = MAXMEMMEGA * 1024^2)

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

## configuration file
if(!file.exists(pa$configfile)){
  logger.error("Please check the config file!")
  print("Please check the config file!")
  parse_args(parser, args = c("--help"))
}


yaml.list <- vector("list")
yaml.list <- yaml.load_file(pa$configfile) 
project_code <- gsub("^\\s+|\\s+$", "", yaml.list$project_code)
project_specie <- tolower(gsub("^\\s+|\\s+$", "", yaml.list$project_specie))
project_dir <- file.path(scrna_project_path, project_code)

source(file.path(seurat_script,'vdj_funcs.R'), encoding = 'UTF-8')
logger.info("----------------------Project: %s --------------------------", project_code)

##--------------keep variables----------------------
## project_work_dir <- file.path(project_dir, paste(project_code, as.character(format(Sys.Date(), "%Y%m%d")), sep = "_"))
project_work_dir <- file.path(project_dir, paste(project_code, 'vdj', 'bcr', 'report', 'single', sep = "_"))
project_work_data_dir <- file.path(project_work_dir, 'Rdata')
makedirs(c(project_work_dir, project_work_data_dir))


## sample list summary
sample_data_src = vector("list")
samplelist <- list.files(file.path(project_dir, 'data', 'bcr'))


for (sample in samplelist){
  print(sample)
  sample_analysis_path <- file.path(project_work_dir, sample, analysis_path)
  sample_data_src[sample] = file.path(project_dir, 'data', 'bcr', sample)
  makedirs(sample_analysis_path)
}


## project work directory
logger.info("project work directory: %s", project_work_dir)


#################################################
# 1_Cellranger
#################################################
print("start 1_Cellranger analysis...")
for(i in seq_along(sample_data_src)){
  sample = names(sample_data_src)[i]
  from_dir = sample_data_src[[sample]]
  sample_cellranger_path <- file.path(project_work_dir, sample, analysis_path[1])
  sapply(list.files(from_dir),function(x){file.copy(paste(from_dir,x,sep="/"),sample_cellranger_path,recursive = T)})
}


#################################################
# 2_data_qc_statistics 
#################################################
print("start 2_Clonotype_abundance analysis...")
source(file.path(seurat_script,'02_Clonotype_abundance.R'), encoding = 'UTF-8')
Clonotype_abundance()


#################################################
# 3_scrna_tcr_analysis
#################################################
print("start 3_scrna_tcr_analysis...")
source(file.path(seurat_script,'03_scrna_tcr_analysis.R'), encoding = 'UTF-8')
scrna_tcr_analysis()


