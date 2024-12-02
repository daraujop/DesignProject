options(repos = c(CRAN = "https://cloud.r-project.org/"))

getOption("repos")

setRepositories(ind = c(1, 2, 3, 4, 5, 6, 8))  # CRAN, BioC software, annotation, experiment, etc.

install.packages("devtools")

library(devtools)

# Then install package with
devtools::install_github("biosurf/cyCombine") # install devtools

library(cyCombine)

# Install flowCore if not already installed
if (!requireNamespace("flowCore", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("flowCore")
}

# Load the library
library(flowCore)

install.packages("tidyverse")  # Installs the entire tidyverse suite
library(tidyverse)             # Loads tidyverse, which includes `%>%`

setwd("C:/Users/dhama/DESIGN_PROJECT")
getwd()  # Confirm the current working directory
list.files()

data_dir <- "Data_own"
list.files(data_dir, pattern = "\\.fcs$")


#This is just to load the data in, but you don't actually have to do this
# Load the fcs file
fcs_file1 <- read.FCS("Data_own/Day1_KO_Control2_QC.fcs", transformation = FALSE)

# View the first few rows of data
head(exprs(fcs_file1))

# Check the column names to see the available parameters
colnames(fcs_file1)

#### Compile FCS files ---- (first function)


#' Compile all .fcs files in a directory to a flowset
#'
#'  A simple function to compile a directory of FCS files into a flowSet object
#'  using flowCore::read.flowSet().
#'  Use the pattern argument to select only a subset of FCS files.
#'
#' @param data_dir Directory containing the .fcs files
#' @param pattern The pattern to use to find the files in the folder
#' @inheritParams flowCore::read.FCS
#' @family dataprep
#' @examples
#' \dontrun{
#' fcs <- compile_fcs(data_dir = "_data/raw", pattern = "\\.fcs")
#' }
#' @export

# flowCore::read.flowSet("Day1_KO_Control2_QC.fcs")

compile_fcs <- function(
    data_dir,
    pattern = "\\.fcs",
    column.pattern = NULL,
    invert.pattern = FALSE)
{
  
  
  # Error checking
  if (data_dir %>% endsWith("/")) {
    data_dir <- stringr::str_sub(data_dir, end = -2)
  }
  
  # Specifying files to use
  files <- list.files(data_dir,
                      pattern = pattern,
                      recursive = FALSE,
                      full.names = TRUE) %>%
    sort()
  if (length(files) == 0) stop("No files found in folder \"", data_dir, "\"")
  
  # Read the data files
  message(paste("Reading", length(files), "files to a flowSet.."))
  fcs_raw <- files %>%
    flowCore::read.flowSet(transformation = FALSE,
                           truncate_max_range = FALSE,
                           emptyValue = FALSE,
                           column.pattern = column.pattern,
                           invert.pattern = invert.pattern)
  
  return(fcs_raw)
}

# Apply the first function
fcs<-compile_fcs(data_dir = "Data_own",pattern = "\\.fcs")

# Needed libraries for the next function
library(dplyr)
library(tibble)
library(stringr)
library(readxl)
library(readr)
library(flowCore)
library(Biobase)
library(purrr)

#' Convert a flowSet into a tibble (second function)
#'
#' Use this function to convert a flowSet into the tibble
#'  object that the remaining functions in cyCombine relies on.
#'  A tibble is a Tidyverse implementation of a data.frame and can be treated as a such.
#'  The majority of arguments revolves adding relevant info from the metadata file/object.
#'  The panel argument is included to adjust the output column names using a panel with channel and antigen columns.
#'  Bear in mind the column names will be altered with the following:
#'  \code{stringr::str_remove_all("^\\d+\[A-Za-z\]+_") %>% stringr::str_remove_all("\[ _-\]")}.
#'
#'
#'
#' @param flowset The flowset to convert
#' @param metadata Optional: Can be either a filename or data.frame of the metadata file. Please give the full path from working directory to metadata file
#' @param sample_ids Optional: If a character, it should be the sample column in the metadata. If its a vector, it should have the same length as the total flowset. If NULL, sample ids will be the file names. If a single value, all rows will be assigned this value.
#' @param batch_ids Optional: If a character, it should be the column in the metadata containing the batch ids. If its a vector, it should have the same length as the total flowset. If a single value, all rows will be assigned this value.
#' @param filename_col Optional: The column in the metadata containing the fcs filenames. Needed if metadata is given, but sample_ids is not
#' @param condition Optional: The column in the metadata containing the condition. Will be used as the covariate in ComBat, but can be specified later. You may use this to add a different column of choice, in case you want to use a custom column in the ComBat model matrix.
#' @param anchor Experimental: The column in the metadata referencing the anchor samples (control references). Will be used as a covariate in ComBat, if specified. Please be aware that this column may be confounded with the condition column. You may use this to add a different column of choice, in case you want to use a custom column in the ComBat model matrix. You may use a custom column name, but it is good practice to add the name to the 'non_markers' object exported by cyCombine, to reduce the risk of unexpected errors.
#' @param down_sample If TRUE, the output will be down-sampled to size sample_size
#' @param sample_size The size to down-sample to. If a non-random sampling type is used and a group contains fewer cells than the sample_size, all cells of that group will be used.
#' @param sampling_type The type of down-sampling to use. "random" to randomly select cells across the entire dataset, "batch_ids" to sample evenly (sample_size) from each batch, or "sample_ids" sample evenly (sample_size) from each sample.
#' @param seed The seed to use for down-sampling
#' @param clean_colnames (Default: TRUE). A logical defining whether column names should be cleaned or not. Cleaning involves removing isotope tags, spaces, dashes, underscores, and all bracket types.
#' @param panel Optional: Panel as a filename or data.frame. Is used to define colnames from the panel_antigen column
#' @param panel_channel Optional: Only used if panel is given. It is the column name in the panel data.frame that contains the channel names
#' @param panel_antigen Optional: Only used if panel is given. It is the column name in the panel data.frame that contains the antigen names
#' @family dataprep
#' @examples
#' \dontrun{
#' df <- convert_flowset(flowset = flowset,
#'  metadata = file.path(data_dir, "metadata.csv"),
#'  filename_col = "FCS_files",
#'  sample_ids = "sample_id",
#'  batch_ids = "batch_ids",
#'  down_sample = FALSE)
#'  }
#' @export
convert_flowset <- function(flowset,
                            metadata = NULL,
                            filename_col = "filename",
                            sample_ids = NULL,
                            batch_ids = NULL,
                            condition = NULL,
                            anchor = NULL,
                            down_sample = TRUE,
                            sample_size = 500000,
                            sampling_type = "random",
                            seed = 473,
                            clean_colnames = TRUE,
                            panel = NULL,
                            panel_channel = "fcs_colname",
                            panel_antigen = "antigen"){
  # Extract information necessary both with and without included metadata
  ## FlowSet row numbers
  nrows <- flowCore::fsApply(flowset, nrow)
  
  ## File names from flowset
  files <- flowset@phenoData %>%
    rownames() %>%
    basename()
  # Add metadata information (long)
  if(!is.null(metadata)){
    # Error handling
    if(is.null(filename_col)){
      stop("Please specify a filename_col.")
    }
    if("character" %in% class(metadata)){
      if(!file.exists(file.path(metadata))){
        stop("File \"", file.path(metadata), "\" was not found")
      }
      
      # Get metadata
      if(endsWith(metadata, suffix = ".xlsx")){
        metadata <- suppressMessages(file.path(metadata) %>%
                                       readxl::read_xlsx())
      } else if(endsWith(metadata, suffix = ".csv")){
        metadata <- suppressMessages(file.path(metadata) %>%
                                       readr::read_csv())
      } else {
        stop("Sorry, the given metadata is not in a supported format. Please use a .xlsx or .csv file.\n",
             "Alternatively, a data.frame of the metadata can be used.")
      }
    }
    
    # Check for errors in metadata columns
    md_cols <- colnames(metadata)
    cyCombine:::check_colname(md_cols, filename_col)
    
    
    # Extract info from metadata
    if(!endsWith(tolower(metadata[[filename_col]][1]), ".fcs")){
      metadata[[filename_col]] <- paste0(metadata[[filename_col]], ".fcs")
    }
    # Check that all metadata rows has a file
    if(any(metadata[[filename_col]] %!in% files)){
      missing_files <- metadata[[filename_col]][metadata[[filename_col]] %!in% files]
      warning("The following samples in the metadata were not found in the provided folder and will be ignored.\n", stringr::str_c(missing_files, collapse = ", "))
      # Remove files from metadata
      metadata <- metadata[metadata[[filename_col]] %in% files,]
    }
    
    # Check that all files are represented in metadata
    if(any(files %!in% metadata[[filename_col]])){
      missing_files <- files[files %!in% metadata[[filename_col]]]
      warning("The samples were not found in the metadata file and will be ignored.\n", stringr::str_c(missing_files, collapse = ", "))
      files <- files[files %in% metadata[[filename_col]]]
      nrows <- nrows[rownames(nrows) %in% files,] %>% as.matrix()
      flowset <- flowset[files]
    }
    
    # Get sample ids
    if (is.null(sample_ids)){
      sample_ids <- files %>%
        stringr::str_remove(".fcs") %>%
        stringr::str_remove(".FCS") %>%
        rep(nrows)
    } else if (length(sample_ids) == 1){
      cyCombine:::check_colname(md_cols, sample_ids)
      sample_ids <- metadata[[sample_ids]][match(files, metadata[[filename_col]])] %>%
        stringr::str_remove(".fcs") %>%
        stringr::str_remove(".FCS") %>%
        rep(nrows)
    }
    
    # Get batch ids
    if (!is.null(batch_ids)){
      if(length(batch_ids) == 1){
        cyCombine:::check_colname(md_cols, batch_ids)
        batch_ids <- metadata[[batch_ids]][match(files, metadata[[filename_col]])] %>%
          as.factor() %>%
          rep(nrows)
      }
    }
    # Get condition
    if(!is.null(condition)){
      if(length(condition == 1)){
        cyCombine:::check_colname(md_cols, condition)
        condition <- metadata[[condition]][match(files, metadata[[filename_col]])] %>%
          as.factor() %>%
          rep(nrows)
      }
    }
    # Get anchor
    if(!is.null(anchor)){
      if(length(anchor == 1)){
        cyCombine:::check_colname(md_cols, anchor)
        anchor <- metadata[[anchor]][match(files, metadata[[filename_col]])] %>%
          as.factor() %>%
          rep(nrows)
      }
    }
    
    
  } else{ # If no metadata given
    if(is.null(sample_ids)){# Get sample ids from filenames
      sample_ids <- files %>%
        stringr::str_remove(".fcs") %>%
        stringr::str_remove(".FCS") %>%
        rep(nrows)
    }
  }
  ## To down sample within fsApply
  tot_nrows <- sum(nrows)
  ids <- 1:tot_nrows
  # Down sampling setup
  if (down_sample){
    set.seed(seed)
    # Sorting here enables major resource savings when down-sampling
    # For non-random sampling, the dplyr::slice allows down-sampling to group size if there are less cells than sample_size.
    if (sampling_type == "random"){
      message("Down sampling to ", sample_size, " cells")
      sample <- sample(ids, sample_size) %>%
        sort()
    } else if ((sampling_type == "batch_ids") && !is.null(batch_ids)){ # even down-sampling from batches
      message(paste("Down sampling to", sample_size, "cells per batch"))
      sample <- tibble::tibble(batch_ids, ids) %>%
        dplyr::group_by(batch_ids) %>%
        dplyr::slice(sample(dplyr::n(), min(sample_size, dplyr::n()))) %>%
        dplyr::pull(ids) %>%
        sort()
      tiny_batches <- table(batch_ids) < sample_size
      if(any(tiny_batches)) message("Please be aware that batches the following batches contained less than ",
                                    sample_size, " cells:\n",
                                    stringr::str_c(names(tiny_batches), collapse = ", "))
    } else{ # Even down-sampling from samples
      message(paste("Down sampling to", sample_size, "cells per sample"))
      sample <- tibble::tibble(sample_ids, ids) %>%
        dplyr::group_by(sample_ids) %>%
        dplyr::slice(sample(dplyr::n(), min(sample_size, dplyr::n()))) %>%
        dplyr::pull(ids) %>%
        sort()
      tiny_samples <- table(sample_ids) < sample_size
      if(any(tiny_samples)) message("Please be aware that batches the following batches contained less than",
                                    sample_size, "cells:\n",
                                    stringr::str_c(names(tiny_samples), collapse = ", "))
    }
    
    
    # Down-sample metadata columns
    ids <- ids[sample]
    if(!is.null(sample_ids) && length(sample_ids) > 1) sample_ids <- sample_ids[sample]
    if(!is.null(batch_ids) && length(batch_ids) > 1) batch_ids <- batch_ids[sample]
    if(!is.null(condition) && length(condition) > 1) condition <- condition[sample]
    if(!is.null(anchor) && length(anchor) > 1) anchor <- anchor[sample]
  }
  
  message("Extracting expression data..")
  fcs_data <- flowset %>%
    purrr::when(down_sample ~ flowCore::fsApply(., cyCombine:::fcs_sample,
                                                sample = sample,
                                                nrows = nrows),
                ~ flowCore::fsApply(., Biobase::exprs)) %>%
    tibble::as_tibble() %>%
    dplyr::mutate(id = ids) %>%
    dplyr::select(id, dplyr::everything())
  
  # Clean column names
  if (!is.null(panel)){
    if("character" %in% class(panel)){
      if(endsWith(panel, suffix = ".xlsx")){
        panel <- suppressMessages(file.path(panel) %>%
                                    readxl::read_xlsx())
      } else if(endsWith(panel, suffix = ".csv")){
        panel <- suppressMessages(file.path(panel) %>%
                                    readr::read_csv())
      } else {
        stop("Sorry, the panel file is not in a supported format. Please use a .xlsx or .csv file.")
      }
    }
    cols <- match(colnames(fcs_data), panel[[panel_channel]]) %>%
      .[!is.na(.)]
    col_names <- panel[[panel_antigen]][cols]
    
    fcs_data <- fcs_data %>%
      dplyr::select(id, dplyr::all_of(panel[[panel_channel]][cols]))
  }else{
    col_names <- flowset[[1]] %>%
      flowCore::parameters() %>%
      Biobase::pData() %>%
      dplyr::pull(desc)
  }
  if(clean_colnames) {
    col_names <- col_names %>%
      stringr::str_remove_all("^\\d+[A-Za-z]+_") %>%
      stringr::str_remove_all("[- _\\[\\](){}]")
  }
  colnames(fcs_data) <- c("id", col_names)
  
  # Add optional columns
  if(!is.null(sample_ids)) fcs_data$sample <- sample_ids
  if(!is.null(batch_ids)) fcs_data$batch <- batch_ids
  if(!is.null(condition)) fcs_data$condition <- condition
  if(!is.null(anchor)) fcs_data$anchor <- anchor
  
  message("Your flowset is now converted into a dataframe.")
  return(fcs_data)
}


# apply second function
df <- convert_flowset(flowset = fcs,
                      #metadata = file.path(data_dir, "metadata.csv"),
                      filename_col = "FCS_files",
                      sample_ids = "sample_id",
                      batch_ids = "batch_ids",
                      down_sample = FALSE)

summary(df)

# I don't think we need to use this function
#' Extract from a flowset given a sample of indices
#' @noRd

fcs_sample <- function(flowframe, sample, nrows, seed = 473){
  
  # Determine which flowframe was given
  ff_name <- flowCore::keyword(flowframe)$FILENAME %>%
    basename()
  ff_number <- which(rownames(nrows) == ff_name)
  
  # Down-sample based on accumulated nrows (ensures the correct rows are extracted from each flowframe)
  nrows_acc <- c(0, nrows %>%
                   purrr::accumulate(`+`))
  ff_sample <- sample - nrows_acc[ff_number]
  ff_sample <- ff_sample[ff_sample > 0]
  ff_sample <- ff_sample[ff_sample <= nrows[ff_number]]
  
  # Extract expression data
  ff <- flowframe %>%
    Biobase::exprs()
  ff <- ff[ff_sample, ]
  
  return(ff)
  
  
  
  #### Data transformation ---- (3th function)
  
  #' Transform data using asinh
  #'
  #' Inverse sine transformation of a tibble.
  #'  This function can also de-randomize data.
  #'
  #' @param df The tibble to transform
  #' @param markers The markers to transform on
  #' @param cofactor The cofactor to use when transforming
  #' @param derand Derandomize. Should be TRUE for CyTOF data, otherwise FALSE.
  #' @param .keep Keep all channels. If FALSE, channels that are not transformed are removed
  #' @param reverse Reverses the asinh transformation if TRUE
  #' @family dataprep
  #' @examples
  #' \dontrun{
  #' uncorrected <- df %>%
  #'   transform_asinh(markers = markers)
  #'   }
  #' @export
  
  transform_asinh <- function(df,
                              markers = NULL,
                              cofactor = 5,
                              derand = TRUE,
                              .keep = FALSE,
                              reverse = FALSE){
    if(is.null(markers)){
      markers <- df %>%
        cyCombine::get_markers()
    }
    # Use global non_markers if available
    if(!is.null(.GlobalEnv$non_markers)) non_markers <- .GlobalEnv$non_markers
    
    if(any(markers %!in% colnames(df))){
      mes <- stringr::str_c("Not all given markers are in the data.\nCheck if the markers contain a _ or -:",
                            stringr::str_c(markers, collapse = ", "),
                            "Columns:",
                            stringr::str_c(colnames(df), collapse = ", "),
                            sep = "\n"
      )
      stop(mes)
    } else if(.keep & any(colnames(df) %!in% unique(colnames(df)))){
      stop("Your data contains non-unique column names. Please ensure they are unique. The column names are: ", stringr::str_c(colnames(df), collapse = ", "))
    }
    message("Transforming data using asinh with a cofactor of ", cofactor, "..")
    transformed <- df %>%
      purrr::when(.keep ~ .,
                  ~ dplyr::select_if(., colnames(.) %in% c(markers, non_markers))) %>%
      # Transform all data on those markers
      dplyr::mutate(dplyr::across(dplyr::all_of(markers),
                                  .fns = function(x){
                                    if(derand & !reverse) x <- ceiling(x)
                                    if(reverse) sinh(x)*cofactor else asinh(x/cofactor)
                                  }))
    return(transformed)
  }
  
  # Trying thrid function, i dont get what the markers should be so i'm stuck on this
  
  
  # Step 1: Verify column names and markers
  parameters <- pData(parameters(fcs[[1]]))
  markers <- parameters$name[!is.na(parameters$desc) & parameters$desc != ""]
  print(colnames(df))
  print(markers)  # Ensure this prints valid marker names
  
  
  # Step 1: Handle NA column names in df
  # Replace NA column names with temporary placeholders
  colnames(df)[is.na(colnames(df))] <- paste0("NA_", seq_along(which(is.na(colnames(df)))))
  print("Updated column names:")
  print(colnames(df))
  
  # Step 1: Handle NA column names
  colnames(df)[is.na(colnames(df))] <- paste0("NA_", seq_along(which(is.na(colnames(df)))))
  
  # Step 2: Standardize column and marker names
  colnames(df) <- gsub("_", "-", colnames(df))
  colnames(df) <- gsub("#", "", colnames(df))
  markers <- gsub("_", "-", markers)
  markers <- gsub("#", "", markers)
  
  # Step 3: Map markers to columns in df
  markers <- markers[markers %in% colnames(df)]
  missing_markers <- setdiff(markers, colnames(df))
  if (length(missing_markers) > 0) {
    warning("The following markers are missing in df: ", paste(missing_markers, collapse = ", "))
  }
  
  # Step 4: Define non-marker columns
  non_markers <- colnames(df)[!colnames(df) %in% markers]
  
  # Step 5: Apply asinh transformation
  `%!in%` <- Negate(`%in%`)
  uncorrected <- df %>%
    transform_asinh(markers = markers, cofactor = 5)
  
  # Inspect the result
  head(uncorrected)
  
  
  
  
  
  #write.csv(uncorrected, "uncorrected_transformed.csv", row.names = FALSE)

  
  #### Wrapper function ----
  
  
  
  #' Prepare a directory of .fcs files
  #'
  #' This is a wrapper function that takes you from a directory of .fcs files or a flowset to a transformed tibble.
  #'
  #'
  #' @param flowset Optional: Prepare a flowset instead of a directory of fcs files
  #' @inheritParams compile_fcs
  #' @inheritParams convert_flowset
  #' @inheritParams transform_asinh
  #' @param transform If TRUE, the data will be transformed; if FALSE, it will not.
  #' @family dataprep
  #' @examples
  #' \dontrun{
  #' uncorrected <- data_dir %>%
  #'   prepare_data(metadata = "metadata.csv",
  #'   markers = markers,
  #'   filename_col = "FCS_name",
  #'   batch_ids = "Batch",
  #'   condition = "condition",
  #'   down_sample = FALSE)
  #'   }
  #' @export
  prepare_data <- function(data_dir = NULL,
                           flowset = NULL,
                           markers = NULL,
                           pattern = "\\.fcs",
                           metadata = NULL,
                           filename_col = "filename",
                           sample_ids = NULL,
                           batch_ids = NULL,
                           condition = NULL,
                           anchor = NULL,
                           down_sample = TRUE,
                           sample_size = 500000,
                           sampling_type = "random",
                           seed = 473,
                           panel = NULL,
                           panel_channel = "fcs_colname",
                           panel_antigen = "antigen",
                           transform = TRUE,
                           cofactor = 5,
                           derand = TRUE,
                           .keep = FALSE,
                           clean_colnames = TRUE){
    
    # Stop if no data is given
    if(is.null(data_dir) & is.null(flowset)) stop("No data given.")
    
    if(!is.null(data_dir)){
      # Remove slash at end of data_dir
      if(data_dir %>% endsWith("/")) data_dir <- data_dir %>% stringr::str_sub(end = -2)
      
      # Compile directory to flowset
      if(is.null(flowset)){
        flowset <- data_dir %>%
          cyCombine::compile_fcs(pattern = pattern)
      }
      
      # Look for metadata in data_dir
      if(!is.null(metadata)){
        if(!"data.frame" %in% class(metadata)){
          if(!file.exists(file.path(metadata)) & file.exists(file.path(data_dir, metadata))) metadata <- file.path(data_dir, metadata)
        }
      }
    }
    # Convert flowset to dataframe
    fcs_data <- flowset %>%
      cyCombine::convert_flowset(metadata = metadata,
                                 filename_col = filename_col,
                                 sample_ids = sample_ids,
                                 batch_ids = batch_ids,
                                 condition = condition,
                                 anchor = anchor,
                                 down_sample = down_sample,
                                 sample_size = sample_size,
                                 sampling_type = sampling_type,
                                 seed = seed,
                                 panel = panel,
                                 panel_channel = panel_channel,
                                 panel_antigen = panel_antigen,
                                 clean_colnames = clean_colnames) %>%
      # Transform dataset with asinh
      purrr::when(transform ~ cyCombine::transform_asinh(., markers = markers,
                                                         cofactor = cofactor,
                                                         derand = derand,
                                                         .keep = .keep),
                  ~ .)
    
    message("Done!")
    return(fcs_data)
  }  
  
  # Applying the wrapper function -> still don't get what the marker is, I don't know how to fill this in
  # Define the directory containing the .fcs files
  data_dir <- "C:/Users/dhama/DESIGN_PROJECT/Data_own"
  
  # Define the markers
  # Ensure `markers` includes only the relevant marker names present in your data
  markers <- c(
    "CD103#BUV395", "CD4#BUV496", "CD8a#BUV615", "CD86#BUV735", "CD45#BUV805",
    "CD62L#BV421", "XCR1#BV510", "Ly6C#BV570", "CD161#BV605", "Ly6G#BV650",
    "CD64#BV711", "CD11c#BV750", "F4/80#BV786", "MHCII#FITC", "CD69#BB700",
    "CD3e#PE", "CD11b#PECF594", "CD19#PECy5", "PD1#PECy7", "CD63#AF647",
    "CD44#RedFluor710", "LD#APCeFluor780"
  )
  
  # Apply the wrapper function
  uncorrected <- prepare_data(
    data_dir = data_dir,  # Path to your .fcs files
    markers = markers,    # The marker names
    filename_col = "filename",  # Adjust to match your metadata
    batch_ids = "batch",  # Column name in metadata representing batch
    condition = NULL,     # Optional: Specify condition column if applicable
    down_sample = FALSE   # Set TRUE if you need down-sampling
  )
  
  # View the transformed data
  summary(uncorrected)
  
  
  
  # Normalization ----
  
  
  #' Batch-wise normalization of data
  #'
  #' This function normalizes the data in a batch-wise manner.
  #'   The purpose is to minimize the impact of batch correction when clustering the data prior to batch correction.
  #'   Three normalisation methods are implemented: Z-score, Rank, and Quantile normalization.
  #'   Z-score is recommended in cases where batches from a single study/experiment is merged.
  #'   Rank is recommend in cases where data from different studies/experiments are merged.
  #'   Quantile is not recommended.
  #'
  #' @param df tibble with expression values
  #' @param markers Markers to normalize. If NULL, markers will be found using the \code{\link{get_markers}} function.
  #' @param norm_method Normalization method. Should be either 'rank', 'scale' or 'qnorm'. Default: 'scale'
  #' @param ties.method The method to handle ties, when using rank. Default: 'average'. See ?rank for other options.
  #' @family batch
  #' @examples
  #' \dontrun{
  #' df_normed <- df %>%
  #'   normalize()
  #'   }
  #' @export
  normalize <- function(df,
                        markers = NULL,
                        norm_method = "scale",
                        ties.method = "average") {
    
    # Remove case-sensitivity
    norm_method <- norm_method %>% stringr::str_to_lower()
    ties.method <- ties.method %>% stringr::str_to_lower()
    
    # Error check
    if (norm_method == "rank" && ties.method %!in% c("average", "first", "last", "random", "max", "min")) {
      stop("When using norm_method = 'rank', please use an available ties.method (average, first, last, random, max, or min).")
    }
    
    # Messaging
    if (norm_method == "rank") {message("Ranking expression data..")
    } else if (norm_method == "scale") {message("Scaling expression data..")
    } else if (norm_method == "qnorm") {
      # message("Quantile normalizing expression data..")
      # Run quantile normalization
      df_normed <- cyCombine:::quantile_norm(df, markers = markers)
      return(df_normed)
    } else stop("Please use either 'scale', 'rank', or 'qnorm' as normalization method." )
    
    if (is.null(markers)) {
      # Get markers
      markers <- df %>%
        cyCombine::get_markers()
    }
    
    # Scale or rank at marker positions individually for every batch
    df_normed <- df %>%
      dplyr::group_by(.data$batch) %>%
      purrr::when(
        norm_method == "rank"  ~ dplyr::mutate(
          ., dplyr::across(dplyr::all_of(markers),
                           .fns = ~ {
                             if(sum(.x) == 0) stop("A marker is 0 for an entire batch. Please remove this marker.")
                             rank(.x, ties.method = ties.method) / length(.x)})),
        norm_method == "scale" ~ dplyr::mutate(., dplyr::across(dplyr::all_of(markers),
                                                                .fns = ~{
                                                                  if(sum(.x) == 0) stop("A marker is 0 for an entire batch. Please remove this marker.")
                                                                  scale(.x)}))
      ) %>%
      dplyr::ungroup()
    return(df_normed)
  }
  
  # Step 1: Define the Markers
  # Extract markers from the processed output (from the wrapper function)
  markers_wrapper <- colnames(uncorrected)[
    colnames(uncorrected) %in% markers
  ]
  
  normalized_wrapper <- normalize(
    df = uncorrected,       # Data to normalize
    markers = markers_wrapper, # Wrapper marker list
    norm_method = "scale",  # Normalization method (Z-score)
    ties.method = "average" # Default for rank normalization
  )
  
  summary(normalized_wrapper)
  
  
  
  #' Batch-wise quantile normalization per marker
  #'
  #' This function quantile normalizes the data in a batch-wise manner.
  #'   The purpose is to minimize the impact of batch correction when clustering the data prior to batch correction.
  #'
  #' @param df Dataframe with expression values
  #' @param markers Markers to correct. If NULL, markers will be found using the \code{\link{get_markers}} function.
  #' @family batch
  #' @examples
  #' \dontrun{
  #' df_qnorm <- preprocessed %>%
  #'   quantile_norm()
  #'   }
  quantile_norm <- function(df, markers = NULL) {
    message("Quantile normalizing expression data..")
    if (is.null(markers)) {
      # Get markers
      markers <- df %>%
        cyCombine::get_markers()
    }
    
    # Determine goal distributions for each marker by getting quantiles across all batches
    refq <- list()
    for (m in markers) {
      # Determine the quantiles
      refq[[m]] <- stats::quantile(unlist(df[,m]), probs=seq(0,1,length.out=5), names = F)
    }
    
    qnormed_expr <- df
    for (b in unique(df$batch)) {
      for (m in markers) {
        qx <- stats::quantile(unlist(df[df$batch == b, m]), probs=seq(0,1,length.out=5), names = F)
        spf <- stats::splinefun(x=qx, y=refq[[m]], method="monoH.FC", ties=min)
        
        # Apply the spline function to adjust quantiles
        qnormed_expr[qnormed_expr$batch == b, m] <- spf(unlist(df[df$batch==b, m]))
        
      }
    }
    
    return(qnormed_expr)
  }
  
    
  # Apply quantile normalization
  df_qnorm <- quantile_norm(df = uncorrected, markers = markers_wrapper)
  
  summary(df_qnorm)

  
  #Needed libraries
  
  library(dplyr)
  library(kohonen)
  
  # Clustering ----
  
  #' Create Self-Organizing Map
  #'
  #' The function uses the kohonen package to create a Self-Organizing Map.
  #'  It is used to segregate the cells for the batch correction to make the correction less affected
  #'  by samples with high abundances of a particular cell type.
  #'
  #' @inheritParams normalize
  #' @param seed The seed to use when creating the SOM.
  #' @param xdim The x-dimension size of the SOM.
  #' @param ydim The y-dimension size of the SOM.
  #' @param rlen Number of times the data is presented to the SOM network
  #' @family batch
  #' @examples
  #' \dontrun{
  #' labels <- uncorrected %>%
  #'   create_som()
  #'   }
  #' @export
  #' @return A vector of clustering labels
  create_som <- function(df,
                         markers = NULL,
                         seed = 473,
                         rlen = 10,
                         xdim = 8,
                         ydim = 8) {
    if (is.null(markers)) {
      # Get markers
      markers <- df %>%
        cyCombine::get_markers()
    }
    
    # SOM grid on overlapping markers, extract clustering per cell
    message("Creating SOM grid..")
    set.seed(seed)
    labels <- df %>%
      dplyr::select(markers) %>%
      as.matrix() %>%
      kohonen::som(grid = kohonen::somgrid(xdim = xdim, ydim = ydim),
                   rlen = rlen,
                   dist.fcts = "euclidean")
    
    labels <- labels$unit.classif
    
    return(labels)
  } 
  
  #I'm going to continue with the normalized_wrapper because the df_qnorm is the same as the obtained after wrapper function
  som_labels_wrapper <- create_som(
    df = normalized_wrapper,
    markers = markers
  )
  
  # Add SOM cluster labels to the normalized dataset
  normalized_wrapper$som_cluster <- som_labels_wrapper
  
  # Inspect the clustering result
  table(normalized_wrapper$som_cluster)
  
  # Visualize SOM clustering distribution
  library(ggplot2)
  
  # Create a data frame with cluster sizes
  cluster_sizes <- as.data.frame(table(normalized_wrapper$som_cluster))
  colnames(cluster_sizes) <- c("Cluster", "Size")
  
  # Plot cluster sizes
  ggplot(cluster_sizes, aes(x = as.numeric(Cluster), y = Size)) +
    geom_bar(stat = "identity", fill = "skyblue") +
    labs(title = "SOM Cluster Sizes",
         x = "Cluster ID",
         y = "Number of Cells") +
    theme_minimal()
  
  
  
  # Batch correction ----
  
  #' Correct data using ComBat
  #'
  #' Compute the batch correction on the data using the ComBat algorithm.
  #'  Define a covariate, either as a character vector or name of tibble column.
  #'  The covariate should preferable be the cell condition types but can be any column that infers heterogeneity in the data.
  #'  The function assumes that the batch information is in the "batch" column and the data contains a "sample" column with sample information.
  #'
  #' @inheritParams normalize
  #' @param label The cluster or cell type label. Either as a column name or vector.
  #' @param covar The covariate ComBat should use. Can be a vector or a column name in the input tibble.
  #'   If NULL, no covar will be used
  #' @param anchor Experimental: A column or vector specifying which samples are replicates and which are not. If specified, this column will be used as a covariate in ComBat. Be aware that it may be confounded with the condition.
  #' @param parametric Default: TRUE. If TRUE, the parametric version of ComBat is used. If FALSE, the non-parametric version is used.
  #' @param method Default: "ComBat". Choose "ComBat" for cytometry data and "ComBat_seq" for bulk RNAseq data.
  #' @param ref.batch Optional. A string of the batch that should be used as the reference for batch adjustments.
  #' @param ... Additional parameters to pass onto ComBat and ComBat_seq
  #' @family batch
  #' @examples
  #' \dontrun{
  #' corrected <- uncorrected %>%
  #'   correct_data(label = labels, covar = "condition")
  #'   }
  #' @export
  correct_data <- function(df,
                           label,
                           markers = NULL,
                           method = c("ComBat", "ComBat_seq"),
                           covar = NULL,
                           anchor = NULL,
                           ref.batch = NULL,
                           parametric = TRUE) {
    method <- match.arg(method)
    message("Batch correcting data..")
    # Check for batch column
    cyCombine:::check_colname(colnames(df), "batch", "df")
    if (is.null(markers)) {
      # Get markers
      markers <- df %>%
        cyCombine::get_markers()
    }
    
    # Add ID column to retain data order
    if("id" %!in% colnames(df)) df$id <- seq_len(nrow(df))
    
    # Add label to df
    if (length(label) == 1) {
      cyCombine:::check_colname(colnames(df), label, "df")
    } else {
      df$label <- label
      label <- "label"
    }
    
    # Add covar to df, if given
    if (!is.null(covar)) {
      if (length(covar) == 1) {
        cyCombine:::check_colname(colnames(df), covar, "df")
        df[[covar]] <- as.factor(df[[covar]])
      } else {
        # Covar was given as a vector
        df$covar <- as.factor(covar)
        covar <- "covar"
      }
      # Ensure there is more than 1 factor level
      if (nlevels(df[[covar]]) == 1) covar <- NULL
    }
    
    # Add anchor to df, if given
    if (!is.null(anchor)) {
      if (length(anchor) == 1) {
        cyCombine:::check_colname(colnames(df), anchor)
        df[[anchor]] <- as.factor(df[[anchor]])
      } else {
        # Anchor was given as a vector
        df$anchor <- as.factor(anchor)
        anchor <- "anchor"
      }
      # Ensure there is more than 1 factor level
      if (nlevels(df[[anchor]]) == 1) anchor <- NULL
    }
    
    # Determine combat method
    combat <- function(x, batch, sample, mod_matrix, parametric, ref.batch, ...) {
      x <- t(x)
      colnames(x) <- sample
      if (method == "ComBat") {
        x <- sva::ComBat(
          x,
          batch = as.character(batch),
          mod = mod_matrix,
          par.prior = parametric,
          ref.batch = ref.batch,
          prior.plots = FALSE
        )
      } else if (method == "ComBat_seq") {
        x <- sva::ComBat_seq(
          x,
          batch = as.character(batch),
          covar_mod = mod_matrix
        )
      }
      return(t(x))
    }
    
    corrected_data <- df %>%
      dplyr::group_by(.data[[label]]) %>%
      # Correct (modify) each label group with ComBat
      dplyr::group_modify(.keep = TRUE, function(df, ...) {
        # Initiate anchor and covar counter
        num_covar <- 1
        num_anchor <- 1
        # Detect if only one batch is present in the node
        num_batches <- df$batch %>%
          factor() %>%
          nlevels()
        lab <- df[[label]][1] # Current label group
        if (num_batches == 1) {
          batch <- df$batch[1]
          message(paste("Label group", lab, "only contains cells from batch", batch))
          df <- df %>% dplyr::select(-label) # Not removed from output, but removed here to prevent bug
          return(df)
        }
        message(paste("Correcting Label group", lab))
        # Calculate number of covars in the node
        if (!is.null(covar)) {
          
          # Only use covar, if it does not confound with batch
          if (!cyCombine:::check_confound(df$batch, stats::model.matrix(~df[[covar]]))) {
            num_covar <- df[[covar]] %>%
              factor() %>%
              nlevels()
            
            # If a node is heavily skewed to a single covar, it should be treated as having only 1 covar.
            # Get number of cells in the condition with most cells
            covar_counts <- df %>%
              dplyr::count(.data[[covar]]) %>%
              dplyr::pull(n)
            
            if (sum(covar_counts) < max(covar_counts) + num_covar*5) {
              message("The label group almost exclusively consists of cells from a single covar. Therefore, covar is ignored for this label group")
              num_covar <- 1
            }
          } else {
            message("Covar is confounded with batch. Ignoring covar in this label group")
          }
        }
        # Do a similar check on anchor
        if (!is.null(anchor)) {
          if (!cyCombine:::check_confound(df$batch, stats::model.matrix(~df[[anchor]]))) {
            num_anchor <- df[[anchor]] %>%
              factor() %>%
              nlevels()
            
            # If a node is heavily skewed to a single anchor, it should be treated as having only 1 covar.
            # Get number of cells in the anchor with most cells
            anchor_counts <- df %>%
              dplyr::count(.data[[anchor]]) %>%
              dplyr::pull(n)
            
            if (sum(anchor_counts) < max(anchor_counts) + num_anchor*5) {
              message("The label group almost exclusively consists of cells from a single anchor group. Therefore, anchor is ignored for this label group")
              num_anchor <- 1
            }
          } else {
            message("Anchor is confounded with batch. Ignoring anchor in this label group")
          }
        }
        if (num_covar > 1 & num_anchor > 1) {
          # If neither covar nor anchor confounds with batch but they do each other, prioritise covar
          if (cyCombine:::check_confound(df$batch, stats::model.matrix(~df[[covar]] + df[[anchor]]))) {
            num_anchor <- 1
            message("Anchor and covar are confounded. Ignoring anchor in this label group")
          }
        }
        
        # Determine model
        if (num_covar > 1 & num_anchor == 1) {
          mod_matrix <- stats::model.matrix(~ df[[covar]])
        } else if (num_covar > 1 & num_anchor > 1) {
          mod_matrix <- stats::model.matrix(~df[[covar]] + df[[anchor]])
        } else if (num_covar == 1 & num_anchor > 1) {
          mod_matrix <- stats::model.matrix(~df[[anchor]])
        } else if (num_covar == 1 & num_anchor == 1) {
          mod_matrix <- NULL # No model matrix needed
        }
        
        
        
        # Compute ComBat correction
        ComBat_output <- df %>%
          dplyr::select(dplyr::all_of(markers)) %>%
          # t() %>%
          # The as.character is to remove factor levels not present in the SOM node
          combat(
            batch = df$batch,
            sample = df$sample,
            mod_matrix = mod_matrix,
            parametric = parametric,
            ref.batch = ref.batch,
            ...) %>%
          # t() %>%
          tibble::as_tibble() %>%
          dplyr::bind_cols(
            dplyr::select(df,
                          -dplyr::all_of(c(markers, label)))) %>%
          # Cap values to range of input data
          dplyr::mutate(dplyr::across(dplyr::all_of(markers),
                                      function(x) {
                                        min <- min(df[[dplyr::cur_column()]])
                                        max <- max(df[[dplyr::cur_column()]])
                                        x <- ifelse(x < min, min, x)
                                        x <- ifelse(x > max, max, x)
                                        return(x)
                                      }))
        # Only add covar column, if it is not null
        if (!is.null(covar)) ComBat_output[[covar]] <- df[[covar]]
        # Only add anchor column, if it is not null
        if (!is.null(anchor)) ComBat_output[[anchor]] <- df[[anchor]]
        
        return(ComBat_output)
      }) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(id) %>%
      dplyr::select(id, dplyr::everything()) %>%
      dplyr::mutate(batch = as.factor(batch))
    return(corrected_data)
  }
  
  #' Alternate correction
  #'
  #' This function allows running ComBat with a custom covar mod matrix.
  #'  A model could look like \code{stats::model.matrix(~df$covar+df$label)}
  #'
  #' @inheritParams correct_data
  #' @param mod Covariate model to use in ComBat.
  #' @examples
  #' \dontrun{
  #' corrected <- uncorrected %>%
  #'   correct_data(mod = stats::model.matrix(~df$covar+df$label))
  #'   }
  correct_data_alt <- function(df,
                               mod,
                               markers = NULL,
                               parametric = TRUE) {
    message("Batch correcting data..")
    if (is.null(markers)){
      # Get markers
      markers <- df %>%
        cyCombine::get_markers()
    }
    
    
    corrected_data <- df %>%
      # Compute ComBat correction
      dplyr::select(dplyr::all_of(markers)) %>%
      t() %>%
      # The as.character is to remove factor levels not present in the SOM node
      sva::ComBat(batch = as.character(df$batch),
                  mod = mod,
                  par.prior = parametric) %>%
      t() %>%
      tibble::as_tibble() %>%
      dplyr::mutate(batch = df$batch,
                    sample = df$sample,
                    id = df$id) %>%
      # Cap values to range of input data
      dplyr::mutate(dplyr::across(dplyr::all_of(markers),
                                  function(x) {
                                    min <- min(df[[dplyr::cur_column()]])
                                    max <- max(df[[dplyr::cur_column()]])
                                    x <- ifelse(x < min, min, x)
                                    x <- ifelse(x > max, max, x)
                                    return(x)
                                  })) %>%
      dplyr::select(id, dplyr::everything())
    return(corrected_data)
  }

  # Load required packages
  library(dplyr)
  library(tibble)
  library(sva)  

  # Specify parameters
  # Use the normalized data from the SOM step
  input_data <- normalized_wrapper
  
  # Use the SOM cluster labels as the `label` for batch correction
  som_labels <- input_data$som_cluster  

  # For simplicity, this example assumes no covariates
  covar <- NULL
  
  # Apply the batch correction function
  corrected_data <- correct_data(
    df = input_data,
    label = som_labels,
    markers = markers,
    covar = covar,
    method = "ComBat",   # Choose between "ComBat" or "ComBat_seq"
    parametric = TRUE    # Use parametric ComBat
  )
  
  # Summarize the corrected data
  summary(corrected_data)  
  
  # Wrapper ----
  
  #' Run batch correction on data
  #'
  #' This is a wrapper function for the cyCombine batch correction workflow.
  #'  To run the workflow manually, type "batch_correct" to see the source code of this wrapper and follow along or read the vignettes on the GitHub page \url{https://github.com/biosurf/cyCombine}.
  #'
  #' @inheritParams create_som
  #' @inheritParams correct_data
  #' @inheritParams normalize
  #' @family batch
  #' @examples
  #' \dontrun{
  #' corrected <- uncorrected %>%
  #'   batch_correct(markers = markers,
  #'   covar = "condition")
  #'   }
  #' @export
  batch_correct <- function(df,
                            label = NULL,
                            xdim = 8,
                            ydim = 8,
                            rlen = 10,
                            parametric = TRUE,
                            method = c("ComBat", "ComBat_seq"),
                            ref.batch = NULL,
                            seed = 473,
                            covar = NULL,
                            anchor = NULL,
                            markers = NULL,
                            norm_method = "scale",
                            ties.method = "average",
                            ...) {
    # A batch column is required
    cyCombine:::check_colname(colnames(df), "batch", "df")
    if (any(is.na(df$batch))) { # Check for NAs
      message("Some batches contain NAs. These will be removed")
      warning("Some batches contain NAs. These will be removed")
      df <- df %>%
        dplyr::filter(!is.na(batch))
    }
    
    # Create SOM on scaled data
    if (is.null(label)) {
      label <- df %>%
        cyCombine::normalize(markers = markers,
                             norm_method = norm_method,
                             ties.method = ties.method) %>%
        cyCombine::create_som(markers = markers,
                              rlen = rlen,
                              seed = seed,
                              xdim = xdim,
                              ydim = ydim)
    }
    
    
    # Run batch correction
    corrected <- df %>%
      cyCombine::correct_data(label = label,
                              covar = covar,
                              anchor = anchor,
                              markers = markers,
                              parametric = parametric,
                              method = method,
                              ref.batch = ref.batch,
                              ...)
    message("Done!")
    return(corrected)
  }
  
  # Run batch correction using 'batch' as the covariate
  corrected_data <- batch_correct(
    df = normalized_wrapper,  # Start with the normalized data
    markers = markers,        # The list of marker columns
    xdim = 8,                 # SOM grid dimensions (adjust if needed)
    ydim = 8,
    rlen = 10,                # Number of iterations for SOM
    parametric = TRUE,        # Use parametric version of ComBat
    method = "ComBat",        # Choose ComBat for cytometry data
    covar = "batch",          # Use 'batch' as the covariate for correction
    norm_method = "scale"     # Apply scaling normalization before correction
  )
  
  
  # View summary of the corrected data
  summary(corrected_data)
  
  
  #EXPORT DATA
  
  #' Exporting a dataframe to SingleCellObject
  #'
  #' Conversion of dataframe into separate data and metadata objects for subsequent transformation.
  #'  Rename variable names to fit the requirements of SCE-based tools
  #'
  #' @param df Tibble with expression values and metadata
  #' @param markers Markers to include in exprs and counts object of SCE. If NULL, markers will be found using the \code{\link{get_markers}} function.
  #' @param non_markers Non-markers to include as colData in SCE. If NULL, non_markers will be based on cyCombine::non_markers.
  #' @param sample_col It is the column name in the df that contains the sample names. Defaults to 'sample'.
  #' @param panel Optional: Panel as a data.frame. Should have colnames Channel, Marker, Type unless otherwise specified in the panel_ args. Should be included if you want to store FCS files
  #' @param panel_channel Optional: Only used if panel is given. It is the column name in the panel data that contains the channel names
  #' @param panel_antigen Optional: Only used if panel is given. It is the column name in the panel data that contains the antigen names
  #' @param panel_type Optional: Only used if panel is given. It is the column name in the panel data that contains the antigen types (none, state, type).
  #'  "none" will be excluded from SCE. Set to NULL to disregard.
  #' @param transform_cofactor The cofactor to use when reverse-transforming to raw counts
  #' @family export
  #' @examples
  #' \dontrun{
  #' sce <- df %>%
  #'   df2SCE(markers = markers, non_markers = NULL, panel = panel)
  #'   }
  #' @export
  df2SCE <- function(
    df,
    markers = NULL, non_markers = NULL,
    scatter = NULL,
    clean_names = FALSE,
    sample_col = "sample",
    panel = NULL, panel_channel = "Channel",
    panel_antigen = "Marker", panel_type = NULL,
    transform_cofactor = 5) {
    
    # Check for packages
    cyCombine:::missing_package("SingleCellExperiment", "Bioc")
    
    message("Converting dataframe to SingleCellExperiment object...")
    
    # Get the non markers and prepare column data
    if (is.null(non_markers)) {
      # Get non markers
      if (sample_col == "sample") {
        non_markers <- cyCombine::non_markers
      } else {
        non_markers <- c(cyCombine::non_markers, dplyr::all_of(sample_col))
      }
    }
    
    # Get the column data (from non markers)
    if (sample_col %in% colnames(df)) {
      colData <- df %>%
        dplyr::select(dplyr::any_of(non_markers))
      
      if (sample_col != 'sample_id') {
        colData <- colData %>%
          dplyr::rename(sample_id = sample_col)
      }
      
    } else {
      stop("Error, none of the non_markers/sample_col are available in the dataframe. You cannot make an SCE without sample names.")
    }
    
    # Get markers and check
    if (is.null(markers)) {
      # Get markers
      markers <- df %>%
        cyCombine::get_markers()
    }
    
    sapply(markers, function(x) {
      cyCombine:::check_colname(colnames(df), x, "df")})
    
    # Extract expression data and transpose to fit SCE format
    exprs <- t(dplyr::select(df, dplyr::all_of(c(scatter, markers))))
    
    
    # Prepare the experiment info table - first identify true meta data columns
    colCounting <- colData %>%
      dplyr::group_split(sample_id) %>%
      lapply(function(x) {(apply(x, 2, function(y) {length(table(y)) == 1}))})
    
    col_stability <- lapply(colCounting, function(z) {names(which(z))}) %>% unlist() %>% table()
    stable_cols <- names(which(col_stability == length(colCounting)))
    
    # Swap ordering
    stable_cols <- c("sample_id", stable_cols[stable_cols != "sample_id"])
    
    experiment_info <- colData %>%
      dplyr::select(dplyr::all_of(stable_cols)) %>%
      dplyr::group_by(sample_id) %>%
      dplyr::mutate(n_cells = dplyr::n()) %>%
      dplyr::distinct(sample_id, .keep_all = TRUE) %>%
      as.data.frame()
    
    
    
    # Prepare row data if available
    if (!is.null(panel) && is.data.frame(panel)) {
      # Check presence of panel's data names
      sapply(c(panel_channel, panel_antigen, panel_type), function(x) {
        cyCombine:::check_colname(colnames(panel), x, "panel")})
      
      rowData <- panel[match(rownames(exprs), panel[[panel_antigen]]), ]
      rm(panel)
      # Exclude none's
      if (!is(panel_type, "NULL")) {
        if ("none" %in% rowData[[panel_type]]){
          rowData <- rowData %>%
            dplyr::filter(.data[[panel_type]] != "none")
        }
      }
      
      
      # Change colnames to fit SCE standard
      if (panel_channel != "channel_name") {
        rowData <- rowData %>%
          dplyr::rename(channel_name = panel_channel)
      }
      if (panel_antigen != "marker_name") {
        rowData <- rowData %>%
          dplyr::rename(marker_name = panel_antigen)
      }
      if (panel_type != "marker_class" && !is(panel_type, "NULL")) {
        rowData <- rowData %>%
          dplyr::rename(marker_class = panel_type)
      }
      
      # Change marker names to exclude spaces and dashes
      if (clean_names) {
        rowData[["marker_name"]] <- rowData[["marker_name"]] %>%
          stringr::str_remove_all("^\\d+[A-Za-z]+_") %>%
          stringr::str_remove_all("[ _-]")
      }
      
      
      # Subset rowData to rows in exprs
      rowData <- rowData %>%
        dplyr::filter(.data[["marker_name"]] %in% rownames(exprs))
      
    } else {
      rowData <- NULL
      warning("To store as FCS files later, you should include panel information at this step.")
    }
    
    exprs <- exprs[c(scatter, markers), ]
    # Creating the SCE
    sce <- SingleCellExperiment::SingleCellExperiment(
      list(exprs = exprs,
           counts = exprs %>%
             t() %>%
             tibble::as_tibble() %>%
             transform_asinh(reverse = TRUE,
                             cofactor = transform_cofactor,
                             derand = FALSE,
                             .keep = TRUE,
                             markers = markers) %>%
             t()),
      colData = colData,
      rowData = rowData,
      metadata = list("experiment_info" = experiment_info)
    )
    message("Your SingleCellExperiment object is now created. The 'counts' assay contains reverse transformed expression values and 'exprs' contains expression values.")
    
    return(sce)
  }

  
  BiocManager::install('SingleCellExperiment')
  
  # Use marker_wrapper directly
  sce <- df2SCE(
    df = corrected_data,
    markers = markers_wrapper,  # Use the markers identified earlier
    sample_col = "sample",    # Sample column name
    transform_cofactor = 5    # Ensure the cofactor matches
  )
  
  # Check the SCE object
  sce
  
  
  #' Convert SingelCellExperiment into flowSet and store as FCS files
  #'
  #' Wrapper function that makes it easier to go from a SCE to flowSet and written FCS files.
  #'   The function uses CATALYST::sce2fcs and flowCore::write.flowSet to store the FCS files.
  #'
  #' @inheritParams CATALYST::sce2fcs
  #' @inheritParams flowCore::write.flowSet
  #' @param sce SingleCellExperiment to write to FCS files
  #' @param outdir If given, the flowSet will be stored in FCS files
  #' @param randomize (Default: FALSE) Logical determining whether counts are randomized for plotting purposes or not. Only works when assay = "counts" (default).
  #' @family export
  #' @examples
  #' \dontrun{
  #'  sce2FCS(sce, outdir = "fcs_files")
  #'   }
  #' @export
  sce2FCS <- function(sce,
                      outdir = NULL,
                      split_by = "sample_id",
                      assay = "counts",
                      keep_dr = TRUE,
                      keep_cd = TRUE,
                      randomize = FALSE) {
    
    
    # Check CATALYST is installed
    cyCombine:::missing_package("CATALYST", repo = "Bioc")
    
    # If no panel information is available
    stopifnot("Your SingleCellExperiment should contain channel information to be stored.\n
            Consider rerunning df2sce with a panel to continue." = !is.null(CATALYST::channels(sce)))
    
    # Randomize counts
    if (randomize && assay == "counts") {
      SummarizedExperiment::assay(sce, "counts") <- SummarizedExperiment::assay(sce, "counts") %>%
        cyCombine:::randomize_matrix()
    } else if (randomize && assay != "counts") {
      warning("Please only use randomization with count values.")
    }
    
    # Convert to flowset
    fcs <- CATALYST::sce2fcs(sce,
                             keep_dr = keep_dr,
                             keep_cd = keep_cd,
                             split_by = split_by,
                             assay = assay)
    
    
    # Write FCS files
    if (!is.null(outdir)) {
      # Create output directory
      cyCombine:::check_make_dir(outdir)
      
      message("Writing FCS files to ", outdir)
      switch(class(fcs)[1],
             "flowFrame" = flowCore::write.FCS(fcs, filename = paste0(outdir, sce$sample_id[1], ".fcs")),
             "flowSet" = flowCore::write.flowSet(fcs, outdir = outdir))
    }
    
    return(fcs)
  }

  BiocManager::install(c("CATALYST", "flowCore"))
  
  fcs <- sce2FCS(
    sce = sce,
    outdir = "fcs_output",      # Directory where FCS files will be saved
    split_by = "sample_id",     # Split FCS files by sample
    assay = "counts",           # Use counts for the export
    keep_dr = TRUE,             # Keep reduced dimensions
    keep_cd = TRUE,             # Keep column data
    randomize = FALSE           # Do not randomize
  )
  