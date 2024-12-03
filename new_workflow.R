# Load packages
library(cyCombine)
library(tidyverse)
library(flowCore)
library(dplyr)
library(tibble)

# Directory with FCS files
data_dir <- "C:/Users/dhama/DESIGN_PROJECT/Data_own"


### EXTRACT PANEL AND METADATA FILE FROM FCS

# List FCS files
fcs_files <- list.files(data_dir, pattern = "\\.fcs$", full.names = TRUE)

# Initialize panel and metadata tibbles
panel <- tibble(Channel = character(), Marker = character(), Type = character())
metadata <- tibble(Filename = character(), batch = character(), condition = character())

# Iterate over FCS files to extract information
for (fcs_path in fcs_files) {
  # Read FCS file
  fcs <- read.FCS(fcs_path, transformation = FALSE)
  
  # Extract marker and channel names
  params <- pData(parameters(fcs))
  channels <- params$name
  markers <- params$desc
  
  # Add to panel (assuming all FCS files share the same markers and channels)
  if (nrow(panel) == 0) {
    panel <- tibble(Channel = channels, Marker = ifelse(is.na(markers), channels, markers), Type = "state")
  }
  
  # Extract metadata from filenames
  metadata <- tibble(Filename = basename(fcs_files)) %>%
    mutate(
      # Extract 'DayX' as batch information
      batch = str_extract(Filename, "Day\\d+"),
      
      # Extract 'KO' or 'WT' as condition
      condition = str_extract(Filename, "(KO|WT)")
    )
  
}


# Save panel and metadata as CSV for future use
write.csv(panel, file.path(data_dir, "panel.csv"), row.names = FALSE)
write.csv(metadata, file.path(data_dir, "metadata.csv"), row.names = FALSE)

# Output panel and metadata summary
print("Generated Panel:")
print(panel)

print("Generated Metadata:")
print(metadata)


### CONTINUATION OF PROVIDED WORKFLOW

# Extract markers from panel
panel_file <- file.path(data_dir, "panel.csv") # Can also be .xlsx
raw_panel <- read.csv(panel_file)
metadata_file <- file.path(data_dir, "metadata.csv") # Can also be .xlsx


# Extract markers of interest
markers <- raw_panel %>%
  dplyr::filter(Type != "none") %>%   #I specified library because otherwise it gave me an error
  dplyr::pull(Marker)

print(markers)

# with this we have obtained some things in out vector that are not markers, e.g. Time, Original_ID, FSC-A...
# I am excluding them manually
markers <- markers[!markers %in% c("FSC-A", "FSC-H", "FSC-W", "SSC-A", "SSC-H", "SSC-W", "Time", "Original_ID")]
print(markers)

# Now we map markers to channels
marker_channel_map <- panel %>%
  dplyr::filter(Marker %in% markers) %>%
  dplyr::select(Marker, Channel)

# Replace Marker names with Channel names
mapped_markers <- marker_channel_map %>%
  dplyr::pull(Channel)

# Validate that all mapped markers are in the FCS data
test_file <- fcs_files[1]
fcs <- read.FCS(test_file, transformation = FALSE)
exprs_data <- as.data.frame(exprs(fcs))

missing_channels <- setdiff(mapped_markers, colnames(exprs_data))
if (length(missing_channels) > 0) {
  print("Mapped markers missing in FCS data:")
  print(missing_channels)
} else {
  print("All mapped markers align with FCS data!")
}

# Filter the markers to only those present in the FCS files
valid_markers <- intersect(mapped_markers, colnames(exprs_data))
print("Valid markers to use:")
print(valid_markers)

# Prepare a tibble from directory of FCS files
uncorrected <- prepare_data(
  data_dir = data_dir,
  metadata = metadata_file, 
  filename_col = "Filename",
  batch_ids = "batch",
  condition = "condition",
  down_sample = FALSE,
  markers = valid_markers
)
print(str(uncorrected))

# Store result
saveRDS(uncorrected, file = file.path(data_dir, "uncorrected.RDS"))