# Loading necessary packages
library(ggplot2)
library(dplyr)
library(stringr)
library(readr)
library(reshape2)
library(pracma)
library(forcats)
library(gtools)
library(purrr)



#' Sorts the peaks detected in a given numeric vector.
#'
#' This function sorts the peaks detected in a given numeric vector.
#'
#' @param x A numeric vector containing the data in which peaks are to be detected.
#' @param ... Additional arguments to be passed to the \code{\link{findpeaks}} function.
#'
#' @return A matrix containing the peaks sorted in ascending order of their index in the input data. Column 1 is peak max intensity, Column 2 is peak max index, Column 3 is peak start index, Column 4 is peak end index.
#'
#' @examples
#' sorted_peaks <- findpeaks.sort(data_vector, minpeakheight = 2, minpeakdistance = 40, npeaks = 7)
#'
#' @importFrom pracma findpeaks
#' @export
findpeaks.sort <- function(x, ...){
  peaks <- pracma::findpeaks(x, ...)
  peaks <- peaks[order(peaks[,2]), ]
  return(peaks)
}



#' Check if peaks exist within a specified range in a matrix column.
#'
#' This function checks if there are peaks within a specified range in a given matrix column.
#'
#' @param matrix_data A matrix where each column represents a dataset, and the target column is to be checked for peaks.
#' @param lower_bound The lower bound of the range to check for peaks.
#' @param upper_bound The upper bound of the range to check for peaks.
#'
#' @return A logical value indicating whether peaks exist within the specified range in the target column.
#'
#' @examples
#' # Create a matrix data
#' data_matrix <- matrix(c(1, 2, 3, 4, 5, 6), ncol = 2)
#' # Check for peaks in the second column within the range [2, 4]
#' peak_exists <- find_peakl(data_matrix, 2, 4)
#' peak_exists
#'
#' @export
find_peakl <- function(matrix_data, lower_bound, upper_bound) {
  row_data <- matrix_data[, 2]
  # Check if any element in the column satisfies the condition
  if (any(row_data >= lower_bound & row_data <= upper_bound)) {
    peak_exist <- TRUE
  } else {
    peak_exist <- FALSE
  }
  return(peak_exist)
}



#' Find indices of rows meeting a specified condition in a matrix column.
#'
#' This function finds the indices of rows in a matrix column that meet a specified condition.
#'
#' @param matrix_data A matrix where each column represents a dataset, and the target column is to be checked for the condition.
#' @param lower_bound The lower bound of the range to check for the condition.
#' @param upper_bound The upper bound of the range to check for the condition.
#'
#' @return A vector containing the indices of rows where the values in the target column satisfy the condition.
#'
#' @examples
#' # Create a matrix data
#' data_matrix <- matrix(c(1, 2, 3, 4, 5, 6), ncol = 2)
#' # Find indices of rows in the second column where values are within the range [2, 4]
#' selected_indices <- find_peaki(data_matrix, 2, 4)
#' selected_indices
#'
#' @export
find_peaki <- function(matrix_data, lower_bound, upper_bound) {
  row_data <- matrix_data[, 2]
  # Check if any element in the column satisfies the condition
  if (any(row_data >= lower_bound & row_data <= upper_bound)) {
    # Find the indices of the rows meeting the condition
    selected_indices <- which(row_data >= lower_bound & row_data <= upper_bound)
  }
  return(selected_indices)
}



#' Fit ladder data and analyze peaks.
#'
#' This function fits ladder data obtained from a directory, detects peaks, and performs a polynomial fit.
#'
#' @param directory The directory containing ladder data file.
#' @param minpeakheight The minimum peak height for peak detection (default = 2).
#' @param minpeakdistance The minimum peak distance for peak detection (default = 40).
#' @param npeaks The number of peaks to detect (default = 7).
#' @param diagnostics If TRUE, print the summary of the linear model fit to the console (default = TRUE).
#' @param plot_ladder If TRUE, plot the ladder data (default = TRUE).
#' @param plot_fit If TRUE, plot the ladder fit with detected peaks (default = TRUE).
#'
#' @return A numeric vector containing the coefficients of the fitted polynomial.
#'
#' @examples
#' # Fit ladder data and analyze peaks
#' coefficients <- ba.fitladder("path/to/directory")
#' coefficients
#'
#' @importFrom dplyr mutate
#' @importFrom ggplot2 ggplot geom_line geom_point stat_function
#' @importFrom stats lm
#' @importFrom utils list.files read.csv
#' @importFrom grDevices rgb
#' @export
ba.fitladder <- function(directory, minpeakheight = 2, minpeakdistance = 40, npeaks = 7, diagnostics = TRUE, plot_ladder = TRUE, plot_fit = TRUE) {
lad_file <- list.files(directory, pattern = 'Ladder.*\\.csv') # Find the ladder file
lad <- read.csv(paste0(directory, '/', lad_file), skip = 17, nrows = 1060) # Read it in
lad <- lad %>% mutate(across(c(Time, Value), as.numeric)) # Convert to numeric
lad_peaks <- findpeaks.sort(lad[,2], minpeakheight = minpeakheight, minpeakdistance = minpeakdistance, npeaks = npeaks) # Call peaks
lad_peaks_df <- data.frame(size_nt = c(25, 200, 500, 1000, 2000, 4000, 6000),
                           time = lad[,1][lad_peaks[,2]],
                           max_fluo = lad_peaks[,1])
lad_fit <- lm(lad_peaks_df$size_nt~poly(lad_peaks_df$time, degree = 2, raw = TRUE ))
if(diagnostics) { # Print the lm() summary to the console
  print(summary(lad_fit))
}
if(plot_ladder) { # plot the ladder (useful if you are getting less than 7 peaks, you may have to lower your threshold)
  print(lad %>%
          ggplot(aes(x = Time, y = Value))+
          geom_line())
}
if(plot_fit) { # plot the datapoints that were called as peaks against their time and the ladder fit
  print(lad_peaks_df %>%
          ggplot(aes(x = time, y = size_nt))+
          geom_point()+
          stat_function(
            fun = function(x) lad_fit$coefficients[3] * x^2 + lad_fit$coefficients[2] * x + lad_fit$coefficients[1],
            geom = "line"))
  
}
return(lad_fit$coefficients) # I only really need to save the coefficients out of this
}



#' Write traces data to CSV files.
#'
#' This function writes the traces data to CSV files, including normalized data.
#'
#' @param df A dataframe containing the columns "Time" and 'size_nt' from previous functions.
#' @param directory The directory where the original data is located, used for writing experiment info.
#' @param output_path The directory path where the output CSV files will be written (does not have to exist).
#'
#' @importFrom dplyr bind_cols select
#' @importFrom purrr map_dbl
#' @importFrom readr write_csv
#' @importFrom stringr str_split_i
#' @export
ba.write_traces <- function(df, directory, output_path) {
  # This dataframe should have the columns "Time" and 'size_nt' from above functions
  norm_df <- bind_cols(df[c('Time', 'size_nt')], as.data.frame(apply(select(df, -Time, -size_nt), 2, function(x) x/max(x))))
  # Create the output directory if it doesn't exist
  if (!dir.exists(output_path)) {
    dir.create(output_path, recursive = TRUE)
  }
  # The directory in this case is the directory that we originally started in. That way we can write experiment info in this file
  write_csv(df, paste0(output_path, '/', paste(str_split_i(directory, '/', -1), 
                                               str_split_i(directory, '/', -2), 
                                               'traces_trim', sep = '_'),
                       ".csv"))
  write_csv(norm_df, paste0(output_path,'/', paste(str_split_i(directory, '/', -1), 
                                                   str_split_i(directory, '/', -2), 
                                                   'traces_trim_norm', sep = '_'),
                            ".csv"))
}



#' Perform calibration and trim CSV files.
#'
#' This function performs calibration and trimming on CSV files containing sample traces.
#'
#' @param directory The directory containing the CSV files.
#' @param t_start The start time for trimming the data (default = 24).
#' @param output_path The directory path where the output CSV files will be written (does not need to exist).
#' @param minpeakheight The minimum peak height for peak detection (default = 2).
#' @param minpeakdistance The minimum peak distance for peak detection (default = 40).
#' @param npeaks The number of peaks to detect (default = 7).
#' @param diagnostics If TRUE, print the summary of the linear model fit to the console (default = FALSE).
#' @param plot_ladder If TRUE, plot the ladder data (default = FALSE).
#' @param plot_fit If TRUE, plot the ladder fit with detected peaks (default = FALSE).
#'
#' @return A dataframe containing normalized trimmed trace data.
#'
#' @examples
#' # Perform calibration and trim CSV files
#' calibrated_data <- ba.cal_trim_csv("path/to/directory", output_path = "output/")
#' calibrated_data
#'
#' @importFrom dplyr inner_join reduce select bind_cols
#' @importFrom readr read_csv
#' @export
ba.cal_trim_csv <- function(directory, t_start = 24, output_path,
                            minpeakheight = 2, minpeakdistance = 40, npeaks = 7, 
                            diagnostics = FALSE, plot_ladder = FALSE, plot_fit = FALSE){
  
  ### Processing metadata file first
  info_file <- list.files(path = directory, pattern = 'info.*\\.csv') # Finding the sample info file
  info_df <- read.csv(paste0(directory, '/', info_file)) # Reading in the sample info file
  colnames(info_df)[1:5] <- c('lane', 'sample_name', 'group', 'target_gene', 'experiment') # Making sure the column names on this info file are correct
  info_df$sample_name <- gsub(pattern = "[^[:alnum:]_]", "_", replacement = '_', x = info_df$sample_name) # Removing special characters
  info_df$sample_name <- ifelse(grepl(pattern = "^[0-9]", x = info_df$sample_name), paste0("s_", info_df$sample_name), info_df$sample_name) # Adding a alpha character if it starts with a number 
  info_df$group <- gsub(pattern = "[^[:alnum:]_]", replacement = '_', x = info_df$group) # Removing special characters
  info_df$group <- ifelse(grepl(pattern = "^[0-9]", x = info_df$group), paste0("g_", info_df$group), info_df$group) # Adding a alpha character if it starts with a number 
  info_df$target_gene <- gsub(pattern = "[^[:alnum:]_]", replacement = '_', x = info_df$target_gene) # Removing special characters
  info_df$target_gene <- ifelse(grepl(pattern = "^[0-9]", x = info_df$target_gene), paste0("t_", info_df$target_gene), info_df$target_gene) # Adding a alpha character if it starts with a number 
  info_df$experiment <- gsub(pattern = "[^[:alnum:]_]", replacement = '_', x = info_df$experiment) # Removing special characters
  info_df$experiment <- ifelse(grepl(pattern = "^[0-9]", x = info_df$experiment), paste0("t_", info_df$experiment), info_df$experiment) # Adding a alpha character if it starts with a number 
  if(!is.numeric(info_df$lane)){
    info_df$lane <- as.numeric(info_df$lane) # Making sure the lane is numeric so that I can use order to sort the table (if the lanes are not in order of the sample names)
  } 
  info_df <- info_df[order(info_df$lane),] # Sorting the sample table by lane so that this all can reference back to the order that the csvs are read in
  
  ### Now dealing with the sample files
  s_files <- mixedsort(list.files(path = directory, pattern = 'Sample[0-9].*\\.csv')) # Defining the sample files, reading them in numeric order
  one_sample <- NULL # Setting a variable to use to detect if there is one or more than one sample of interest
  df_list <- list() # Making a list to use to store dataframes if there is more than one sample of interest
  for (i in s_files) { # Looping through these sample files
    sample_number <- str_split_i(i, '_Sample', 2) %>% str_split_i(., '\\.', 1) # Getting the sample number out of the file name
    if (sample_number %in% as.character(info_df$lane)) { # Testing if it is a sample of interest for us
      if (length(info_df$lane) == 1) { # If there is only one sample of interest in our metadata table, this should only be true once
        one_sample <- i # Setting the one_sample variable to the file of the one sample we are interested in
      } else { #If there is more than one sample of interest
        tmp_df <- read.csv(paste0(directory, '/', i), skip = 17, nrows = 1060) # Read it in. You may have do adjust the "skip = 17" based on the exact export format from the instrument
        df_list[[length(df_list) + 1]] <- tmp_df # Add it to the list of sample dataframes
      }
    }
  }
  
  if (is.null(one_sample)) { # If the 'one_sample' variable is still NULL, there is more than one sample
    trace_df <- df_list %>% reduce(inner_join, by = "Time") # combining all dataframes into a big dataframe with a single time column
  } else #If there is a single sample of interest, the 'one_sample' variable should now carry this file name
    trace_df <- read.csv(paste0(directory, '/', one_sample), skip = 17, nrows = 1060)
  
  # Now processing the dataframes with the traces
  colnames(trace_df)[-1] <- info_df$sample_name # Adding the names as the column names
  trim_trace_df <- trace_df[trace_df$Time >= t_start, ] # Trimming the dataframe
  fit_coef <- ba.fitladder(directory = directory, 
                           minpeakheight = minpeakheight, minpeakdistance = minpeakdistance, npeaks = npeaks, 
                           diagnostics = diagnostics, plot_ladder = plot_ladder, plot_fit = plot_fit)
  trim_trace_df$size_nt <- fit_coef[1]+ fit_coef[2]*trim_trace_df$Time + fit_coef[3] * (trim_trace_df$Time)^2
  trim_trace_df <- trim_trace_df %>% select(Time, size_nt, everything())
  # Create the output directory if it doesn't exist
  if (!dir.exists(output_path)) {
    dir.create(output_path, recursive = TRUE)
  }
  ba.write_traces(trim_trace_df, directory = directory, output_path = output_path)
  trim_trace_df_norm <- bind_cols(trim_trace_df[c('Time', 'size_nt')], as.data.frame(apply(select(trim_trace_df, -Time, -size_nt), 2, function(x) x/max(x))))
  return(trim_trace_df_norm)
}

#' Analyze 28S rRNA integrity in zebrafish total RNA.
#'
#' This function analyzes peaks in rRNA samples from the provided input data frame. It is time-based to keep the number of datapoints fixed for each sample.
#' The indices used to define each peak position will depend on how the raw traces are trimmed. All default parameters are interdependent on those of previous functions (i.e. assume using t_start = 24 in the "ba.cal_trim_csv()" function).
#' With consistent data, the default indices used to define the rRNA peaks can be employed across many chips. 
#' With more variability between chips, it is advantageous to set custom indices for each chip
#'
#' @param input_df Either a file path to a CSV or an R dataframe containing the data to analyze.
#' @param t_start The start time for trimming the data (default = 35).
#' @param t_end The end time for trimming the data (default = 55).
#' @param minpeakheight The minimum peak height for peak detection (default = 0.15).
#' @param minpeakdistance The minimum peak distance for peak detection (default = 40).
#' @param npeaks The number of peaks to detect (default = 4).
#' @param range_28s The range of indices searched for the 28S rRNA peak (default = c(460,500)).
#' @param range_18s The range of indices searched for the 18S rRNA peak (default = c(330,375)).
#' @param range_deg_l The range of indices searched for the high molecular weight 28S rRNA degradation product peak (default = c(405,450)).
#' @param range_deg_s The range of indices searched for the low molecular weight 28S rRNA degradation product peak (default = c(265,300)).
#'
#' @return A dataframe containing the analyzed rRNA peaks.
#'
#' @examples
#' # Analyze 28S rRNA integrity in zebrafish total RNA from a single chip
#' analyzed_peaks <- ba.analyze_rrna("path/to/input.csv")
#' analyzed_peaks
#' 
#' # Analyze 28S rRNA integrity in zebrafish total RNA from multiple chips with varied profiles
#' ## Finding the trimmed, normalized sample files
#' input_files <- list.files(path = 'processed_data/tables',pattern = 'traces_trim_norm.csv', recursive = TRUE, full.names = TRUE)
#' 
#' ## Employing "ba.analyze_rrna_list_only()" to see where the peaks are in each file
#' select_range_list <- list()
#' for (i in input_files){
#' select_range_list[[i]] <- ba.analyze_rrna_list_only(i, npeaks = 6, minpeakheight = 0.10)
#' }
#' select_range_list
#' 
#' ## Making a list based on ranges for each chip
#' m_ranges <- list(
#' list(r28s = c(505,520), r18s = c(380,395), rdl = c(465,480), rds = c(295,310)), # chip 1
#' list(r28s = c(490,505), r18s = c(365,380), rdl = c(450,465), rds = c(275,290)), # chip 2
#' list(r28s = c(505,520), r18s = c(380,395), rdl = c(465,480), rds = c(295,310)), # chip 3
#' list(r28s = c(490,505), r18s = c(365,380), rdl = c(450,465), rds = c(275,290)) # chip 4
#' )
#' 
#' ## Merging the list with the name of their input files
#' search_list <- setNames(m_ranges, input_files)

#' ## Making a list to hold the output dataframes
#' output_list <- list()
#' ## Making a for loop to do the analysis
#' for (file in input_files){
#' fhead <- str_split_i(file, '/', 3)
#' print(fhead)
#' r_sublist <- search_list[[file]]
#' rrna <- ba.analyze_rrna(file, npeaks = 6, minpeakheight = 0.10, range_28s = r_sublist$r28s, range_18s = r_sublist$r18s, range_deg_l = r_sublist$rdl, range_deg_s = r_sublist$rds) # Here I lowered the threshold a little to detect real peaks. Guessing these samples are a bit more concentrated as they have a flatter baseline
#' rrna$chip <- rep(fhead, length(rrna$sample_name))
#' output_list[[fhead]] <- rrna
#' }
#' ## collapsing this list
#' output_df <- do.call(bind_rows, output_list)
#' output_df
#' 
#'
#' @importFrom stats is.character is.data.frame
#' @importFrom dplyr bind_cols rownames
#' @export
ba.analyze_rrna <- function(input_df, t_start = 35, t_end = 55, minpeakheight = 0.15, minpeakdistance = 40, npeaks = 4,
                            range_28s = c(460,500), range_18s = c(330,375), range_deg_l = c(405,450), range_deg_s = c(265,300)){
  if (is.character(input_df)) {  # Check if input is a file path
    if (file.exists(input_df)) {
      df <- read.csv(input_df)
    } else {
      stop("File not found.")
    }
  } else if (is.data.frame(input_df)) {  # Check if input is a dataframe
    df <- input_df
    rownames(df) <- NULL # Reset the index so that this is the original reference to which to re-index
  } else {
    stop("Input must be either a file path to a CSV or a dataframe.")
  }
  if (ncol(df) > 3) { # Checking if there is more than one sample
    peaks_list <- apply(df[df$Time >= t_start & df$Time <= t_end,3:ncol(df)], 2, 
                        function(x) findpeaks.sort(x, minpeakheight = minpeakheight, minpeakdistance = minpeakdistance, npeaks = npeaks), 
                        simplify = FALSE)
  } else { # If there is only one sample, it should have only 3 columns
    peaks_list <- list()
    peaks_list[[colnames(df)[3]]] <- findpeaks.sort(df[df$Time >= t_start & df$Time <= t_end,3], minpeakheight = minpeakheight, minpeakdistance = minpeakdistance, npeaks = npeaks)
  }
  peaks_list_ri <- lapply(peaks_list, function(m) cbind(m[,1], m[,2:4] + (as.numeric(rownames(df[df$Time >= t_start & df$Time <= t_end, 1:3])[1]) -1)))
  
  # 28S peak first as the model
  ana_28s_list <- list()
  for(i in seq_along(names(peaks_list_ri))){
    n = names(peaks_list_ri[i]) # Getting the sample name
    m = peaks_list_ri[[n]] # reading in the peak matrix
    if(find_peakl(m, lower_bound = range_28s[1], upper_bound = range_28s[2])){
      rnum = find_peaki(m, lower_bound = range_28s[1], upper_bound = range_28s[2])
      ind_max = m[rnum,2] # getting the peak index
      ind_start = m[rnum,3] # getting the index of the start of the peak
      ind_end = m[rnum,4] # getting the index of the end of the peak
      p_max = df[ind_max, n] # Getting the peak intensity at max
      p_max_size = df[ind_max, 'size_nt'] # Getting the peak size at max intensity
      p_start = df[ind_start, n] # Getting the peak intensity at start
      p_start_size = df[ind_start, 'size_nt'] # Getting the peak size at start 
      p_end = df[ind_end, n] # Getting the peak intensity at start
      p_end_size = df[ind_end, 'size_nt'] # Getting the peak size at start 
      p_halfmax = p_max/2 # Getting the intensity of the peak at half max
      # Now I need to make a line that calculates the width (in bp) from this halfmax measurement
      p_df = df[ind_start:ind_end, c('size_nt', n)]
      p_df_half = p_df[p_df[,2] >= p_halfmax,]
      p_width_halfmax = p_df_half[nrow(p_df_half),1] - p_df_half[1,1]
    } else {
      p_max = mean(df[range_28s[1]:range_28s[2], n])
      p_max_size = NA
      p_start = NA
      p_start_size = NA
      p_end = NA
      p_end_size = NA 
      p_halfmax = NA
      p_width_halfmax = 0}
    results = c(n, p_max, p_max_size, p_start, p_start_size, 
                p_end, p_end_size, p_halfmax, p_width_halfmax, "rrna_28S")
    ana_28s_list[[i]] <- results
  }
  
  # 18S peak as the second
  ana_18s_list <- list()
  for(i in seq_along(names(peaks_list_ri))){
    n = names(peaks_list_ri[i]) # Getting the sample name
    m = peaks_list_ri[[n]] # reading in the peak matrix
    if(find_peakl(m, lower_bound = range_18s[1], upper_bound = range_18s[2])){
      rnum = find_peaki(m, lower_bound = range_18s[1], upper_bound = range_18s[2])
      ind_max = m[rnum,2] # getting the peak index
      ind_start = m[rnum,3] # getting the index of the start of the peak
      ind_end = m[rnum,4] # getting the index of the end of the peak
      p_max = df[ind_max, n] # Getting the peak intensity at max
      p_max_size = df[ind_max, 'size_nt'] # Getting the peak size at max intensity
      p_start = df[ind_start, n] # Getting the peak intensity at start
      p_start_size = df[ind_start, 'size_nt'] # Getting the peak size at start 
      p_end = df[ind_end, n] # Getting the peak intensity at start
      p_end_size = df[ind_end, 'size_nt'] # Getting the peak size at start 
      p_halfmax = p_max/2 # Getting the intensity of the peak at half max
      # Now I need to make a line that calculates the width (in bp) from this halfmax measurement
      p_df = df[ind_start:ind_end, c('size_nt', n)]
      p_df_half = p_df[p_df[,2] >= p_halfmax,]
      p_width_halfmax = p_df_half[nrow(p_df_half),1] - p_df_half[1,1]
    } else {
      p_max = mean(df[range_18s[1]:range_18s[2], n])
      p_max_size = NA
      p_start = NA
      p_start_size = NA
      p_end = NA
      p_end_size = NA 
      p_halfmax = NA
      p_width_halfmax = 0}
    results = c(n, p_max, p_max_size, p_start, p_start_size, 
                p_end, p_end_size, p_halfmax, p_width_halfmax, "rrna_18S")
    ana_18s_list[[i]] <- results
  }
  
  # The larger fragment of 28s degradation as the third peak
  ana_28s_deg_long_list <- list()
  for(i in seq_along(names(peaks_list_ri))){
    n = names(peaks_list_ri[i]) # Getting the sample name
    m = peaks_list_ri[[n]] # reading in the peak matrix
    if(find_peakl(m, lower_bound = range_deg_l[1], upper_bound = range_deg_l[2])){
      rnum = find_peaki(m, lower_bound = range_deg_l[1], upper_bound = range_deg_l[2])
      ind_max = m[rnum,2] # getting the peak index
      ind_start = m[rnum,3] # getting the index of the start of the peak
      ind_end = m[rnum,4] # getting the index of the end of the peak
      p_max = df[ind_max, n] # Getting the peak intensity at max
      p_max_size = df[ind_max, 'size_nt'] # Getting the peak size at max intensity
      p_start = df[ind_start, n] # Getting the peak intensity at start
      p_start_size = df[ind_start, 'size_nt'] # Getting the peak size at start 
      p_end = df[ind_end, n] # Getting the peak intensity at start
      p_end_size = df[ind_end, 'size_nt'] # Getting the peak size at start 
      p_halfmax = p_max/2 # Getting the intensity of the peak at half max
      # Now I need to make a line that calculates the width (in bp) from this halfmax measurement
      p_df = df[ind_start:ind_end, c('size_nt', n)]
      p_df_half = p_df[p_df[,2] >= p_halfmax,]
      p_width_halfmax = p_df_half[nrow(p_df_half),1] - p_df_half[1,1]
    } else {
      p_max = mean(df[range_deg_l[1]:range_deg_l[2], n])
      p_max_size = NA
      p_start = NA
      p_start_size = NA
      p_end = NA
      p_end_size = NA 
      p_halfmax = NA
      p_width_halfmax = 0}
    results = c(n, p_max, p_max_size, p_start, p_start_size, 
                p_end, p_end_size, p_halfmax, p_width_halfmax, "deg_28s_long")
    ana_28s_deg_long_list[[i]] <- results
  }
  
  # The smaller fragment of 28s degradation as the final peak
  ana_28s_deg_short_list <- list()
  for(i in seq_along(names(peaks_list_ri))){
    n = names(peaks_list_ri[i]) # Getting the sample name
    m = peaks_list_ri[[n]] # reading in the peak matrix
    if(find_peakl(m, lower_bound = range_deg_s[1], upper_bound = range_deg_s[2])){
      rnum = find_peaki(m, lower_bound = range_deg_s[1], upper_bound = range_deg_s[2])
      ind_max = m[rnum,2] # getting the peak index
      ind_start = m[rnum,3] # getting the index of the start of the peak
      ind_end = m[rnum,4] # getting the index of the end of the peak
      p_max = df[ind_max, n] # Getting the peak intensity at max
      p_max_size = df[ind_max, 'size_nt'] # Getting the peak size at max intensity
      p_start = df[ind_start, n] # Getting the peak intensity at start
      p_start_size = df[ind_start, 'size_nt'] # Getting the peak size at start 
      p_end = df[ind_end, n] # Getting the peak intensity at start
      p_end_size = df[ind_end, 'size_nt'] # Getting the peak size at start 
      p_halfmax = p_max/2 # Getting the intensity of the peak at half max
      # Now I need to make a line that calculates the width (in bp) from this halfmax measurement
      p_df = df[ind_start:ind_end, c('size_nt', n)]
      p_df_half = p_df[p_df[,2] >= p_halfmax,]
      p_width_halfmax = p_df_half[nrow(p_df_half),1] - p_df_half[1,1]
    } else {
      p_max = mean(df[range_deg_s[1]:range_deg_s[2], n])
      p_max_size = NA
      p_start = NA
      p_start_size = NA
      p_end = NA
      p_end_size = NA 
      p_halfmax = NA
      p_width_halfmax = 0}
    results = c(n, p_max, p_max_size, p_start, p_start_size, 
                p_end, p_end_size, p_halfmax, p_width_halfmax, "deg_28s_short")
    ana_28s_deg_short_list[[i]] <- results
  }
  
  allpeaks_df <- as.data.frame(do.call(rbind, c(ana_28s_list, ana_18s_list,
                                                ana_28s_deg_long_list, ana_28s_deg_short_list)))
  colnames(allpeaks_df) = c('sample_name', 'max_int', 'peak_position_nt', 'start_int', 'start_position_nt',
                            'end_int', 'end_position_nt', 'halfmax_int', 'width_halfmax_nt', 'id_peak' )
  allpeaks_df[,2:9] <- lapply(allpeaks_df[,2:9], function(x) as.numeric(x))
  allpeaks_df[, c(1,10)] <- lapply(allpeaks_df[, c(1,10)], function(x) as.factor(x) )
  return(allpeaks_df)
}



#' Return peaks from total RNA input
#'
#' This function analyzes peaks in rRNA samples from the provided input data frame and returns a list of peaks.
#' 
#' This is a good function for troubleshooting the pipeline in case it is giving weird results.
#' 
#' You can check the indices of the 18S and 28S rRNA against the indices in the functions and to see how consistent the data is across chips.
#'
#' @param input_df Either a file path to a CSV or a dataframe containing the data to analyze.
#' @param t_start The start time for trimming the data (default = 35).
#' @param t_end The end time for trimming the data (default = 55).
#' @param minpeakheight The minimum peak height for peak detection (default = 0.15).
#' @param minpeakdistance The minimum peak distance for peak detection (default = 40).
#' @param npeaks The number of peaks to detect (default = 4). If the data is visually noisy, increasing the npeaks here should help find the rRNA peaks
#'
#' @return A list of matrices, each containing the analyzed peaks.
#'
#' @examples
#' # Analyze rRNA peaks in vitro and return peak list
#' peak_list <- ba.analyze_rrna_list_only("path/to/input.csv")
#' peak_list
#'
#' @importFrom stats is.character is.data.frame
#' @importFrom dplyr bind_cols rownames
#' @export
ba.analyze_rrna_list_only <- function(input_df, t_start = 35, t_end = 55, minpeakheight = 0.15, minpeakdistance = 40, npeaks = 4){
  if (is.character(input_df)) {  # Check if input is a file path
    if (file.exists(input_df)) {
      df <- read.csv(input_df)
    } else {
      stop("File not found.")
    }
  } else if (is.data.frame(input_df)) {  # Check if input is a dataframe
    df <- input_df
    rownames(df) <- NULL # Reset the index so that this is the original reference to which to re-index
  } else {
    stop("Input must be either a file path to a CSV or a dataframe.")
  }
  if (ncol(df) > 3) { # Checking if there is more than one sample
    peaks_list <- apply(df[df$Time >= t_start & df$Time <= t_end,3:ncol(df)], 2, 
                        function(x) findpeaks.sort(x, minpeakheight = minpeakheight, minpeakdistance = minpeakdistance, npeaks = npeaks), 
                        simplify = FALSE)
  } else { # If there is only one sample, it should have only 3 columns
    peaks_list <- list()
    peaks_list[[colnames(df)[3]]] <- findpeaks.sort(df[df$Time >= t_start & df$Time <= t_end,3], minpeakheight = minpeakheight, minpeakdistance = minpeakdistance, npeaks = npeaks)
  }
  peaks_list_ri <- lapply(peaks_list, function(m) cbind(m[,1], m[,2:4] + (as.numeric(rownames(df[df$Time >= t_start & df$Time <= t_end,])[1]) -1)))
  return(peaks_list_ri)
}