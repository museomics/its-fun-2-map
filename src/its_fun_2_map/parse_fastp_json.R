if (!require("jsonlite")) install.packages("jsonlite")
library(jsonlite)

json_parse <- function(json_paths){  
  data <- fromJSON(json_paths)

  df <- data.frame(
    file = basename(json_paths),

    before_total_reads = as.numeric(data$summary$before_filtering$total_reads),
    before_total_bases = as.numeric(data$summary$before_filtering$total_bases),
    before_q30_bases = as.numeric(data$summary$before_filtering$q30_bases),
    before_q30_rate = as.numeric(data$summary$before_filtering$q30_rate),
    before_gc_content = as.numeric(data$summary$before_filtering$gc_content),

    after_total_reads = as.numeric(data$summary$after_filtering$total_reads),
    after_total_bases = as.numeric(data$summary$after_filtering$total_bases),
    after_q30_bases = as.numeric(data$summary$after_filtering$q30_bases),
    after_q30_rate = as.numeric(data$summary$after_filtering$q30_rate),
    after_read1_mean_length = as.numeric(data$summary$after_filtering$read1_mean_length),

    after_read2_mean_length = ifelse(length(data$summary$after_filtering$read2_mean_length)!=0,
                                     as.numeric(data$summary$after_filtering$read2_mean_length), NA),
    after_gc_content = as.numeric(data$summary$after_filtering$gc_content),

    passed_filter_reads = as.numeric(data$filtering_result$passed_filter_reads),
    low_quality_reads = as.numeric(data$filtering_result$low_quality_reads),
    too_short_reads = as.numeric(data$filtering_result$too_short_reads),

    duplication_rate = as.numeric(data$duplication$rate),
    stringsAsFactors = TRUE
  )
}

json2csv <- function(json_dir, out_dir){
  json_paths <- list.files(path = json_dir, pattern = ".json", full.names = T)
  tables <- lapply(json_paths, json_parse)
  combined.df <- do.call(rbind , tables)
  
  # Split into 2 dataframes, trimmed data and merged data
  trimmed_df <- combined.df[grep("trim", combined.df$file),]
    trimmed_df$ID <- gsub("_trim.json", "", trimmed_df$file) 
  merged_df <- combined.df[grep("merge", combined.df$file),]
    merged_df$ID <- gsub("_merge.json", "", merged_df$file) 
  
  trimmed_out <- paste0(out_dir, "/trimmed_json_summary.csv")
  merged_out <- paste0(out_dir, "/merged_json_summary.csv")

  # Write out CSVs
  write.csv(trimmed_df, trimmed_out, row.names = FALSE)
  write.csv(merged_df, merged_out, row.names = FALSE)
  
}

## Example usage: 
## json2csv(json_dir= "fastp_outdir", out_name = "json2csv_output.csv")

