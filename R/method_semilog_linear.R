



find_log_cutoff <- function(x, y1, y2, x_max){
  x1 = log(500)
  x2 = x_max
  m = (y1-y2)/(log(500)-x_max)
  b = y1 - m*x1
  y_pred = (m*x) + b
  return(y_pred)
}

get_custom_slope <- function(y1, y2, x_max){
  x1 = log(500)
  x2 = x_max
  m = (y1-y2)/(log(500)-x_max)
  b = y1 - m*x1
  return(c("slope"=m, "y-intercept"=b))
}

## df must be mouse percent, graft cutoff, host cutoff columns only
assign_new_call <- function(df){
  names(df) <- c("host_percent", "graft_cut", "host_cut")
  df$call <- NA
  for(i in 1:dim(df)[1]){
    #print(df[i,"host_percent"])
    if(df[i,"host_percent"] <= df[i,"graft_cut"])
    {df[i,"call"] = "GRCh38"}
    else if(df[i,"host_percent"] >= df[i,"host_cut"])
    {df[i,"call"] = "mm10"}
    else
    {df[i,"call"] = "Multiplet"}
  }
  return(df$call)
}








gem_df <- prep_gem_counts(gem_data)

find_semilog_threshold <- function(gem_df){

  log_max = max(gem_df$totalReadsLog)
  gem_df$graft_cutpoint_log <- sapply(gem_df$totalReadsLog, find_log_cutoff, y1=.25, y2=.10, x_max=log_max)
  gem_df$host_cutpoint_log <- sapply(gem_df$totalReadsLog, find_log_cutoff, y1=.75, y2=.90, x_max=log_max)
  gem_df$SemiLog_Call <- assign_new_call(gem_df[,c("percentMouse", "graft_cutpoint_log", "host_cutpoint_log")])

  graft_line = get_custom_slope(y1=.25, y2=.10, x_max=log_max)
  host_line = get_custom_slope(y1=.75, y2=.90, x_max=log_max)

  return(gem_df)
}





gem_df$barcode <- row.names(gem_df)
for_gem_file <- gem_df[,c("barcode", "GRCh38_count", "mm10_count", "MESSY_Call")]
names(for_gem_file) <- c("barcode", "GRCh38", "mm10", "call")

write.csv(gem_df, file=paste0("MultipletCalls_Matrix_", title, ".csv"), quote = FALSE, row.names = FALSE)
write.csv(for_gem_file, file=paste0(save_dir, "gem_classification.csv"), quote = FALSE, row.names = FALSE)

# write out Loupe Annotations
write.csv(gem_df[,c("barcode","x10X_Call")], file=paste0(title, "_10x_gem_class_Annotations.csv"), quote = FALSE, row.names = FALSE)
write.csv(gem_df[,c("barcode","MESSY_Call")], file=paste0(title, "_MESSY_class_Annotations.csv"), quote = FALSE, row.names = FALSE)

