#!/usr/bin/env Rscript
# ============================================================================
# View Work Log Summary
# Display recent work log entries
# ============================================================================

cat("=== FLVCR1 Work Log Viewer ===\n\n")

log_file <- "WORKLOG.md"

if (!file.exists(log_file)) {
  cat("No work log found (WORKLOG.md)\n")
  cat("Run: Rscript scripts/log_daily.R to create your first entry\n")
  quit(status = 1)
}

# Read the log file
log_content <- readLines(log_file)

# Get command line arguments for filtering
args <- commandArgs(trailingOnly = TRUE)

if (length(args) > 0 && args[1] == "--today") {
  # Show only today's entries
  today <- as.character(Sys.Date())
  in_today <- FALSE
  today_lines <- c()
  
  for (line in log_content) {
    if (grepl(paste0("^## ", today), line)) {
      in_today <- TRUE
    }
    if (in_today) {
      today_lines <- c(today_lines, line)
      if (grepl("^---$", line) && length(today_lines) > 5) {
        break
      }
    }
  }
  
  if (length(today_lines) > 0) {
    cat("=== Today's Log ===\n")
    cat(paste(today_lines, collapse = "\n"))
  } else {
    cat("No log entries for today yet.\n")
    cat("Run: Rscript scripts/log_daily.R to create an entry\n")
  }
  
} else if (length(args) > 0 && args[1] == "--last") {
  # Show last N days (default 7)
  n_days <- ifelse(length(args) > 1, as.numeric(args[2]), 7)
  
  cat(paste("=== Last", n_days, "Days ===\n\n"))
  
  # Find all date headers
  date_pattern <- "^## [0-9]{4}-[0-9]{2}-[0-9]{2}"
  date_lines <- grep(date_pattern, log_content)
  
  # Get last n_days entries
  recent_start <- max(1, length(date_lines) - n_days + 1)
  recent_dates <- date_lines[recent_start:length(date_lines)]
  
  if (length(recent_dates) > 0) {
    start_line <- recent_dates[1]
    cat(paste(log_content[start_line:length(log_content)], collapse = "\n"))
  } else {
    cat("No recent entries found.\n")
  }
  
} else if (length(args) > 0 && args[1] == "--stats") {
  # Show statistics
  cat("=== Work Log Statistics ===\n\n")
  
  # Count entries
  date_pattern <- "^## [0-9]{4}-[0-9]{2}-[0-9]{2}"
  n_entries <- length(grep(date_pattern, log_content))
  
  # Find date range
  dates <- gsub("^## (\\S+).*", "\\1", grep(date_pattern, log_content, value = TRUE))
  
  if (length(dates) > 0) {
    cat("Total log entries:", n_entries, "\n")
    cat("First entry:", dates[1], "\n")
    cat("Latest entry:", dates[length(dates)], "\n")
    
    # Count task types
    completed <- length(grep("^- ✅|^✅", log_content))
    in_progress <- length(grep("^- 🔄|^🔄", log_content))
    issues <- length(grep("^- ⚠️|^⚠️", log_content))
    
    cat("\nTask counts:\n")
    cat("  Completed: ", completed, "\n")
    cat("  In progress:", in_progress, "\n")
    cat("  Issues:    ", issues, "\n")
  } else {
    cat("No entries found.\n")
  }
  
} else {
  # Show full log
  cat("=== Full Work Log ===\n\n")
  cat(paste(log_content, collapse = "\n"))
  cat("\n\n=== Usage ===\n")
  cat("Rscript scripts/view_log.R              # View full log\n")
  cat("Rscript scripts/view_log.R --today      # View today's entries\n")
  cat("Rscript scripts/view_log.R --last 7     # View last 7 days\n")
  cat("Rscript scripts/view_log.R --stats      # View statistics\n")
}

cat("\n")

