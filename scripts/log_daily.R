#!/usr/bin/env Rscript
# ============================================================================
# Daily Work Log Helper Script
# Automatically creates timestamped log entries
# ============================================================================

# Usage: Rscript scripts/log_daily.R "Your log message here"

cat("=== FLVCR1 Daily Log Helper ===\n\n")

# Get today's date
today <- Sys.Date()
timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

# Get log message from command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  cat("Interactive mode - Enter your daily log:\n")
  cat("(Press Enter after each item, type 'done' when finished)\n\n")
  
  completed <- c()
  in_progress <- c()
  issues <- c()
  notes <- c()
  
  # Collect completed tasks
  cat("\n✅ Completed tasks (type 'done' when finished):\n")
  repeat {
    task <- readline(prompt = "  - ")
    if (tolower(task) == "done" || task == "") break
    completed <- c(completed, task)
  }
  
  # Collect in-progress tasks
  cat("\n🔄 In progress (type 'done' when finished):\n")
  repeat {
    task <- readline(prompt = "  - ")
    if (tolower(task) == "done" || task == "") break
    in_progress <- c(in_progress, task)
  }
  
  # Collect issues
  cat("\n⚠️ Issues/Blockers (type 'done' when finished):\n")
  repeat {
    issue <- readline(prompt = "  - ")
    if (tolower(issue) == "done" || issue == "") break
    issues <- c(issues, issue)
  }
  
  # Collect notes
  cat("\n📝 Additional notes (type 'done' when finished):\n")
  repeat {
    note <- readline(prompt = "  - ")
    if (tolower(note) == "done" || note == "") break
    notes <- c(notes, note)
  }
  
  # Build log entry
  log_entry <- paste0(
    "\n---\n\n",
    "## ", today, "\n\n",
    "**Logged at:** ", timestamp, "\n\n"
  )
  
  if (length(completed) > 0) {
    log_entry <- paste0(log_entry, "### ✅ Completed\n")
    for (task in completed) {
      log_entry <- paste0(log_entry, "- ", task, "\n")
    }
    log_entry <- paste0(log_entry, "\n")
  }
  
  if (length(in_progress) > 0) {
    log_entry <- paste0(log_entry, "### 🔄 In Progress\n")
    for (task in in_progress) {
      log_entry <- paste0(log_entry, "- ", task, "\n")
    }
    log_entry <- paste0(log_entry, "\n")
  }
  
  if (length(issues) > 0) {
    log_entry <- paste0(log_entry, "### ⚠️ Issues/Blockers\n")
    for (issue in issues) {
      log_entry <- paste0(log_entry, "- ", issue, "\n")
    }
    log_entry <- paste0(log_entry, "\n")
  }
  
  if (length(notes) > 0) {
    log_entry <- paste0(log_entry, "### 📝 Notes\n")
    for (note in notes) {
      log_entry <- paste0(log_entry, "- ", note, "\n")
    }
    log_entry <- paste0(log_entry, "\n")
  }
  
} else {
  # Quick mode - just log the message
  message <- paste(args, collapse = " ")
  log_entry <- paste0(
    "\n---\n\n",
    "## ", today, "\n\n",
    "**", timestamp, ":** ", message, "\n\n"
  )
}

# Append to WORKLOG.md
log_file <- "WORKLOG.md"
if (!file.exists(log_file)) {
  cat("WORKLOG.md not found. Creating new file...\n")
  writeLines("# FLVCR1 Research Work Log\n", log_file)
}

# Append the entry
cat(log_entry, file = log_file, append = TRUE)

cat("\n✓ Log entry added to WORKLOG.md\n")
cat("\nYour entry:\n")
cat(log_entry)

# Optional: Auto-commit to git
auto_commit <- readline(prompt = "\nCommit this log to git? (y/n): ")
if (tolower(auto_commit) == "y") {
  system2("git", c("add", log_file))
  commit_msg <- paste0("Update work log: ", today)
  system2("git", c("commit", "-m", shQuote(commit_msg)))
  cat("✓ Changes committed\n")
  
  push_changes <- readline(prompt = "Push to GitHub? (y/n): ")
  if (tolower(push_changes) == "y") {
    system2("git", c("push", "origin", "main"))
    cat("✓ Changes pushed to GitHub\n")
  }
}

cat("\nDone!\n")

