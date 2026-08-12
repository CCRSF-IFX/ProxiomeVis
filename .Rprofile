if (
  tolower(Sys.getenv("RENV_CONFIG_AUTOLOADER_ENABLED", unset = "false")) %in% c("true", "t", "1") &&
    file.exists("renv/activate.R")
) {
  source("renv/activate.R")
}
