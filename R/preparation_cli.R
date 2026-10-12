parse_preparation_cli <- function(args, options, positional_counts, usage) {
  positional <- character()
  supplied <- character()
  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]
    if (startsWith(arg, "--")) {
      name <- sub("^--", "", arg)
      if (!name %in% names(options)) stop("Unknown option: ", arg, "\n", usage, call. = FALSE)
      if (name %in% supplied) stop("Repeated option: ", arg, call. = FALSE)
      if (i == length(args) || startsWith(args[[i + 1L]], "--") || !nzchar(args[[i + 1L]])) {
        stop("Supply a value after ", arg, call. = FALSE)
      }
      supplied <- c(supplied, name)
      options[[name]] <- args[[i + 1L]]
      i <- i + 2L
    } else {
      positional <- c(positional, arg)
      i <- i + 1L
    }
  }
  if (!length(positional) %in% positional_counts) stop(usage, call. = FALSE)
  list(positional = positional, options = options)
}
