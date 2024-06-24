.onLoad <- function(libname, pkgname) { #nolint
  backports::import(pkgname)
}

utils::globalVariables(c("LSU", "SSU", "col_character", "col_double", "cols",
                         "end", "is.flag", "pos", "start"))

# Dereplicate a collection of fastq.gz files,
# Use ITSx to locate ITS1, ITS2, and LSU regions
# Output separate files for each input file and region.

# The goal is to avoid running ITSx many times on identical sequences.

TF <- function(x) { #nolint
  if (isTRUE(x)) {
    return("T")
  } else if (isFALSE(x)) {
    return("F")
  }
  # raise a nice error.  We don't do this first because we expect to get a
  #valid argument
  assertthat::assert_that(assertthat::is.flag(x))
}

#' @importFrom magrittr "%>%"
#' @export
magrittr::`%>%`

#' Extract Ribosomal RNA Gene Regions from Eukaryotic DNA
#'
#' Calls the external program ITSx, which must be installed and on the path. For
#' more information on installation, algorithms, or options, see
#' http://microbiology.se/software/itsx/.
#'
#' For basic usage, this function only checks the arguments and gernerates the
#' function call.  However, to simplify integration with R-based workflows, it
#' is also possible to specify the \code{in_file} as one of several R classes
#' that hold DNA sequences , and by specifying a \code{read_function}, the
#' output file(s) can be automatically read into R and returned as members of a
#' \code{list}.
#'
#'
#' @param in_file Name of a fasta file to read sequences from. Alternatively, a
#'        member of classes \link[ShortRead]{ShortRead},
#'        \link[Biostrings]{DNAStringSet} or \link[seqinr]{SeqFastadna}, in
#'        which case the sequences will be written to a temporary .fasta file.
#' @param out_root Root name for all output files.  The default will put these
#'        in a pseudorandomly generated file name in the R temp directory, which
#'        is probably not what you want if you intend to access them later.
#' @param taxon The taxonomic group(s) to attempt to find rDNA from.  See the
#'        manual for ITSx for available options.  In contrast to the default for
#'        ITSx ("all") the default for this function is "fungi".
#' @param e_value E-value cutoff for inclusion in results.
#' @param s_value Score cutoff for inclusion in results.
#' @param n_value Number of domains required for inclusion in results.
#' @param selection_priority Priority for determining sequence origin.
#' @param search_eval Actual E-value cutoff for HMMER search.
#'        Only one of search_eval and search_score may be supplied.
#' @param search_score Actual score cutoff for HMMER search.
#'        Only one of search_eval and search_score may be supplied.
#' @param allow_single_domain Allow inclusion of sequences where only one domain
#'        match is found.  Either FALSE or a (double, integer) pair giving
#'        inclusion criteria for e_value and Score.
#' @param allow_reorder Allow inclusion of sequences where the domains do not
#'        occur in the expected order.
#' @param complement Search the reverse complement of each sequence also.
#' @param cpu Number of threads to use.
#' @param multi_thread Whether to use multiple threads or not.
#' @param heuristics Use heuristic filtering to speed up HMMER search.
#' @param nhmmer Use nhmmer instead of hmmsearch.
#' @param summary Output summary results in \code{[out_root].summary.txt}.
#' @param graphical Output graphical results in \code{[out_root].graph}.
#' @param fasta Output full ITS (ITS1 + 5.8S + ITS2) results in
#'        \code{[out_root].full.fasta}.
#' @param preserve Use original sequence headers instead of writing new ones.
#' @param save_regions Additional reagions to save output for; Options are one
#'        or more of \code{"SSU"}, \code{"ITS1"}, \code{"5.8S"}, \code{"ITS2"},
#'        \code{"LSU"}, or \code{"all"} or \code{"none"}. There will be output
#'        in \code{[out_root].SSU.fasta}, \code{[out_root].ITS1.fasta}, etc.
#' @param anchor Number of extra bases included at the beginning and end of each
#'        region. Only one of \code{anchor} and \code{require_anchor} may be
#'        given.
#' @param require_anchor As \code{anchor}, but the anchor bases are required to
#'        be present for the region to be included in output.
#'        Only one of \code{anchor} and \code{require_anchor} may be given.
#' @param only_full Limit output to full length ITS1 and ITS2 regions.
#' @param partial Save additional files for partial regions.  The argument give
#'        the minimum number of bases required.  These will be saved as
#'        \code{[out_root].full_and_partial.fasta} for the full ITS region, and
#'        as \code{[out_root].SSU.full_and_partial.fasta},
#'        \code{[out_root].ITS1.full_and_partial.fasta}, etc. for the other
#'        regions.
#' @param concat Output concatenated ITS1 and ITS2 regions in
#'        \code{[out_root].concat.fasta}.
#' @param minlen Minimum length for ITS regions to be included in \code{concat}.
#' @param positions Output a table of positions where HMMER matches were
#'        detected in \code{[out_root].positions.txt}.
#' @param table Output a table of ITS results in \code{[out_root].hmmer.table}.
#' @param detailed_results Output detailed results in
#'        \code{[out_root].extraction.results}.
#' @param not_found  Output a list of sequences for which no match was found in
#'        \code{[out_root]_no_detections.txt}.
#' @param truncate Remove ends of ITS sequences if they extend beyond the ITS
#'        region.
#' @param silent Supress printing of information to screen.
#' @param graph_scale Sets the scale of the graphical output.
#' @param save_raw Save raw data from searches in directory
#'        \code{[out_root]_ITSx_raw_output}.
#' @param read_function A function which can be used to read a fasta format
#'        file, such as \code{\link[ShortRead:readFasta]{ShortRead::readFasta}},
#'   \code{\link[Biostrings:readDNAStringSet]{Biostrings::readDNAStringSet}}, or
#'        \code{\link[seqinr:read.fasta]{seqinr::read.fasta}}. If a value is
#'        given, then all output files will be read and then deleted, using the
#'        given function to read fasta files.  If no fasta output is requested,
#'        then this behavior can be triggered by giving any other value, for
#'        example \code{TRUE} (but also \code{FALSE}!)
#'
#' @return The return value of the ITSx program, or if \code{read_function} is
#'         given, a list containing the contents of all files which were created
#'         (except the raw output).
#' @export
itsx <- function(in_file, out_root = tempfile("itsx"),
                 taxon = "fungi",
                 e_value = 1e-5,
                 s_value = 0,
                 n_value = 2,
                 selection_priority = c("score", "sum", "domains", "eval"),
                 search_eval = 0.01,
                 search_score = NULL,
                 allow_single_domain = c(1e-9, 0),
                 allow_reorder = FALSE,
                 complement = TRUE,
                 cpu = 1, multi_thread = cpu > 1, heuristics = FALSE,
                 nhmmer = FALSE, summary = TRUE,
                 graphical = TRUE, fasta = TRUE, preserve = FALSE,
                 save_regions = c("ITS1", "ITS2"),
                 anchor = 0, require_anchor = 0,
                 only_full = FALSE,
                 partial = 0, concat = FALSE, minlen = 0,
                 positions = TRUE, table = FALSE, detailed_results = FALSE,
                 not_found = TRUE, truncate = TRUE, silent = FALSE,
                 graph_scale = 0, save_raw = FALSE,
                 read_function = NULL) {

  # check arguments

  selection_priority <- match.arg(selection_priority)
  cat("selection priority:", selection_priority, "\n")
  save_regions <- match.arg(
     save_regions,
     c("ITS1", "ITS2", "SSU", "LSU", "5.8S", "all", "none"),
     several.ok = TRUE
  )

  assertthat::assert_that(
    assertthat::is.number(e_value),
    assertthat::is.count(s_value) ||
      assertthat::is.number(s_value) && s_value == 0,
    assertthat::is.count(n_value) ||
      assertthat::is.number(s_value) && s_value == 0,
    assertthat::is.string(selection_priority),
    missing(search_score) || missing(search_eval),
    missing(search_eval) ||
      (assertthat::is.number(search_eval) && search_eval >= e_value),
    missing(search_score) ||
      (assertthat::is.count(search_score) && search_score <= s_value),
    isFALSE(allow_single_domain) ||
      (is.numeric(allow_single_domain) &&
         length(allow_single_domain) == 2 &&
         assertthat::is.number(allow_single_domain[1]) &&
         allow_single_domain[1] <= e_value &&
         (assertthat::is.count(allow_single_domain[2]) ||
            allow_single_domain[2] == 0) &&
         allow_single_domain[2] >= s_value),
    assertthat::is.flag(allow_reorder),
    assertthat::is.flag(complement),
    assertthat::is.count(cpu) || assertthat::is.number(cpu) && cpu == 0,
    assertthat::is.flag(multi_thread),
    assertthat::is.flag(heuristics),
    assertthat::is.flag(nhmmer),
    assertthat::is.flag(summary),
    assertthat::is.flag(graphical),
    assertthat::is.flag(fasta),
    assertthat::is.flag(preserve),
    is.character(save_regions),
    if ("all" %in% save_regions) length(save_regions) == 1 else TRUE,
    if ("none" %in% save_regions) length(save_regions) == 1 else TRUE,
    missing(anchor) || assertthat::is.count(anchor) ||
      (assertthat::is.number(anchor) && anchor == 0) ||
      identical(anchor, "HMM"),
    missing(require_anchor) || assertthat::is.count(require_anchor) ||
      (assertthat::is.number(require_anchor) && require_anchor == 0) ||
      identical(require_anchor, "HMM"),
    assertthat::is.flag(only_full),
    assertthat::is.count(partial) ||
      assertthat::is.number(partial) && partial == 0,
    assertthat::is.flag(concat),
    assertthat::is.count(minlen) ||
      assertthat::is.number(minlen) && minlen == 0,
    assertthat::is.flag(positions),
    assertthat::is.flag(table),
    assertthat::is.flag(detailed_results),
    assertthat::is.flag(not_found),
    assertthat::is.flag(truncate),
    assertthat::is.flag(silent),
    assertthat::is.count(graph_scale) ||
      assertthat::is.number(graph_scale) && graph_scale == 0,
    assertthat::is.flag(save_raw)
  )
  if (methods::is(in_file, "ShortRead")) {
    .data <- in_file
    in_file <- paste0(out_root, ".fasta")
    ShortRead::writeFasta(.data, in_file = paste0(out_root, ".fasta"))
    on.exit(file.remove(in_file))
  } else if (methods::is(in_file, "DNAStringSet")) {
    .data <- in_file
    in_file <- paste0(out_root, ".fasta")
    Biostrings::writeXStringSet(.data, in_file)
    on.exit(file.remove(in_file))
  } else if (is.list(in_file) &&
             all(purrr::map_lgl(in_file, methods::is, "SeqFastadna"))) {
    .data <- in_file
    in_file <- paste0(out_root, ".fasta")
    seqinr::write.fasta(sequences = .data, names = seqinr::getName(.data),
                        file.out = in_file)
    on.exit(file.remove(in_file))
  } else {
    assertthat::assert_that(file.exists(in_file),
                            file.access(in_file, 4) == 0)
  }

  #required command line
  cl <- glue::glue("ITSx -i \"{in_file}\" -o \"{out_root}\"",
                   "-t \"{paste(taxon, collapse = \",\")}\"",
                   "-E {format(e_value, scientific = TRUE)}",
                   "-S {s_value} -N {n_value}",
                   "--selection_priority {selection_priority}",
                   if (missing(search_eval)) {
                     ""
                   } else {
                     "--search_eval {format(search_eval, scientific = TRUE)}"
                   },
                   if (missing(search_score)) {
                     ""
                   } else {
                     "--search_score {search_score}"
                   },
                   "--allow_single_domain",
                   if (isFALSE(allow_single_domain)) {
                     "F"
                   } else {
                     glue::double_quote(paste(allow_single_domain,
                                              collapse = ","))
                   },
                   "--allow_reorder {TF(allow_reorder)}",
                   "--complement {TF(complement)}",
                   "--cpu {cpu}",
                   "--multi_thread {TF(multi_thread)}",
                   "--heuristics {TF(heuristics)}",
                   "--nhmmer {TF(nhmmer)}",
                   "--summary {TF(summary)}",
                   "--graphical {TF(graphical)}",
                   "--fasta {TF(fasta)}",
                   "--preserve {TF(preserve)}",
                   "--save_regions {paste(save_regions, collapse = \",\")}",
                   if (missing(anchor)) {
                     ""
                   } else {
                     "--anchor {anchor}"
                   },
                   if (missing(require_anchor)) {
                     ""
                   } else {
                     "--require_anchor {require_anchor}"
                   },
                   "--only_full {TF(only_full)}",
                   "--partial {partial}",
                   "--concat {TF(concat)}",
                   "--minlen {minlen}",
                   "--positions {TF(positions)}",
                   "--table {TF(table)}",
                   "--detailed_results {TF(detailed_results)}",
                   "--not_found {TF(not_found)}",
                   "--truncate {TF(truncate)}",
                   "--silent {TF(silent)}",
                   "--graph_scale {graph_scale}",
                   "--save_raw {TF(save_raw)}",
                   .sep = " "
  )
  out <- system(cl)
  if (out > 0) stop("ITSx failed with return code ", out, ".")
  if (missing(read_function)) return(out)
  itsx_read_output(out_root, summary, graphical, fasta, save_regions, partial,
                   concat, positions, table, detailed_results, not_found,
                   read_function)
}

itsx_read_output <-  function(out_root, summary = TRUE,
                              graphical = TRUE, fasta = TRUE,
                              save_regions = c("ITS1", "ITS2"),
                              partial = 0, concat = FALSE,
                              positions = TRUE, table = FALSE,
                              detailed_results = FALSE,
                              not_found = TRUE,
                              read_function) {
  out <- list()
  do_partial <- partial > 0
  if (isTRUE(fasta)) {
    readfile <- paste0(out_root, ".full.fasta")
    if (file.exists(readfile)) {
      out$full <- read_function(readfile)
      file.remove(readfile)
    } else {
      out$full <- NULL
    }
    if (isTRUE(do_partial)) {
      readfile <- paste0(out_root, ".full_and_partial.fasta")
      if (file.exists(readfile)) {
        out$full_partial <- read_function(readfile)
        file.remove(readfile)
      } else {
        out$full_partial <- NULL
      }
    }
  }

  if ("all" %in% save_regions) {
    out_regions <- c("ITS1", "ITS2", "SSU", "LSU", "5.8S")
  } else if ("none" %in% save_regions) {
    out_regions <- character(0)
  } else {
    out_regions <- save_regions
  }
  for (r in out_regions) {
    readfile <-
       paste0(out_root, ".", stringr::str_replace_all(r, "\\.", "_"), ".fasta")
    if (file.exists(readfile)) {
      out[[r]] <- read_function(readfile)
      file.remove(readfile)
    } else {
      out[[r]] <- NULL
    }
    if (isTRUE(do_partial)) {
      readfile <- paste0(out_root, ".full_and_partial.", r, ".fasta")
      if (file.exists(readfile)) {
        out[[paste0(r, ".full.partial")]] <- read_function(readfile)
        file.remove(readfile)
      } else {
        out[[paste0(r, ".full.partial")]] <- NULL
      }
    }
  }
  if (isTRUE(summary)) {
    readfile <- paste0(out_root, ".summary.txt")
    if (file.exists(readfile)) {
      out$summary <- readLines(readfile)
      file.remove(readfile)
    } else {
      out$summary <- NULL
    }
  }
  if (isTRUE(graphical)) {
    readfile <- paste0(out_root, ".graph")
    if (file.exists(readfile)) {
      out$graphical <- readLines(readfile)
      file.remove(readfile)
    } else {
      out$graphical <- NULL
    }
  }
  if (isTRUE(concat)) {
    readfile <- paste0(out_root, ".concat.fasta")
    if (file.exists(readfile)) {
      out$concat <- read_function(readfile)
      file.remove(readfile)
    } else {
      out$concat <- NULL
    }
  }
  if (isTRUE(positions)) {
    readfile <- paste0(out_root, ".positions.txt")
    if (file.exists(readfile)) {
      out$positions <-
        readr::read_tsv(readfile,
                        col_names = c("seq_id", "length", "SSU", "ITS1",
                                      "5_8S", "ITS2", "LSU", "comment"),
                        col_types = readr::cols(
                           seq_id = readr::col_character(),
                           length = readr::col_character(),
                           SSU = readr::col_character(),
                           ITS1 = readr::col_character(),
                           `5_8S` = readr::col_character(),
                           ITS2 = readr::col_character(),
                           LSU = readr::col_character(),
                           comment = readr::col_character()
                        )) %>%
        tidyr::gather(key = "region", value = "pos", SSU:LSU) %>%
        dplyr::mutate_at("pos", stringr::str_extract, pattern = "\\d+-\\d+") %>%
        tidyr::extract(
           pos,
           into = c("start", "end"),
           regex = "(\\d+)-(\\d+)"
        ) %>%
        dplyr::mutate_at(dplyr::vars(start:end), as.integer)
      file.remove(readfile)
    } else {
      out$positions <- NULL
    }
  }
  if (isTRUE(table)) {
    readfile <- paste0(out_root, ".hmmer.table")
    if (file.exists(readfile)) {
      out$table <- readr::read_tsv(readfile)
      file.remove(readfile)
    } else {
      out$table <- NULL
    }
  }
  if (isTRUE(detailed_results)) {
    readfile <- paste0(out_root, ".extraction.results")
    if (file.exists(readfile)) {
      out$detailed_results <- readr::read_tsv(readfile)
      file.remove(readfile)
    } else {
      out$detailed_results <- NULL
    }
  }
  if (isTRUE(not_found)) {
    readfile <- paste0(out_root, "_no_detections.txt")
    if (file.exists(readfile)) {
      out$not_found <- readLines(readfile)
      file.remove(readfile)
    } else {
      out$not_found <- NULL
    }
  }
  out
}
