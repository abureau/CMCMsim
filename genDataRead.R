genDataRead = function (file.in = stop("Filename must be given!", call. = FALSE), 
    file.out = NULL, dir.out = ".", format = stop("Format parameter is required!"), 
    header = FALSE, n.vars, cov.file.in, cov.header, map.file, 
    allele.sep = ";", na.strings = "NA", col.sep = "", overwrite = NULL) 
{
#    if (!file.exists(file.in)) {
#        stop("The given file (", file.in, ") doesn't exist! Check and try again.", 
#            call. = FALSE)
#    }
    if (!missing(cov.file.in)) {
        if (!file.exists(cov.file.in)) {
            stop("The given file (", cov.file.in, ") doesn't exist! Check and try again.", 
                call. = FALSE)
        }
    }
    if (!missing(map.file)) {
        if (!file.exists(map.file)) {
            stop("The given map.file (", map.file, ") doesn't exist! Check and try again.", 
                call. = FALSE)
        }
    }
    map.file <- NULL
    files.list <- Haplin:::f.make.out.filename(file.in, file.out, dir.out = dir.out, 
        overwrite = overwrite)
    if (!(format %in% c("haplin", "ped"))) {
        stop("Given format (", format, ") not recognized! Format has to be one of: ", 
            paste(c("haplin", "ped"), collapse = ", "), call. = FALSE)
    }
    if (format == "haplin") {
        if ((na.strings == col.sep) | ((na.strings != "") & (na.strings == 
            allele.sep))) {
            stop("The \"na.strings\" argument should be different from \"col.sep\" and \"allele.sep\" arguments!", 
                call. = F)
        }
        split <- F
        if ((col.sep != "" & col.sep == allele.sep) | allele.sep == 
            " ") {
            split <- T
        }
    }
    if (!missing(cov.header)) {
        if (!is.character(cov.header)) {
            stop("'cov.header' is specified, but it's not a character vector!", 
                call. = FALSE)
        }
    }
    skip.first <- 0
    if (header) {
        header.line <- scan(file = file.in, what = "character", 
            nlines = 1, strip.white = TRUE, sep = col.sep, na.strings = na.strings)
        skip.first <- 1
    }
    first.line <- scan(file = file.in, what = "character", nlines = 1, 
        strip.white = TRUE, sep = col.sep, na.strings = na.strings, 
        skip = skip.first)
    nb.cols <- length(first.line)
    if (missing(n.vars)) {
        if (format == "haplin") {
            stop("The value of \"n.vars\", i.e. the number\n of covariate columns ahead of genetic data,\n should be specified when format = \"haplin\".\n Set n.vars = 0 if there are no covariates in haplin file.", 
                call. = FALSE)
        }
        if (format == "ped") {
            n.vars <- 6
        }
        cat("'n.vars' was not given explicitly and will be set to ", 
            n.vars, " based on the format given.\n", sep = "")
    }
    else if (n.vars < 0 | n.vars >= nb.cols) {
        stop("The 'n.vars' value (", n.vars, ") is not valid!", 
            call. = FALSE)
    }
    if (format == "haplin" & n.vars > 0 & missing(cov.header)) {
        cat("The format of the file is 'haplin' with covariate data but no names of the covariate data is given. Will generate dummy names.\n")
    }
    nb.lines.per.chunk <- ceiling(1e+08/nb.cols)
#    in.file <- file(file.in, "r")
    in.file <- file.in
    gen.data.in.ffdf <- list()
    cov.data.in <- c()
    gen.levels <- c()
    nb.rows.tot <- 0
    i <- 1
    cat("Reading the data in chunks...\n")
    cat(" -- chunk ", i, " --\n", sep = "")
    cur.chunk <- matrix(scan(in.file, what = "character", nlines = nb.lines.per.chunk, 
        sep = col.sep, na.strings = na.strings), ncol = nb.cols, 
        byrow = TRUE)
    while (length(cur.chunk) != 0) {
        if (n.vars > 0) {
            cov.data.in <- rbind(cov.data.in, cur.chunk[, 1:n.vars, 
                drop = F])
            cur.chunk <- cur.chunk[, -(1:n.vars), drop = F]
        }
        if (format == "haplin") {
            if (!split) {
                cur.chunk <- Haplin:::f.split.matrix(cur.chunk, split = allele.sep)
            }
            else {
                d <- dim(cur.chunk)
                .data.gen1 <- cur.chunk[, seq(1, d[2], 2), drop = F]
                .data.gen2 <- cur.chunk[, seq(2, d[2], 2), drop = F]
                .is.na <- is.na(.data.gen1) | is.na(.data.gen2)
                .data.gen1[.is.na] <- NA
                .data.gen2[.is.na] <- NA
                .d <- dim(.data.gen1) * c(1, 2)
                cur.chunk <- matrix(NA_character_, nrow = d[1], 
                  ncol = d[2])
                cur.chunk[, seq(1, d[2], 2)] <- .data.gen1
                cur.chunk[, seq(2, d[2], 2)] <- .data.gen2
                rm(.data.gen1)
                rm(.data.gen2)
            }
        }
        cur.levels <- unique(as.vector(cur.chunk))
        recode.na <- FALSE
        if (0 %in% cur.levels) {
            recode.na <- TRUE
            na.symbol <- 0
        }
        else if ("0" %in% cur.levels) {
            recode.na <- TRUE
            na.symbol <- "0"
        }
        if (recode.na) {
            cur.chunk[cur.chunk == na.symbol] <- NA
            cur.levels[cur.levels == na.symbol] <- NA
        }
        gen.levels <- as.character(union(gen.levels, cur.levels))
        tmp.ff <- ff::as.ff(cur.chunk, vmode = Haplin:::.haplinEnv$.vmode.gen.data, 
            levels = gen.levels)
        gen.data.in.ffdf <- c(gen.data.in.ffdf, list(tmp.ff))
        nb.rows.tot <- nb.rows.tot + nrow(tmp.ff)
        rm(cur.chunk, tmp.ff)
        i <- i + 1
        cat(" -- chunk ", i, " -- \n", sep = "")
        cur.chunk <- matrix(scan(in.file, what = "character", 
            nlines = nb.lines.per.chunk), ncol = nb.cols, byrow = TRUE)
    }
    cat("... done reading.\n")
    close(in.file)
    nb.cols.per.chunk <- get(".nb.cols.per.chunk", envir = Haplin:::.haplinEnv)
    nb.cols.gen.data <- ncol(gen.data.in.ffdf[[1]])
    nb.col.chunks <- ceiling(nb.cols.gen.data/nb.cols.per.chunk)
    gen.list.length <- length(gen.data.in.ffdf)
    gen.data.col.wise <- list()
    design <- "cc"
    if (format == "haplin") {
        design <- "triad"
    }
    tot.gen.ncol <- ncol(gen.data.in.ffdf[[1]])
    gen.data.colnames <- Haplin:::f.create.snp.names(map.file, ncol = tot.gen.ncol, 
        format = format, design = design)
    marker.names <- gen.data.colnames$marker.names
    gen.data.colnames <- gen.data.colnames$gen.data.colnames
    cat("Preparing data...\n")
    for (i in 1:nb.col.chunks) {
        cur.cols <- ((i - 1) * nb.cols.per.chunk + 1):(min(i * 
            nb.cols.per.chunk, nb.cols.gen.data))
        tmp.gen.data <- ff::ff(vmode = Haplin:::.haplinEnv$.vmode.gen.data, 
            levels = gen.levels, dim = c(nb.rows.tot, min(nb.cols.per.chunk, 
                max(cur.cols) - min(cur.cols) + 1)))
        prev.rows <- 0
        for (j in 1:gen.list.length) {
            cur.rows <- (prev.rows + 1):(prev.rows + nrow(gen.data.in.ffdf[[j]]))
            tmp.gen.data[cur.rows, ] <- gen.data.in.ffdf[[j]][, 
                cur.cols]
            prev.rows <- max(cur.rows)
        }
        colnames(tmp.gen.data) <- gen.data.colnames[cur.cols]
        gen.data.col.wise <- c(gen.data.col.wise, list(tmp.gen.data))
        rm(tmp.gen.data)
    }
    cat("... done preparing\n")
    rm(gen.data.in.ffdf)
    cov.data.colnames <- c()
    if (n.vars > 0) {
        if (format == "ped") {
            cov.data.colnames <- get(".cov.data.colnames", envir = Haplin:::.haplinEnv)
        }
        else if (!missing(header)) {
            cov.data.colnames <- header.line[1:n.vars]
        }
        else {
            cov.data.colnames <- paste0("cov.", 1:n.vars)
        }
    }
    cov.n.vars <- 0
    if (!missing(cov.file.in)) {
        cat("Reading covariate file... \n")
        if (missing(cov.header)) {
            cat("    'cov.header' not given - assuming the first line is the header...\n")
            cov.file.header <- TRUE
        }
        else {
            cov.file.header <- FALSE
        }
        cov.add.data <- read.table(cov.file.in, header = cov.file.header, 
            stringsAsFactors = FALSE)
        if (missing(cov.header)) {
            cov.header <- colnames(cov.add.data)
        }
        cov.data.colnames <- c(cov.data.colnames, cov.header)
        if (nrow(cov.add.data) != nb.rows.tot) {
            stop("The number of rows in the additional covariate data (", 
                nrow(cov.add.data), ") doesn't match the number of rows in the main file with genetic data (", 
                nb.rows.tot, ")!", call. = FALSE)
        }
        if (!cov.file.header & (ncol(cov.add.data) != length(cov.header))) {
            stop("The length of given 'cov.header' names (", 
                length(cov.header), ") doesn't match the number of all the covariate data columns (", 
                ncol(cov.add.data), ")! Check and try again.")
        }
        if (length(cov.data.in) != 0) {
            cov.data.in <- cbind(cov.data.in, cov.add.data)
        }
        else {
            cov.data.in <- cov.add.data
        }
        cov.n.vars <- ncol(cov.add.data)
        cat("...done\n")
    }
    else if (!missing(cov.header)) {
        if (format == "haplin") {
            cov.n.vars <- 0
            if (n.vars != length(cov.header)) {
                stop("The length of given 'cov.header' names (", 
                  length(cov.header), ") doesn't match the number of all the covariate data columns (", 
                  n.vars, ")! Check and try again.")
            }
            else {
                cov.data.colnames <- cov.header
            }
        }
        else if (format == "ped") {
            cov.n.vars <- n.vars - 6
            if (cov.n.vars != length(cov.header)) {
                stop("The length of given 'cov.header' names (", 
                  length(cov.header), ") doesn't match the number of all the covariate data columns (", 
                  cov.n.vars, ")! Check and try again.")
            }
            else {
                cov.data.colnames[7:n.vars] <- cov.header
            }
        }
    }
    n.vars <- n.vars + cov.n.vars
    cat("Saving data...\n")
    cur.names <- c()
    for (i in 1:length(gen.data.col.wise)) {
        cur.name <- paste(get(".gen.cols.name", envir = Haplin:::.haplinEnv), 
            i, sep = ".")
        assign(cur.name, gen.data.col.wise[[i]])
        cur.names <- c(cur.names, cur.name)
    }
    save.list <- cur.names
    if (n.vars > 0) {
        if (n.vars == 1) {
            cov.data.in <- matrix(cov.data.in, ncol = 1)
        }
        colnames(cov.data.in) <- cov.data.colnames
        save.list <- c(save.list, "cov.data.in")
    }
    else {
        cov.data.in <- NULL
    }
    .info <- list()
    .info$filename <- c(file.in = file.in)
    if (!missing(cov.file.in)) {
        .info$filename <- c(.info$filename, cov.file.in = cov.file.in)
    }
    .info$filespecs$n.vars <- n.vars
    .info$filespecs$sep <- col.sep
    .info$filespecs$allele.sep <- allele.sep
    .info$filespecs$na.strings <- na.strings
    .info$filespecs$format <- format
    aux <- list(info = .info, class = "haplin.data")
    aux$marker.names <- marker.names
    aux$map.filename <- map.file
    save.list <- c(save.list, "aux")
    ff::ffsave(list = save.list, file = file.path(dir.out, files.list$file.out.base))
    cat("... saved to files: ", files.list$file.out.ff, ", ", 
        files.list$file.out.aux, "\n", sep = "")
    data.out <- list(cov.data = cov.data.in, gen.data = gen.data.col.wise, 
        aux = aux)
    class(data.out) <- aux$class
    return(data.out)
}
