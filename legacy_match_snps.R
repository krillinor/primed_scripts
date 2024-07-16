library(data.table)

strand_flip = function(x) {
    fcase(
        x == "A" , "T",
        x == "C" , "G",
        x == "T" , "A",
        x == "G" , "C"
    )
}

match_snps = function(path_1kg, path_bim, build = "hg38") {
    allowed_builds = c("hg19", "hg38")
    if (!build %in% allowed_builds) {
        stop(paste("Build has to be ", paste(allowed_builds, collapse = " or "), sep = ""))
    }

    snp_col = paste("snp", build, sep = "_")

    matched_variants = fread(path_1kg)
    matched_variants[, c("chr", "pos", "a1", "a2") := tstrsplit(get(snp_col), ":", fixed = TRUE)]
    matched_variants[, pos := as.numeric(pos)]

    bim_col_names = c("chr", "id", "cm", "pos", "a1", "a2")
    bim = fread(path_bim, col.names = bim_col_names)
    bim[, chr := as.character(chr)]
    cols = c("a1", "a2")
    flip_cols = c("a1_flip", "a2_flip")
    bim[, (flip_cols) := lapply(.SD, strand_flip), .SDcols = cols]

    merged1 <- bim[matched_variants, on = .(chr, pos, a1, a2), nomatch = 0]
    merged2 <- bim[matched_variants, on = c("chr", "pos", a1 = "a2", a2 = "a1"), nomatch = 0]
    merged3 <- bim[matched_variants, on = c("chr", "pos", a1_flip = "a1", a2_flip = "a2"), nomatch = 0]
    merged4 <- bim[matched_variants, on = c("chr", "pos", a2_flip = "a1", a1_flip = "a2"), nomatch = 0]

    bim_merged = rbindlist(list(merged1, merged2, merged3, merged4))
    snp_target_col = paste("snp", "target", build, sep = "_")
    bim_merged[, (snp_target_col) := paste(chr, pos, a1, a2, sep = ":")]

    final_cols = c(snp_col, snp_target_col)
    bim_merged[, ..final_cols]
}
