
ld_snplist <- function (variants, with_alleles = TRUE, pop = "EUR",
                       opengwas_jwt = get_opengwas_jwt()) {

  res <- ieugwasr::api_query("ld/matrix", query = list(rsid = variants,pop = pop),
                   opengwas_jwt = opengwas_jwt) %>%
                   get_query_content()

  if (all(is.na(res)))
    stop("None of the requested SNPs were found")

  variants2 <- res$snplist
  res <- res$matrix
  res <- matrix(as.numeric(res), nrow(res), ncol(res))
  variants3 <- do.call(rbind, strsplit(variants2, split = "_"))

  if (with_alleles) {
    rownames(res) <- variants2
    colnames(res) <- variants2
  } else {
    rownames(res) <- variants3[, 1]
    colnames(res) <- variants3[, 1]
  }

  missing <- variants[!variants %in% variants3[, 1]]

  if (length(missing) > 0) {
    warning("The following variants are not present in the LD reference panel\n",
            paste(missing, collapse = "\n"))
    }

  ord <- match(variants3[, 1], variants)
  res <- res[order(ord), order(ord)]
  return(res)
}
