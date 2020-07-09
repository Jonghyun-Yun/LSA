coda_df <- function(coda.object,
                    parameters = NULL) {

    if (!coda::is.mcmc(coda.object) && !coda::is.mcmc.list(coda.object))
        stop("Not an mcmc or mcmc.list object")

    mat     <- as.matrix(coda.object, iter = TRUE, chain = TRUE)
    df      <- as.data.frame(mat)

    names(df)[names(df) == "CHAIN"] <- "chain"
    names(df)[names(df) == "ITER"]  <- "iter"

    if(is.null(parameters))
        out.df <- df

    if(!is.null(parameters))
        out.df <- subset(df, select = c("chain", "iter", parameters))

    out.df
}
