"SScamp" <-
structure(function (input, Ths, alp, scal) 
{
    .expr1 <- ifelse(alp * input < 1, 1, alp * input)
    .expr3 <- .expr1^(-(scal))
    .value <- Ths * .expr3
    .actualArgs <- as.list(match.call()[c("Ths", "alp", "scal")])
    if (all(unlist(lapply(.actualArgs, is.name)))) {
        .grad <- array(0, c(length(.value), 3), list(NULL, c("Ths", 
            "alp", "scal")))
        .grad[, "Ths"] <- .expr3
        .grad[, "alp"] <- -(Ths * ((.expr1^(-((scal) + 1))) * 
            ((scal) * input)))
        .grad[, "scal"] <- -(Ths * (.expr3 * (log(.expr1))))
        dimnames(.grad) <- list(NULL, .actualArgs)
        attr(.value, "gradient") <- .grad
    }
    .value
}, initial = function (mCall, data, LHS) 
{
    xy <- data.frame(sortedXyData(mCall[["input"]], LHS, data))
    if (nrow(xy) < 5) {
        stop("Too few distinct input values to fit a van Genuchten model")
    }
    ndistinct <- nrow(xy)
    nlast <- max(3, round(ndistinct/2))
    Ths <- xy[1:(ndistinct - nlast), ][["y"]][1]
    dmid <- xy[2:(ndistinct - 1), ]
    pars2 <- coef(lm(log(x) ~ log(y), data = dmid))
    alp <- exp(pars2[1])
    scal <- -1/(pars2[2])
    value <- c(Ths = Ths, alp = alp, scal = scal)
    names(value) <- mCall[c("Ths", "alp", "scal")]
    value
}, pnames = c("Ths", "alp", "scal"), class = "selfStart")
