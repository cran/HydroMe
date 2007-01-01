"SSvan" <-
structure(function (input, Thr, Ths, alp, scal) 
{
    .expr1 <- Ths - Thr
    .expr3 <- ifelse((exp(alp)) * input <= 0, 1, (exp(alp)) * 
        input)
    .expr4 <- .expr3^scal
    .expr5 <- 1 + .expr4
    .expr7 <- 1 - (1/scal)
    .expr8 <- .expr5^.expr7
    .expr11 <- 1/.expr8
    .expr14 <- .expr5^(.expr7 - 1)
    .expr22 <- .expr8^2
    .value <- Thr + (.expr1/.expr8)
    .actualArgs <- as.list(match.call()[c("Thr", "Ths", "alp", 
        "scal")])
    if (all(unlist(lapply(.actualArgs, is.name)))) {
        .grad <- array(0, c(length(.value), 4), list(NULL, c("Thr", 
            "Ths", "alp", "scal")))
        .grad[, "Thr"] <- 1 - .expr11
        .grad[, "Ths"] <- .expr11
        .grad[, "alp"] <- -((.expr1 * (.expr14 * (.expr7 * ((.expr3^(scal - 
            1)) * (scal * .expr3)))))/.expr22)
        .grad[, "scal"] <- -((.expr1 * ((.expr14 * (.expr7 * 
            (.expr4 * (log(.expr3))))) + (.expr8 * ((log(.expr5)) * 
            (1/(scal^2))))))/.expr22)
        dimnames(.grad) <- list(NULL, .actualArgs)
        attr(.value, "gradient") <- .grad
    }
    .value
}, initial = function (mCall, data, LHS) 
{
    xy <- data.frame(sortedXyData(mCall[["input"]], LHS, data))
    if (nrow(xy) < 4) {
        stop("Too few distinct input values to fit a van Genuchten model")
    }
    ndistinct <- nrow(xy)
    nlast <- max(3, round(ndistinct/2))
    dfirst <- xy[1, ][["y"]][1]
    dlast <- xy[nrow(xy), ][["y"]][1]
    Thr1 <- ifelse(((xy[1, ][["y"]] - xy[ndistinct, ][["y"]]) < 
        (xy[1, ][["y"]])/2), 0, xy[ndistinct, ][["y"]])
    Ths1 <- dfirst
    dmid <- xy[(ndistinct/2 - 2):(ndistinct/2 + 1), ]
    pars2 <- coef(lm(y ~ log(x), data = dmid))
    ymid <- xy[1:max(3, round(nrow(xy)/2)), ][["y"]][max(3, round(nrow(xy)/2))]
    ax <- (ymid - pars2[1])/pars2[2]
    slopep <- pars2[2]/(dfirst - dlast)
    m1 <- ifelse(abs(slopep) < 1, (1 - exp(-0.8 * (abs(slopep)))), 
        (1 - (0.5755/(abs(slopep))) + (0.1/(abs(slopep))^2) + 
            (0.025/(abs(slopep))^3)))
    scal1 <- 1/(1 - m1)
    alp1 <- (((2^(1/m1)) - 1)^(1 - m1))/exp(ax)
    pars <- c(Thr = Thr1, Ths = Ths1, alp = alp1, scal = scal1)
    val <- c(pars[1], pars[2], pars[3], pars[4])
    names(val) <- mCall[c("Thr", "Ths", "alp", "scal")]
    val
}, pnames = c("Thr", "Ths", "alp", "scal"), class = "selfStart")
