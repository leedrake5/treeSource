###Copy and paste this script into the R console on your computer. It will be compatible with Windows, Mac, and Linux.

###All sentences beginning with "###" will be invisible to the software, and will caption and describe each step of the analysis and figures for [CITATION]

###Erase everything that comes before
rm(list = ls(all = TRUE))

###Compatibility
if(.Platform$OS.type=="windows") {
    quartz<-function() windows()
}

###IMPORTANT NOTE: R uses packages to facilitate the analysis of data and production of figures. If you do not have the TTR, bcp, or ggplot2 packages installed, the following three lines of text will do it for you - all you have to do is delete the "###" that precedes the commands

###The command below will bring up a list of download sites. Pick one closest to you to speed up the download process
###chooseCRANmirror()

###The script below will then install the TTR package (for moving averages), bcp package (for Bayesian Change-Point analysis), and ggplot2 (for generating data plots)
###Note: Installation of packages may take up to an hour, depending upon the speed of your internect connection.
###install.packages("TTR", dependencies = TRUE)
###install.packages("bcp", dependencies = TRUE)
###install.packages("ggplot2", dependencies = TRUE)

###Activate the packages

library(ggplot2)
library(gridExtra)
library(dplR)
library(pbapply)
library(reshape)
library(reshape2)
library(Biobase)
library(xlsx)
library(forecast)
library(ggmap)
library(plyr)
library(akima)


#####Functions


####Function to organize plots in a window
layOut = function(...) {
    
    require(grid)
    
    x <- list(...)
    n <- max(sapply(x, function(x) max(x[[2]])))
    p <- max(sapply(x, function(x) max(x[[3]])))
    pushViewport(viewport(layout = grid.layout(n, p)))
    
    for (i in seq_len(length(x))) {
        print(x[[i]][[1]], vp = viewport(layout.pos.row = x[[i]][[2]],
        layout.pos.col = x[[i]][[3]]))
    }
}

sourcePrepare <- function(source.dataframe) {
    
    year <- source.dataframe$Year
    
    
}


lm.dat <- function (formula, data, subset, weights, na.action, method = "qr",
model = TRUE, x = FALSE, y = FALSE, qr = TRUE, singular.ok = TRUE,
contrasts = NULL, offset, ...)
{
    dat.fram <- data.frame(x, y)
    dat.fram <- dat.fram[complete.cases(dat.fram),]
    x <- dat.fram$x
    y <- dat.fram$y
    ret.x <- x
    ret.y <- y
    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action",
    "offset"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())
    if (method == "model.frame")
    return(mf)
    else if (method != "qr")
    warning(gettextf("method = '%s' is not supported. Using 'qr'",
    method), domain = NA)
    mt <- attr(mf, "terms")
    y <- model.response(mf, "numeric")
    w <- as.vector(model.weights(mf))
    if (!is.null(w) && !is.numeric(w))
    stop("'weights' must be a numeric vector")
    offset <- as.vector(model.offset(mf))
    if (!is.null(offset)) {
        if (length(offset) != NROW(y))
        stop(gettextf("number of offsets is %d, should equal %d (number of observations)",
        length(offset), NROW(y)), domain = NA)
    }
    if (is.empty.model(mt)) {
        x <- NULL
        z <- list(coefficients = if (is.matrix(y)) matrix(, 0,
        3) else numeric(), residuals = y, fitted.values = 0 *
        y, weights = w, rank = 0L, df.residual = if (!is.null(w)) sum(w !=
        0) else if (is.matrix(y)) nrow(y) else length(y))
        if (!is.null(offset)) {
            z$fitted.values <- offset
            z$residuals <- y - offset
        }
    }
    else {
        x <- model.matrix(mt, mf, contrasts)
        z <- if (is.null(w))
        lm.fit(x, y, offset = offset, singular.ok = singular.ok,
        ...)
        else lm.wfit(x, y, w, offset = offset, singular.ok = singular.ok,
        ...)
    }
    class(z) <- c(if (is.matrix(y)) "mlm", "lm")
    z$na.action <- attr(mf, "na.action")
    z$offset <- offset
    z$contrasts <- attr(x, "contrasts")
    z$xlevels <- .getXlevels(mt, mf)
    z$call <- cl
    z$terms <- mt
    if (model)
    z$model <- mf
    if (ret.x)
    z$x <- x
    if (ret.y)
    z$y <- y
    if (!qr)
    z$qr <- NULL
    z
}

unlist.tree <- function(temp, myfiles){
    n <- length(temp)
    for (i in n) {
        temp[i] <- myfiles[[i]]
    }
}

ig.na <- function(x) {
    length(na.omit(x))
}

nonNAs <- function(x) {
    n <- as.vector(apply(x, 1, function(x) length(which(!is.na(x)))))
    return(n)
}

readRWL.simp <- function(file) {
    raw <- read.rwl(file)
    years <- rownames(raw)
    non.total <- data.frame(years, raw)
    colnames(non.total)[1] <- "Year"
    return(non.total)
}

readRWLArima <- function(file) {
    raw <- read.rwl(file)
    raw <- read.rwl(file)
    years <- rownames(raw)
    detrended <- data.frame(detrend(raw, method="Spline", nyrs=50))
    list.arima <- pbapply(X=detrended, MARGIN=2, FUN=auto.arima)
    arima.data <- data.frame(subListExtract(L=list.arima, name="residuals"))
    arima.total <- data.frame(years, arima.data)
    colnames(arima.total)[1] <- "Year"
    return(arima.total)
}

readDataArima <- function(file) {
    raw <- read.csv(file)
    n <- length(raw)
    years <- raw[,1]
    data <- raw[,2:n]
    detrended <- data.frame(detrend(data, method="Spline", nyrs=50))
    arima.data <- pbapply(X=detrended, MARGIN=2, function(x) FUN=auto.arima(x)$residuals)
    arima.total <- data.frame(years, arima.data)
    colnames(arima.total)[1] <- "Year"
    return(arima.total)
}


readDataArimaFit <- function(file) {
    raw <- read.csv(file)
    n <- length(raw)
    years <- raw[,1]
    data <- raw[,2:n]
    detrended <- data.frame(detrend(data, method="Spline", nyrs=50))
    list.arima <- pbapply(X=detrended, MARGIN=2, function(x) FUN=fitted(auto.arima(x)))
    arima.total <- data.frame(years, list.arima)
    colnames(arima.total)[1] <- "Year"
    return(arima.total)
}


readDataArima4 <- function(file) {
    raw <- read.csv(file)
    n <- length(raw)
    years <- raw[,1]
    data <- raw[,2:n]
    detrended <- data.frame(detrend(data, method="Spline"))
    list.arima <- pbapply(X=detrended, MARGIN=2, FUN=auto.arima)
    arima.data <- data.frame(subListExtract(L=list.arima, name="x"))
    arima.total <- data.frame(years, arima.data)
    colnames(arima.total)[1] <- "Year"
    return(arima.total)
}


meanSequence <- function(tree.dataframe, name) {
    n <- length(tree.dataframe)
    sequ <- rowMeans(tree.dataframe[2:n], na.rm=TRUE)
    results.frame <- data.frame(as.numeric(as.vector(tree.dataframe$Year)), sequ)
    colnames(results.frame) <- c("Year", name)
    return(results.frame)
}

treeHypothesis <- function(timemin, timemax, tree.dataframe.1, tree.dataframe.2) {
    
    
    
    tree.a <- subset(tree.dataframe, !(tree.dataframe[,1] < timemin | tree.dataframe[,1] > timemax))
    tree.b <- subset(tree.source, !(tree.source[,1] < timemin | tree.source[,1] > timemax))
    
    
    
    df <- data.frame(time, tree.a$Mean, tree.b$Mean, tree.a$SD, tree.b$SD, tree.a$N, tree.b$N)
    colnames(df) <- c("Year", "FirstMean", "SecondMean", "FirstSD", "SecondSD", "FirstN", "SecondN")
    
    
    df$Ttest <- c(abs(df$FirstMean-df$SecondMean)/(sqrt((df$FirstSD^2)/df$FirstN + (df$SecondSD^2)/df$SecondN)))
    
    
    df$DF <- c(((((df$FirstSD^2)/df$FirstN) +  ((df$SecondSD^2)/df$SecondN))^2)/((df$FirstSD^4)/((df$FirstN^2)*(df$FirstN-1)) + (df$SecondSD^4)/((df$SecondN^2)*(df$SecondN-1))))
    
    
    df$pvalue <- c((2*pt(df$Ttest, df$DF, lower=FALSE)))
    df$Significant <- rep("Yes", length(df$Year))
    df <- transform(df, Significant = ifelse(pvalue > 0.05, "No", Significant))
    
    
    return(df)
}

sigCount <- function(tree.hypothesis.test.results) {
    Yes <- subset(tree.hypothesis.test.results$pvalue, tree.hypothesis.test.results$Significant=="Yes")
    No <- subset(tree.hypothesis.test.results$pvalue, tree.hypothesis.test.results$Significant=="No")
    
    Yess <- length(Yes)
    Nos <- length(No)
    
    results <- data.frame(mean(Yes), mean(No), Yess, Nos)
    colnames(results) <- c("p-value diff", "p-value same", "p < 0.05", "p > 0.05")
    return(results)
    
}


treeCorTest <- function(timemin, timemax, tree.object, tree.source){
    
    
    
    tree.a <- subset(tree.dataframe, !(tree.dataframe[,1] < timemin | tree.dataframe[,1] > timemax))
    tree.b <- subset(tree.source, !(tree.source[,1] < timemin | tree.source[,1] > timemax))
    
    tree.a <- tree.a[complete.cases(tree.a), ]
    tree.b <- tree.b[complete.cases(tree.b), ]
    
    tree.a.arima <-arima(tree.a[,2], order=c(1, 0 ,1))
    tree.a.n <- tree.a.arima$residuals
    
    #tree.b.arima <-arima(tree.b[,2], order=c(1, 0 ,1))
    #tree.b.n <- tree.b.arima$residuals
    
    tree.a.frame <- data.frame(tree.a$Year, tree.a.n)
    colnames(tree.a.frame) <- c("Year", "A")
    tree.b.frame <- data.frame(tree.b$Year, tree.b[,2])
    colnames(tree.b.frame) <- c("Year", "B")
    
    tree.a.re <- tree.a.frame$A[tree.a.frame$Year %in% tree.b$Year]
    tree.b.re <- tree.b.frame$B[tree.b.frame$Year %in% tree.a$Year]
    
    
    
    trees.grid <- data.frame(tree.a.re, tree.b.re)
    colnames(trees.grid) <- c("First", "Second")
    
    
    trees.lm <- lm(trees.grid$First~trees.grid$Second)
    trees.s.lm <- summary(trees.lm)
    
    turn.to.t <- function(x.lm) {
        x.s.lm <- summary(x.lm)
        r.sq <- x.s.lm$r.squared
        just.r <- sqrt(r.sq)
        t <- (just.r*sqrt(length(x.lm$residuals)-2))/sqrt(1-r.sq)
        return(t)
    }
    
    
    trees.t <- turn.to.t(trees.lm)
    
    result.frame <- data.frame(trees.t,sqrt(trees.s.lm$r.squared), length(trees.lm$residuals))
    colnames(result.frame) <- c("t", "r", "overlap")
    return(result.frame)
}


treeCorTestMultiple <- function(timemin, timemax, tree.object, tree.sources){
    
    tree.sources.n <- length(tree.sources)
    
    tree.a <- subset(tree.dataframe, !(tree.dataframe[,1] < timemin | tree.dataframe[,1] > timemax))
    tree.b <- subset(tree.source, !(tree.source[,1] < timemin | tree.source[,1] > timemax))
    
    tree.a <- tree.a[complete.cases(tree.a), ]
    tree.b <- tree.b[complete.cases(tree.b), ]
    
    tree.a.arima <-arima(tree.a[,2], order=c(1, 0 ,1))
    tree.a.n <- tree.a.arima$residuals
    
    #tree.b.arima <-arima(tree.b[,2], order=c(1, 0 ,1))
    #tree.b.n <- tree.b.arima$residuals
    
    tree.a.frame <- data.frame(tree.a$Year, tree.a.n)
    colnames(tree.a.frame) <- c("Year", "A")
    tree.b.frame <- tree.b
    
    tree.a.re <- tree.a.frame$A[tree.a.frame$Year %in% tree.b$Year]
    tree.b.re <- semi_join(tree.b.frame, tree.a.frame, by="Year")
    
    
    tree.sources.frame <- tree.b.re[2:tree.sources.n]
    colnames(tree.sources.frame) <- names(tree.sources[2:tree.sources.n])
    tree.total.frame <- data.frame(tree.a.re, tree.sources.frame)
    colnames(tree.total.frame) <- c("to.test", names(tree.sources.frame))
    
    trees.r2 <- apply(tree.sources.frame, 2, function(x) summary(lm(x~tree.a.re))$r.squared)
    trees.n <- apply(tree.sources.frame, 2, function(x) length(summary(lm(x~tree.a.re))$residuals))
    
    
    turn.to.t <- function(trees.rsquared, trees.residual.n) {
        x.s.lm <- summary(x.lm)
        r.sq <- x.s.lm$r.squared
        just.r <- sqrt(r.sq)
        t <- (just.r*sqrt(length(x.lm$residuals)-2))/sqrt(1-r.sq)
        return(t)
    }
    
    trees.t <- sqrt(trees.r2)*sqrt(trees.n-2)/sqrt(1-trees.r2)
    trees.r <- sqrt(trees.r2)
    
    result.frame <- data.frame(trees.t, trees.r, trees.n)
    colnames(result.frame) <- c("t", "r", "overlap")
    return(format(result.frame, digits=3))
}


treeJackKnife <- function(timemin, timemax,  tree.dataframe, tree.source) {
    
    time <- seq(from=timemin, to=timemax, by=1)
    
    #tree.a <- tree.dataframe[match(time, tree.dataframe$Year, nomatch=0),]
    #tree.b <- tree.source[match(time, tree.source$Year, nomatch=0),]
    
    tree.a <- tree.dataframe[,colSums(is.na(tree.dataframe))<nrow(tree.dataframe)]
    tree.b <- tree.source
  
  #tree.dataframe <- data.frame(tree.dataframe)
  #tree.source <- data.frame(tree.source)
    
    #tree.a <- subset(tree.dataframe, !(tree.dataframe[,1] < timemin | tree.dataframe[,1] > timemax))
    #tree.b <- subset(tree.source, !(tree.source[,1] < timemin | tree.source[,1] > timemax))
    
    tree.a.mod <- tree.a[, colSums(is.na(tree.a)) != nrow(tree.a)]
    
    samp.n <- length(names(tree.a.mod))
    tree.names <- names(tree.a.mod[2:samp.n])
    
    tree.a.re <- tree.a.mod[tree.a.mod$Year %in% tree.b$Year, ]
    tree.b.re <- tree.b[tree.b$Year %in% tree.a.mod$Year, ]
    
    tree.a.re.re <- tree.a.re[2:samp.n]
    
    source <- tree.b.re[,2]
    
    
    n <- length(ls(tree.a.re))
    
    
    group.lm.r2 <- apply(tree.a.re.re, 2, function(x) as.vector(summary(lm(x~source))$r.squared))
    
    group.lm.r <- sqrt(group.lm.r2)
    
    group.lm.res.n <- apply(tree.a.re.re, 2, function(x) as.numeric(length(summary(lm(x~source))$residuals)))
    
    group.t <- sqrt(group.lm.r2)*sqrt(group.lm.res.n-2)/sqrt(1-group.lm.r2)
    
    
    result.frame <- data.frame(group.t, group.lm.r, group.lm.res.n)
    colnames(result.frame) <- c("t-value", "r-value", "Sample Overlap")
    return(format(result.frame, digits=3))
    
    
}


treeJackKnifeHist <- function(timemin, timemax,  tree.dataframe, tree.source) {
    

    tree.a <- tree.dataframe[,colSums(is.na(tree.dataframe))<nrow(tree.dataframe)]
    tree.b <- tree.source
    
    time <- seq(from=timemin, to=timemax, by=1)
    
    tree.a <- tree.dataframe[match(time, tree.dataframe$Year, nomatch=0),]
    tree.b <- tree.source[match(time, tree.source$Year, nomatch=0),]
    
    tree.a <- tree.a[,colSums(is.na(tree.a))<nrow(tree.a)]


    #tree.dataframe <- data.frame(tree.dataframe)
    #tree.source <- data.frame(tree.source)
    
    #tree.a <- subset(tree.dataframe, !(tree.dataframe[,1] < timemin | tree.dataframe[,1] > timemax))
    #tree.b <- subset(tree.source, !(tree.source[,1] < timemin | tree.source[,1] > timemax))
    
    tree.a.mod <- tree.a[, colSums(is.na(tree.a)) != nrow(tree.a)]
    
    samp.n <- length(names(tree.a.mod))
    tree.names <- names(tree.a.mod[2:samp.n])
    
    tree.a.re <- tree.a.mod[tree.a.mod$Year %in% tree.b$Year, ]
    tree.b.re <- tree.b[tree.b$Year %in% tree.a.mod$Year, ]
    
    tree.a.re.re <- tree.a.re[2:samp.n]
    
    source <- tree.b.re[,2]
    
    
    n <- length(ls(tree.a.re))
    
    
    group.lm.r2 <- apply(tree.a.re.re, 2, function(x) as.vector(summary(lm(x~source))$r.squared))
    
    group.lm.r <- sqrt(group.lm.r2)
    
    group.lm.res.n <- apply(tree.a.re.re, 2, function(x) as.numeric(length(summary(lm(x~source))$residuals)))
    
    group.t <- sqrt(group.lm.r2)*sqrt(group.lm.res.n-2)/sqrt(1-group.lm.r2)
    
    
    result.frame <- data.frame(group.t, group.lm.r, group.lm.res.n)
    colnames(result.frame) <- c("t-value", "r-value", "Sample Overlap")
    return(format(result.frame, digits=3))
    
    
}




treeJackKnifeAlt <- function(timemin, timemax,  tree.dataframe) {
    
    
    
    tree.a <- subset(tree.dataframe, !(tree.dataframe[,1] < timemin | tree.dataframe[,1] > timemax))
    
    tree.a.mod <- tree.a[, colSums(is.na(tree.a)) != nrow(tree.a)]
    
    samp.n <- length(names(tree.a.mod))
    tree.names <- names(tree.a.mod[2:samp.n])
    
    tree.a.re.re <- tree.a.mod[2:samp.n]
    
    
    n <- length(tree.a.re.re)
    
    group.lm.r2 <- do.call("rbind", sapply(1:n, FUN = function(i) summary(lm(tree.a.re.re[,i]~as.vector(rowMeans(tree.a.re.re[,-i], na.rm=TRUE))))$r.squared, simplify=FALSE))
    
    group.lm.res.n <- do.call("rbind", sapply(1:n, FUN = function(i) length(summary(lm(tree.a.re.re[,i]~as.vector(rowMeans(tree.a.re.re[,-i], na.rm=TRUE))))$residuals), simplify=FALSE))
    
    
    group.lm.r <- sqrt(group.lm.r2)
    
    group.t <- sqrt(group.lm.r2)*sqrt(group.lm.res.n-2)/sqrt(1-group.lm.r2)
    
    
    result.frame <- data.frame(group.t, group.lm.r, group.lm.res.n)
    colnames(result.frame) <- c("t-value", "r-value", "Sample Overlap")
    return(format(result.frame, digits=3))
    
    
}

treeJackKnifeMultiple <- function(timemin, timemax, tree.dataframe, tree.source.list, return = c("t-value", "r-value", "Sample-Overlap")) {
    
    source.name.list <- sapply(tree.source.list, names)
    source.names <- source.name.list[2,]
    
    all.group.t <- pblapply(tree.source.list, function(tree.source.list) treeJackKnife(timemin, timemax, tree.dataframe, tree.source.list))
    
    t.value <- data.frame(subListExtract(L=all.group.t, name="t-value"))
    colnames(t.value) <- source.names
    t.value$Source <- colnames(t.value)[apply(t.value,1,which.max)]
    
    
    r.value <- data.frame(subListExtract(L=all.group.t, name="r-value"))
    colnames(r.value) <- source.names
    
    samp.over <- data.frame(subListExtract(L=all.group.t, name="Sample Overlap"))
    colnames(samp.over) <- source.names
    
    result <- if (return=="t-value") {
        t.value
    } else if (return=="r-value") {
        r.value
    } else if (return=="Sample-Overlap") {
        samp.over
    }
    
    return(format(result, digits=3))
    
}


treeJackKnifeMultipleSourceSig <- function(timemin, timemax, tree.dataframe, tree.source.list, return = c("t-value", "r-value", "Sample-Overlap")) {
    
    tree.names <- names(tree.dataframe)
    
    cat(gettext(tree.names))
    
    source.name.list <- sapply(tree.source.list, names)
    source.names <- source.name.list[2,]
    
    all.group.t <- pblapply(tree.source.list, function(tree.source.list) treeJackKnife(timemin=timemin, timemax=timemax, tree.dataframe=tree.dataframe, tree.source=tree.source.list))
    
    t.value <- as.data.frame(subListExtract(L=all.group.t, name="t-value"), stringsAsFactors=TRUE)
    colnames(t.value) <- source.names
    n <- length(names(t.value))
    t.value <- as.data.frame(lapply(t.value, as.numeric))
    scaled.t <- t(apply(t.value, 1, function(x) scale(x)[,1]))
    t.value$Mean <- rowMeans(t.value)
    scaled.mean <- rowMeans(scaled.t)
    t.value$SD <- apply(t.value, 1, sd)
    scaled.sd <- apply(scaled.t, 1, sd)
    t.value$SourceValue <- apply(t.value, 1, max)
    scaled.max.value <- apply(scaled.t, 1, max)
    t.value$Source <- colnames(t.value)[apply(t.value,1,which.max)]
    
    
    
    t.value$ZScore <- (scaled.max.value-scaled.mean)/scaled.sd
    t.value$pvalue <- pnorm(-abs(t.value$ZScore))
    
    
    t.value$Difference <- rep("Yes", length(t.value$Mean))
    t.value <- transform(t.value, Difference = ifelse(pvalue > 0.05, "No", Difference))
    t.value.names <- names(t.value)
    
    
    
    r.value <- data.frame(subListExtract(L=all.group.t, name="r-value"))
    colnames(r.value) <- source.names
    
    samp.over <- data.frame(subListExtract(L=all.group.t, name="Sample Overlap"))
    colnames(samp.over) <- source.names
    
    result <- if (return=="t-value") {
        t.value
    } else if (return=="r-value") {
        r.value
    } else if (return=="Sample-Overlap") {
        samp.over
    }
    
    return(format(result, digits=3))
    
    
}



treeJackKnifeMultipleSourceSigHist <- function(timemin, timemax, tree.dataframe, tree.source.list, return = c("t-value", "r-value", "Sample-Overlap")) {
    
    tree.names <- names(tree.dataframe)
    
    cat(gettext(tree.names))
    
    source.name.list <- sapply(tree.source.list, names)
    source.names <- source.name.list[2,]
    
    all.group.t <- pblapply(tree.source.list, function(tree.source.list) treeJackKnifeHist(timemin=timemin, timemax=timemax, tree.dataframe=tree.dataframe, tree.source=tree.source.list))
    
    t.value <- as.data.frame(subListExtract(L=all.group.t, name="t-value"), stringsAsFactors=TRUE)
    colnames(t.value) <- source.names
    n <- length(names(t.value))
    t.value <- as.data.frame(lapply(t.value, as.numeric))
    scaled.t <- t(apply(t.value, 1, function(x) scale(x)[,1]))
    t.value$Mean <- rowMeans(t.value)
    scaled.mean <- rowMeans(scaled.t)
    t.value$SD <- apply(t.value, 1, sd)
    scaled.sd <- apply(scaled.t, 1, sd)
    t.value$SourceValue <- apply(t.value, 1, max)
    scaled.max.value <- apply(scaled.t, 1, max)
    t.value$Source <- colnames(t.value)[apply(t.value,1,which.max)]
    
    
    
    t.value$ZScore <- (scaled.max.value-scaled.mean)/scaled.sd
    t.value$pvalue <- pnorm(-abs(t.value$ZScore))
    
    
    t.value$Difference <- rep("Yes", length(t.value$Mean))
    t.value <- transform(t.value, Difference = ifelse(pvalue > 0.05, "No", Difference))
    t.value.names <- names(t.value)
    
    
    
    r.value <- data.frame(subListExtract(L=all.group.t, name="r-value"))
    colnames(r.value) <- source.names
    
    samp.over <- data.frame(subListExtract(L=all.group.t, name="Sample Overlap"))
    colnames(samp.over) <- source.names
    
    result <- if (return=="t-value") {
        t.value
    } else if (return=="r-value") {
        r.value
    } else if (return=="Sample-Overlap") {
        samp.over
    }
    
    return(format(result, digits=3))
    
    
}



treeJackKnifeMultipleSourceSigSamp <- function(timemin, timemax, tree.dataframe, tree.source.list, return = c("t-value", "r-value", "Sample-Overlap")) {
    
    tree.n <- length(tree.dataframe)
    tree.names <- names(tree.dataframe[2:tree.n])
    
    source.name.list <- sapply(tree.source.list, names)
    source.names <- source.name.list[2,]
    
    all.group.t <- pblapply(tree.source.list, function(tree.source.list) treeJackKnife(timemin=timemin, timemax=timemax, tree.dataframe=tree.dataframe, tree.source=tree.source.list))
    
    t.value <- as.data.frame(subListExtract(L=all.group.t, name="t-value"), stringsAsFactors=TRUE)
    colnames(t.value) <- source.names
    n <- length(names(t.value))
    t.value <- as.data.frame(lapply(t.value, as.numeric))
    scaled.t <- t(apply(t.value, 1, function(x) scale(x)[,1]))
    t.value$Mean <- rowMeans(t.value)
    scaled.mean <- rowMeans(scaled.t)
    t.value$SD <- apply(t.value, 1, sd)
    scaled.sd <- apply(scaled.t, 1, sd)
    t.value$SourceValue <- apply(t.value, 1, max)
    scaled.max.value <- apply(scaled.t, 1, max)
    t.value$Source <- colnames(t.value)[apply(t.value,1,which.max)]
    
    
    
    t.value$ZScore <- (scaled.max.value-scaled.mean)/scaled.sd
    t.value$pvalue <- pnorm(-abs(t.value$ZScore))
    
    t.value$Difference <- rep("Yes", length(t.value$Mean))
    t.value <- transform(t.value, Difference = ifelse(pvalue > 0.05, "No", Difference))
    t.value.names <- names(t.value)
    t.value <- data.frame(tree.names, t.value)
    colnames(t.value) <- c("Specimen", t.value.names)
    
    
    r.value <- data.frame(subListExtract(L=all.group.t, name="r-value"))
    colnames(r.value) <- source.names
    
    
    samp.over <- data.frame(subListExtract(L=all.group.t, name="Sample Overlap"))
    colnames(samp.over) <- source.names
    
    result <- if (return=="t-value") {
        t.value
    } else if (return=="r-value") {
        r.value
    } else if (return=="Sample-Overlap") {
        samp.over
    }
    
    return(format(result, digits=3))
    
    
}

percent <- function(x, digits = 2, format = "f", ...) {
    paste0(formatC(100 * x, format = format, digits = digits, ...), "%")
}

populationSig <- function(x) {
    Yes <- length(subset(x$pvalue, x$Difference=="Yes"))
    No <- length(subset(x$pvalue, x$Difference=="No"))
    Ratio <- Yes/(sum(Yes, No))
    
    result.frame <- data.frame(Yes, No, Ratio)
    return(format(result.frame, digits=3))
    
}

populationSigDefinition <- function(x, source.hypothesis) {
    Yes <- length(subset(x$pvalue, x$Difference=="Yes" & x$Source == source.hypothesis))
    No <- length(x$pvalue)-Yes
    Ratio <- Yes/length(x$pvalue)
    
    result.frame <- data.frame(Yes, No, Ratio)
    return(format(result.frame, digits=3))
    
}

populationSigDefinitionMultiple <- function(x, source.hypotheses) {
    
    hold <- rep(0, length(source.hypotheses))
    hold.frame <- data.frame(t(hold))
    colnames(hold.frame) <- source.hypotheses
    
    n <- length(x$Difference)
    
    x.subset <- subset(x, x$Difference=="Yes")
    x.source <- table(x.subset$Source)
    x.values <- as.numeric(paste(x.source))
    x.frame <- data.frame(t(x.values))
    colnames(x.frame) <- names(x.source)
    Yes <- length(subset(x$Source, x$Difference=="Yes"))
    None.n <- length(subset(x$Source, x$Difference=="No"))
    None <- None.n/n
    Ratio <- Yes/(Yes+None.n)
    
    results.frame <- merge(hold.frame, x.frame, all=TRUE)
    results.frame[is.na(results.frame)] <- 0
    
    
    result.frame <- data.frame(results.frame[2,]/n, None, Ratio)
    return(format(result.frame, digits=2))
    
}

multiplePopulationSig <- function(t.table.list) {
    
    yes.no.table <- as.data.frame(sapply(t.table.list, FUN=populationSig, USE.NAMES=TRUE))
    return(yes.no.table)
}




multiplePopulationSig <- function(x,...) {
    
    yes.no.table <- as.data.frame(sapply(t.table.list, FUN=populationSig, USE.NAMES=TRUE))
    return(yes.no.table)
}

sourceGrid <- function(fileDirectory, sourceDefinitions, timemin, timemax, source.hypotheses) {
    setwd(fileDirectory)
    temp <- list.files(pattern="*.txt")
    myRWL <- pblapply(temp, readRWLArima)
    myJackKnife <- lapply(myRWL,  function(x) treeJackKnifeMultipleSourceSigHist(timemin, timemax, x, source.list, return="t-value"))
    mySigDef <- lapply(myJackKnife, function(x) populationSigDefinitionMultiple(x, source.hypotheses))
    myResults <- ldply (mySigDef, data.frame)
    rownames(myResults) <- gsub(".txt", "", temp, )
    return(myResults)
}

results <- read.csv("data/myResults.csv", sep=",")
