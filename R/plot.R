#' Fit a loess model to the p-value data for later use in generating heatmaps.
#'
#' Given a vector of p-values and two population parameters, this will fit a
#' \code{\link{loess}} model to the data, with the model of the form z ~ x * y,
#' where z is the predicted p-value.  This can later be used to make a heatmap,
#' and is separated from the \code{\link{heatmap2}} function because of the
#' long amount of time it can take for this function to run.
#' 
#' 
#' @param pval vector of p-values from simulated runs
#' @param paramX population parameter to condition \code{z} on (e.g., a column
#' from a matrix created using \code{\link{createParamMx}}).  This will appear
#' on the x-axis of the heatmap.
#' @param paramY same as \code{paramX}, but for the Y-axis
#' @param minx minimum x value
#' @param maxx maximum x value
#' @param miny minimum y value
#' @param maxy maximum y value
#' @param svec a logical vector to select out only certain rows from
#' \code{pval}, \code{paramX} and \code{paramY}
#' @param span a value to be passed along to the \code{\link{loess}} function
#' that determines the smoothing level
#' @return A fitted \code{\link{loess}} model.
#' @seealso \code{\link{loess}}, \code{\link{heatmap2}}
#' @examples
#' 
#' pmx <- createParamMx(100)
#' 
#' result <- apply(pmx, 1, function(x) {fitanova(mkDf(24,24,mcr.params=x,FALSE),FALSE)})
#' lpred <- loessPred(pval=result["pmax",],
#'                 paramX=pmx[,"t11"], paramY=pmx[,"w11"])
#' 
#' heatmap2(lpred, "t11", "w11", "F1+F2")
#' 
#' @export loessPred
loessPred <- function(pval, paramX, paramY,
                      minx=0.01,maxx=2.99,miny=.01,maxy=2.99,
                      svec=rep(TRUE, length(paramX)),span=.9) {
                                        # calculate predictions (Type I error or Power)
                                        # get rid of NAs, and use svec to select out observations
  lvec <- (!is.na(pval)) & svec
  ff <- pval[lvec]
  p.sig <- ifelse(ff<=.05,1,0)
  x <- paramX[lvec]; y <- paramY[lvec]
  p.map <- loess(p.sig ~ x * y, span=span)
  return(p.map)
}



#' Create a heatmap from a given loess fit.
#' 
#' Creates a graph of estimated p-values from a loess fit to stored simulation
#' data.  The loess fit should be created using \code{\link{loessPred}}.
#' 
#' 
#' @param lofit fitted loess object
#' @param labelX label for the x-axis variable.  This can be a character string
#' or an expression.  If one of "t00", "t11", "w00", "w11", will be properly
#' formatted using \code{\link{plotmath}}.
#' @param labelY same as \code{labelX}, but for the y axis.
#' @param title title for the plot
#' @param mxseq sequence for the x variable, determining points at which
#' predictions are calculated
#' @param myseq same as \code{mxseq}, but for the y variable
#' @param powermap if plotting power (TRUE) or type I error (FALSE); allows for
#' different grayscale sensitivities.
#' @seealso \code{\link{loessPred}}
#' @examples
#' 
#' pmx <- createParamMx(100)
#' 
#' result <- apply(pmx, 1, function(x) {fitanova(mkDf(24,24,mcr.params=x,FALSE),FALSE)})
#' lpred <- loessPred(pval=result["pmax",],
#'                 paramX=pmx[,"t11"], paramY=pmx[,"w11"])
#' 
#' heatmap2(lpred, "t11", "w11", "F1+F2")
#' 
#' @export heatmap2
heatmap2 <- function(lofit, labelX, labelY, title,
                     mxseq=seq(.01,2.99,length.out=300),
                     myseq=seq(.01,2.99,length.out=300),
                     powermap=FALSE) {
  p.grid <- expand.grid(x=mxseq,y=myseq)
  lpred.mx <- predict(lofit, newdata=p.grid)
  rownames(lpred.mx) <- mxseq
  colnames(lpred.mx) <- myseq
                                        # color scheme for heatmap
  makergb <- function(x,ming,maxg) {
    rge.x <- max(x)-min(x)
    px <- (x-min(x))/rge.x
    rge.g <- abs(ming-maxg)
    pg <- ming + px*rge.g*(ifelse(ming>maxg,-1,1))
    return(unlist(lapply(pg, function(y) {rgb(y,y,y)})))
  }
  if (!powermap) {
    mycol1 <- makergb(seq(.05,.09,.01), .95, .70)
    mycol2 <- makergb(seq(.10,.20,.01), .70, .50)
#    mycol2 <- makergb(seq(.10,.20,.01), .70, .40)    
    mycol3 <- makergb(seq(.21,.70,.01), .50, .25)
#    mycol3 <- makergb(seq(.21,.70,.01), .40, .15)
    mycol <- c(mycol1,mycol2,mycol3)
    mybreaks <- c(seq(.05,.70,.01),2.0)
  } else {
    mycol1 <- makergb(seq(.05,.09,.01), .95, .90)
    mycol2 <- makergb(seq(.10,.20,.01), .90, .70)
    mycol3 <- makergb(seq(.21,.70,.01), .70, .25)
    mycol <- c(mycol1,mycol2,mycol3)
    mybreaks <- c(seq(.05,.70,.01),2.0)
  }
  lab2expr <- function(x) { # map label to expression
    ee <- x
    if (!is.expression(ee)) {
      if (x=="t00") {
        ee <- expression({tau[0][0]}^2)
      } else if (x == "t11") {
        ee <- expression({tau[1][1]}^2)
      } else if (x == "w00") {
        ee <- expression({omega[0][0]}^2)
      } else if (x == "w11") {
        ee <- expression({omega[1][1]}^2)
      } else if (x == "subj_logratio") {
        ee <- expression(paste("log (", {tau[1][1]}/{tau[0][0]}, ")", sep=""))
      } else if (x == "item_logratio") {
        ee <- expression(paste("log (", {omega[1][1]}/{omega[0][0]}, ")", sep=""))
      } else if (x == "subj_cr") {
        ee <- expression(paste("log (", {tau[1][1]}/{omega[0][0]}, ")", sep=""))
      } else if (x == "item_cr") {
        ee <- expression(paste("log (", {omega[1][1]}/{tau[0][0]}, ")", sep=""))
      }
    } else {}
    return(ee)
  } 
  e1 <- lab2expr(labelX); e2 <- lab2expr(labelY) 
  
  mxseq <- as.numeric(rownames(lpred.mx))
  myseq <- as.numeric(colnames(lpred.mx))
  minx <- floor(min(mxseq)); maxx <- ceiling(max(mxseq))
  miny <- floor(min(myseq)); maxy <- ceiling(max(myseq))
  image(z=lpred.mx,
        x=mxseq, y=myseq,
        xlab=e1,ylab=e2,
        breaks=mybreaks, col=mycol,
        main=title, cex.lab=1.4,
        xaxt='n',yaxt='n',useRaster=TRUE)
  
  contour(lpred.mx, x=mxseq,y=myseq, add=TRUE, levels=c(seq(0,.2,by=.01),seq(.25,.5,by=.05),seq(.6,1,.1)))
  axis(1,at=round(seq(minx,maxx,length.out=7),2))
  axis(2,at=round(seq(miny,maxy,length.out=7),2))  
}
