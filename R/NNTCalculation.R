#' @title Calculate the NNT based on the Kaplan-Meier estimated
#' survival rates between the treatment and control groups
#'
#' @description For survival endpoints, the NNT-KM is computed as the reciprocal of
#'  the absolute risk reduction (ARR),
#'  which is the difference in Kaplan-Meier estimated survival rates or
#'  the difference in cumulative incidences at a time point of clinical
#'  interest between the treatment and control groups.
#'
#' @param time The time to event or censor.
#' @param status The indicator of the event or censor at the end of the follow-up.
#' @param arm The variable indicates the treatment (arm = 1) and control (arm = 0) groups.
#' @param tau The chosen time point of clinical interest.
#' @param confint The percentile of confidence interval. The default value is
#'   \code{confint = 0.95}.
#' @param digits The decimal of the results. The default value is \code{digits = 3}.
#'
#' @return A matrix contains the KM-NNT and its confidence interval.
#'
#' @examples
#' library(survival)
#' dat <- pbc[!is.na(pbc$trt),]
#' time <- dat$time/365.25
#' status <- (dat$status == 2) + 0
#' arm <- (dat$trt == 2) + 0
#' KM2NNT(time, status, arm, tau = NULL, confint = 0.95, digits = 3)
#'
#' @import "survival","stats"
#'
#' @references
#' 1. Altman DG, Andersen PK: Calculating the number needed to treat for trials
#' where the outcome is time to an event. BMJ 319:1492-5, 1999
#'
#' 2. Altman DG: Confidence intervals for the number needed to treat. BMJ
#' 317:1309-12, 1998
#'
#' @export
#'

KM2NNT <- function(time, status, arm, tau = NULL, confint = 0.95, digits = 3) {
    idx <- arm == 0; tt <- time[idx]; tau0max <- max(tt)
    idx <- arm == 1; tt <- time[idx]; tau1max <- max(tt)
    tauMax = round(min(tau0max, tau1max),digits)
    if(!is.null(tau)){
        if(tau > tauMax){
            tau <- tauMax
        }
    }
    if(is.null(tau)){
        tau <- tauMax
    }
    alpha <- 1 - confint
    KM.fit <- survfit(Surv(time, status) ~ arm, data = NULL)
    KMsurv <- summary(KM.fit, censor = TRUE, times = tau)
    KMres <- data.frame(
        Time = KMsurv$time,
        Surv = round(KMsurv$surv*100,digits),
        SE = round(KMsurv$std.err*100,digits),
        Surv95LCI = round(KMsurv$lower*100,digits),
        Surv95UCI = round(KMsurv$upper*100,digits),
        Group = KMsurv$strata
    )
    KMres0 <- KMres[KMres$Group == "arm=0",]
    names(KMres0) <- paste(names(KMres0), "0", sep = "")
    KMres1 <- KMres[KMres$Group == "arm=1",]
    names(KMres1) <- paste(names(KMres1), "U", sep = "")
    KMResU <- cbind(KMres0, KMres1)
    KMResU$ARRU <- round((KMResU$SurvU - KMResU$Surv0)/100,digits)
    KMResU$seSRRU <- round(sqrt((KMResU$SEU/100)^2 + (KMResU$SE0/100)^2), digits)
    KMResU$NNTU <- 1/KMResU$ARRU
    KMResU$NNTULower <- 1/(KMResU$SurvU/100 - KMResU$Surv0/100 -
                               qnorm(1 - alpha/2)*KMResU$seSRRU)
    KMResU$NNTUUpper <- 1/(KMResU$SurvU/100 - KMResU$Surv0/100 +
                               qnorm(1 - alpha/2)*KMResU$seSRRU)
    Res <- matrix(NA, nrow = 1, ncol = 6)
    Res[1,1] <- as.numeric(KMResU["Time0"])
    Res[1,2] <- round(as.numeric(KMResU["ARRU"]),digits)
    Res[1,3] <- round(as.numeric(KMResU["ARRU"]) - qnorm(1 - alpha/2)*as.numeric(KMResU["seSRRU"]), digits)
    Res[1,4] <- round(as.numeric(KMResU["ARRU"]) + qnorm(1 - alpha/2)*as.numeric(KMResU["seSRRU"]), digits)
    Res[1,5] <- round(1/as.numeric(KMResU["ARRU"]), digits)
    if (KMResU$NNTULower <= 0 && KMResU$NNTUUpper > 0) {
        Res[1,6] <- paste(round(abs(1/round(Res[1,4],digits)), digits),
                          " to Inf to ",
                          abs(round(1/round(Res[1,3],digits), digits)) ,sep = "")
    } else {
        Res[1,6] <- paste(round(abs(1/round(Res[1,4],digits)), digits),
                          " to ",
                          round(abs(1/round(Res[1,3],digits)), digits) ,sep = "")
    }
    colnames(Res) <- c("TimePoint", "ARR",
                       paste("ARR(",confint*100,"%LowerCI)", sep = ""),
                       paste("ARR(",confint*100,"%UpperCI)", sep = ""),
                       "NNT", "NNT(Benefit) -> NNT(Harm)")
    row.names(Res) <- "KM(Arm1-Arm0)"
    return(Res)
}



#' @title Calculate the NNT based on the restricted mean survival times between
#'  the treatment and control groups
#'
#' @description For survival endpoints, the NNT-RMST is defined as the RMST in
#'   the control group divided by the difference in RMSTs between the treatment
#'   and control groups up to a chosen time \code{t}.
#'
#' @param time The time to event or censor.
#' @param status The indicator of the event or censor at the end of the follow-up.
#' @param arm The variable indicates the treatment (arm = 1) and control (arm = 0) groups.
#' @param tau The chosen time point of clinical interest.
#' @param confint The percentile of confidence interval. The default value is
#'   \code{confint = 0.95}.
#' @param digits The decimal of the results. The default value is \code{digits = 3}.
#'
#' @return A matrix contains the RMST-NNT and its confidence interval.
#'
#' @examples
#' library(survival)
#' dat <- pbc[!is.na(pbc$trt),]
#' time <- dat$time/365.25
#' status <- (dat$status == 2) + 0
#' arm <- (dat$trt == 2) + 0
#' RM2NNT(time, status, arm, tau = NULL, confint = 0.95, digits = 3)
#'
#' @import "survRM2","stats"
#'
#' @references
#' 1. An alternative approach for estimating the number needed to treat for survival endpoints.
#'  PLoS One. 2019 Oct 18;14(10):e0223301. doi: 10.1371/journal.pone.0223301.
#'
#' @export
#'
RM2NNT <- function(time, status, arm, tau = NULL, confint = 0.95, digits = 3) {
    idx <- arm == 0; tt <- time[idx]; tau0max <- max(tt)
    idx <- arm == 1; tt <- time[idx]; tau1max <- max(tt)
    tauMax = min(tau0max, tau1max)
    if(!is.null(tau)){
        if(tau > tauMax){
            tau <- tauMax
        }
    }
    if(is.null(tau)){
        tau <- tauMax
    }
    RMST1 <- survRM2::rmst2(time = time, status = status,
                            arm = arm, tau = tau,
                            alpha = 1 - confint)
    DeathTimeRMST <- RMST1$RMST.arm0$rmst[1]
    RMST <- matrix(NA, nrow = 2, ncol = 6)
    RMST[1,1] <- round(tau,digits)
    RMST[1,2] <- sum(RMST1$RMST.arm1$fit$n.event)
    RMST[1,3:6] <- RMST1$RMST.arm1$rmst
    RMST[2,1] <- round(tau,digits)
    RMST[2,2] <- sum(RMST1$RMST.arm0$fit$n.event)
    RMST[2,3:6] <- RMST1$RMST.arm0$rmst
    row.names(RMST) <- c("Arm1", "Arm0")
    colnames(RMST) <- c("TimePoint", "N", "RMST", "SE(RMST)", "RMSTLowerCI", "RMSTUpperCI")
    rmstNNT <- matrix(NA, nrow = 1, ncol = 7)
    colnames(rmstNNT) <- c("TimePoint", "RMST(1-0)", "NNT",
                           "LCI(NNT)", "UCI(NNT)", "NNT(Benefit) -> NNT(Harm)",
                           "1-DeathTime(RMST0)")
    rmstNNT[1,1] <- round(tau,digits)
    rmstNNT[1,2] <- round(RMST1$unadjusted.result[1,1],digits)
    rmstNNT[1,3] <- round(DeathTimeRMST/RMST1$unadjusted.result[1,1],digits)
    rmstNNT[1,4] <- 1/round(RMST1$unadjusted.result[2,3] - 1,digits)
    rmstNNT[1,5] <- 1/round(RMST1$unadjusted.result[2,2] - 1,digits)
    if (rmstNNT[1,5] <= 0 && rmstNNT[1,4] > 0) {
        rmstNNT[1,6] <- paste(formatC(abs(rmstNNT[1,4]), format = "f", digits = digits),
                              " to Inf to ",
                              formatC(abs(rmstNNT[1,5]), format = "f", digits = digits)
                              , sep = "")
    } else {
        rmstNNT[1,6] <- paste(formatC(abs(rmstNNT[1,4]), format = "f", digits = digits),
                              " to ",
                              formatC(abs(rmstNNT[1,5]), format = "f", digits = digits)
                              , sep = "")
    }
    rmstNNT[1,7] <- round(DeathTimeRMST, digits)
    Res <- list(RMST = RMST, rmstNNT = rmstNNT)
    return(Res)
}


#' @title Compare the performance between the NNT-RMST and NNT-KM through
#'   the average life gain per patient
#'
#' @description For the NNT-RMST, the average life gain per patient is the area
#' between the survival curves, which is the instrinsic treatment
#'  benefit in survival time during the t-period follow-up. For the NNT-KM, the
#'  average life gain per patient is defined as the ratio between the average
#'  survival time of one death in patients and the NNT-KM up to t.
#'
#' @param time The time to event or censor.
#' @param status The indicator of the event or censor at the end of the follow-up.
#' @param arm The variable indicates the treatment (arm = 1) and control (arm = 0) groups.
#' @param tau The chosen time point of clinical interest.
#' @param confint The percentile of confidence interval. The default value is
#'   \code{confint = 0.95}.
#' @param digits The decimal of the results. The default value is \code{digits = 3}.
#'
#' @return A list contains:
#' @return \item{RMSTNNT}{The RMST-NNT and its confidence interval.}
#' @return \item{KMNNT}{The KM-NNT and its confidence interval.}
#' @return \item{LifeGain}{The average life gain per patient based on the RMST-NNT and KM-NNT.}
#'
#' @examples
#' library(survival)
#' dat <- pbc[!is.na(pbc$trt),]
#' time <- dat$time/365.25
#' status <- (dat$status == 2) + 0
#' arm <- (dat$trt == 2) + 0
#' RMvsKM(time, status, arm, tau = NULL, confint = 0.95, digits = 3)
#'
#' @references
#' 1. An alternative approach for estimating the number needed to treat for survival endpoints.
#'  PLoS One. 2019 Oct 18;14(10):e0223301. doi: 10.1371/journal.pone.0223301.
#'
#' @export
#'
RMvsKM <- function(time, status, arm, tau = NULL, confint = 0.95, digits = 3) {
    idx <- arm == 0; tt <- time[idx]; tau0max <- max(tt)
    idx <- arm == 1; tt <- time[idx]; tau1max <- max(tt)
    tauMax = min(tau0max, tau1max)
    if(!is.null(tau)){
        if(tau > tauMax){
            tau <- tauMax
        }
    }
    if(is.null(tau)){
        tau <- tauMax
    }
    RMNNT <- RM2NNT(time, status, arm, tau, confint, digits)
    KMNNT <- KM2NNT(time, status, arm, tau, confint, digits)
    lifeGain <- matrix(NA, nrow = 1, ncol = 5)
    lifeGain[1,1] <- round(tau,digits)
    lifeGain[1,2] <- round(as.numeric(RMNNT$rmstNNT[,"NNT"]), digits = digits)
    lifeGain[1,3] <- as.numeric(RMNNT$rmstNNT[,"RMST(1-0)"])
    lifeGain[1,4] <- round(as.numeric(KMNNT[,"NNT"]), digits)
    lifeGain[1,5] <- round(round(RMNNT$RMST["Arm0","RMST"], digits)/as.numeric(KMNNT[,"NNT"]), digits = digits)
    colnames(lifeGain) <- c("TimePoint", "NNT-RMST", "LifeGain/Patient(RMST)",
                            "NNT-KM", "LifeGain/Patient(KM)")
    Z <- list()
    Z$Note <- "NNT-RMST vs. NNT-KM"
    Z$RMSTNNT <- RMNNT$rmstNNT
    Z$KMNNT <- KMNNT
    Z$LifeGain <- lifeGain
    return(Z)
}


