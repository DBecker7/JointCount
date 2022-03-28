library(tictoc)
library(lubridate)
library(ggplot2)
library(tidyr)
library(rjags)
library(splines)
library(dplyr)

select <- dplyr::select
filter <- dplyr::filter


bcoutline <- readRDS("Data/Modified/bcoutline2.rds")
FireCentreShape <- readRDS("Data/Modified/shape3.rds")
#bc3 <- readRDS("Data/Modified/bc3.rds")


load.module("glm")
load.module("dic")

ql <- function(x) quantile(x, probs = 0.055, na.rm = TRUE)
qh <- function(x) quantile(x, probs = 0.945, na.rm = TRUE)

# Sets global variables - bad, but useful
set_timer <- function(maxiter = NULL, units = "mins"){
    .maxiter <<- maxiter
    .units <<- units
    .time_at_start <<- Sys.time()
    .last_time <<- .time_at_start
    .record_of_times <<- c()
}

get_timer <- function(current_loop){
    elapsed <- difftime(Sys.time(), .last_time, units = .units)
    .record_of_times <<- c(.record_of_times, elapsed)
    if(!is.null(.maxiter)){
        est_remain <- mean(.record_of_times)*(.maxiter - current_loop)
        cat(paste0("Loop ", current_loop, 
            " of ", .maxiter, 
            " took ", round(elapsed, 3), 
            " ", .units, 
            ", estimated ", round(est_remain, 3), " ", .units,
            " remaining.\n"))
    } else {
        cat(paste0("Loop ", current_loop, 
            " took ", round(elapsed, 3), 
            " ", units, ".\n"))
    }
    .last_time <<- Sys.time()
}

end_timer <- function(clean = TRUE){
    total_time <- difftime(Sys.time(), .time_at_start, units = .units)
	if(length(.record_of_times) > 1){
		cat(paste0("Done. ", length(.record_of_times), " loops took ",
			round(total_time, 3), " ", .units, 
			", average loop took ", round(mean(.record_of_times), 3), 
			" ", .units, ".\n"))
	} else {
		cat(paste0("Done. ", round(total_time, 3), " ", .units,
			" elapsed.\n"))
	}
    if(clean) rm(list = c(".time_at_start", ".maxiter", ".record_of_times",
            ".last_time", ".units"), 
        envir = .GlobalEnv)
}

# Usage:
#set_timer(5); for(i in 1:5){
#    Sys.sleep(0.5)
#    get_timer(i, "secs")
#}; end_timer()



melt.mcmc.list.rename <- function(mcmclist){
    mcmcdf <- as.data.frame(mcmclist[[1]])
    if(length(mcmclist) > 1){
        for(i in 2:length(mcmclist)){
            mcmcdf <- rbind(mcmcdf, 
                as.data.frame(mcmclist[[i]]))
        }
    }
    mcmcdf$chain = factor(rep(1:length(mcmclist), 
        each = nrow(mcmclist[[1]])))
    atts <- attr(mcmclist[[1]], "mcpar")
    
    n.iter <- atts[2] - atts[1] + atts[3]
    n.thin <- atts[3]
    iter <- 1:nrow(mcmclist[[1]])
    
    mcmcdf$iter <- rep(iter*n.thin, length(mcmclist))
    
    names(mcmcdf) <- gsub("\\[", "_", names(mcmcdf))
    names(mcmcdf) <- gsub("]", "", names(mcmcdf))
    names(mcmcdf) <- gsub(",", ".", names(mcmcdf))
    
    mcmcdf
}


mix_data_prep <- function(fires, units = "metric", 
    mycov = c("ISI", "FWI", "Slope", "WIND", "TEMP", "RH", "Elevation")){
    
    fire.season <<- yday("2003-04-01"):yday("2003-10-01")
    fires <- filter(fires, DayYear %in% fire.season, !is.na(FireCentre2))
    
    
    Ni <- fires %>% 
        group_by(DayYear, FireCentre2) %>% 
        summarise(.groups = "drop", Ni = n()) %>% #, 
        #FWI = mean(FWI, na.rm = TRUE), 
        #WIND = mean(WIND, na.rm = TRUE),
        #TEMP = mean(TEMP, na.rm = TRUE),
        #RH = mean(RH, na.rm = TRUE))%>% 
        right_join(expand.grid(DayYear = 1:365, FireCentre2 = 1:6),
            by = c("DayYear", "FireCentre2")) %>% 
        mutate(Ni = ifelse(is.na(Ni), 0, Ni)) %>%
        filter(DayYear %in% fire.season) %>% 
        as.data.frame()
    
    
    B <- bs(fire.season, knots = quantile(fire.season, 0:10/10))
    J <- ncol(B)
    A <- matrix(runif(J^2)*2-1, ncol=J) 
    Sigma <- t(A) %*% A
    
    ic_switch <- if(units == "metric"){
        matrix(ncol = 3, byrow = TRUE,
            data = c(0, 0, 0.1,
                0.1, 0, 0.25,
                0.5, 0.25, 0.75,
                1, 0.75, 1.25,
                1.5, 1.25, 1.75, 
                2, 1.75, 2.25))
    } else {
        matrix(ncol = 3, byrow = TRUE,
            data = c(0, 0, 0.05,
                0.1, 0, 0.25,
                0.4, 0.2, 0.6,
                0.8, 0.6, 1,
                1.2, 1, 1.4, 
                1.6, 1.4, 1.8))
    }
    
    
    oldx <- matrix(data = rep(fires$FinalControlSize, 2), ncol = 2, byrow = FALSE)
    newx <- oldx
    
    for(i in 1:nrow(newx)){
        if(newx[i, 1] %in% ic_switch[, 1]){
            this <- which(ic_switch[, 1] == newx[i, 1])
            newx[i, ] <- ic_switch[this, 2:3]
        } else {
            newx[i, ] <- newx[i, ] + c(-0.05, 0.05)
        }
        
        oldx[i, ] <- oldx[i, ] + c(-0.05,0.05)
        if(oldx[i,1] < 0){
            oldx[i,1] <- 0
        }
    }
    
    xcovars <- fires[, mycov]
    xcovars <- apply(fires[, mycov], 2, function(x){
        (x - mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE)
    })
    
    
    jdata <- list(
        ndays = length(fire.season),
        idays = nrow(Ni),
        Ni = Ni$Ni,
        jfires = nrow(fires),
        Xijr1 = rep(1, nrow(fires)),
        nxcovar = length(mycov),
        xcovars = xcovars,
        #nNcovar = 4,
        DayYear = fires$DayYear - min(fire.season) + 1,
        DayYearN = Ni$DayYear - min(fire.season) + 1,
        FireCentre2 = fires$FireCentre2,
        FireCentre2N = Ni$FireCentre2,
        z = (Ni$Ni > 0)*1,
        ones = rep(1, nrow(Ni)),
        B = B, J = ncol(B),
        lambeta_posdef = Sigma,
        zero_vec = rep(0, J),
        dint = cbind(oldx, newx)
    )
    
    idata <- list(
        Xijr_star = fires$FinalControlSize + 0.0001
    )
    
    mycovnames <- mycov
    names(mycovnames) <- 1:length(mycovnames)
    
    return(list(jdata = jdata, idata = idata, mycov = mycovnames))
}

round_data_prep <-  function(fires, units = "metric", 
    mycov = c("ISI", "FWI", "Slope", "WIND", "TEMP", "RH", "Elevation")){

	mix_data <- mix_data_prep(fires = fires, units = units, mycov = mycov)
	mix_data$jdata$dint <- mix_data$jdata$dint[,1:2]
	mix_data
}


mix_data_prep_2 <- function(fires, ncovar,
    units = "metric", 
    mycov = c("ISI", "FWI", "Slope", "WIND", "TEMP", "RH", "Elevation"),
    myncov = c("meanPrecip", "meanTemp")){
    
    fire.season <<- yday("2003-04-01"):yday("2003-10-01")
    fires <- filter(fires, DayYear %in% fire.season, !is.na(FireCentre2))
    
    ncovar <- filter(ncovar, Year == fires$FireYear[1])
    
    Ni <- fires %>% 
        group_by(DayYear, FireCentre2) %>% 
        summarise(.groups = "drop", Ni = n()) %>% #, 
        #FWI = mean(FWI, na.rm = TRUE), 
        #WIND = mean(WIND, na.rm = TRUE),
        #TEMP = mean(TEMP, na.rm = TRUE),
        #RH = mean(RH, na.rm = TRUE))%>% 
        right_join(expand.grid(DayYear = 1:365, FireCentre2 = 1:6),
            by = c("DayYear", "FireCentre2")) %>% 
        mutate(Ni = ifelse(is.na(Ni), 0, Ni)) %>%
        filter(DayYear %in% fire.season) %>% 
        as.data.frame()
    
    ncovar <- ncovar %>%
        select(all_of(myncov), FireCentre2, DayYear) %>%
        right_join(Ni, by = c("FireCentre2", "DayYear")) %>% 
        select(all_of(myncov))
    
    
    B <- bs(fire.season, knots = quantile(fire.season, 0:10/10))
    J <- ncol(B)
    A <- matrix(runif(J^2)*2-1, ncol=J) 
    Sigma <- t(A) %*% A
    
    ic_switch <- if(units == "metric"){
        matrix(ncol = 3, byrow = TRUE,
            data = c(0, 0, 0.1,
                0.1, 0, 0.25,
                0.5, 0.25, 0.75,
                1, 0.75, 1.25,
                1.5, 1.25, 1.75, 
                2, 1.75, 2.25))
    } else {
        matrix(ncol = 3, byrow = TRUE,
            data = c(0, 0, 0.05,
                0.1, 0, 0.25,
                0.4, 0.2, 0.6,
                0.8, 0.6, 1,
                1.2, 1, 1.4, 
                1.6, 1.4, 1.8))
    }
    
    
    oldx <- matrix(data = rep(fires$FinalControlSize, 2), ncol = 2, byrow = FALSE)
    newx <- oldx
    
    for(i in 1:nrow(newx)){
        if(newx[i, 1] %in% ic_switch[, 1]){
            this <- which(ic_switch[, 1] == newx[i, 1])
            newx[i, ] <- ic_switch[this, 2:3]
        } else {
            newx[i, ] <- newx[i, ] + c(-0.05, 0.05)
        }
        
        oldx[i, ] <- oldx[i, ] + c(-0.05,0.05)
        if(oldx[i,1] < 0){
            oldx[i,1] <- 0
        }
    }
    
    xcovars <- fires[, mycov, drop = FALSE]
    xcovars <- apply(fires[, mycov, drop = FALSE], 2, function(x){
        (x - mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE)
    })
    
    
    jdata <- list(
        ndays = length(fire.season),
        idays = nrow(Ni),
        Ni = Ni$Ni,
        jfires = nrow(fires),
        Xijr1 = rep(1, nrow(fires)),
        nxcovar = length(mycov),
        xcovars = xcovars,
        ncovars = as.matrix(ncovar),
        nNcovar = ncol(ncovar),
        DayYear = fires$DayYear - min(fire.season) + 1,
        DayYearN = Ni$DayYear - min(fire.season) + 1,
        FireCentre2 = fires$FireCentre2,
        FireCentre2N = Ni$FireCentre2,
        z = (Ni$Ni > 0)*1,
        ones = rep(1, nrow(Ni)),
        B = B, J = ncol(B),
        lambeta_posdef = Sigma,
        zero_vec = rep(0, J),
        dint = cbind(oldx, newx)
    )
    
    idata <- list(
        Xijr_star = fires$FinalControlSize + 0.0001
    )
    
    mycovnames <- mycov
    names(mycovnames) <- 1:length(mycovnames)
    
    return(list(jdata = jdata, idata = idata, mycov = mycovnames))
}



get_waic <- function(jdf, jdata){
    jdf$iter2 <- 1:nrow(jdf)
    
    likNi <- jdf %>% 
        select(starts_with("logLik"), iter2) %>%
        gather(name, logLik, -iter2) %>% 
        separate(name, sep = "\\_", into = c("name", "iday")) %>% 
        select(-name) %>% 
        mutate(iday = as.numeric(as.character(iday)),
            lik = exp(logLik),
            region = jdata$jdata$FireCentre2N[iday],
            day = jdata$jdata$DayYearN[iday]) %>% 
        group_by(iter2, day) %>% 
        summarise(lik = prod(lik), .groups = "drop") #%>% pull(lik) %>% unique()
    
    likXij <- jdf %>% 
        select(starts_with("XLik"), iter2) %>% 
        gather(name, logLik, -iter2) %>% 
        separate(name, sep = "\\_", into = c("name", "jfire")) %>% 
        select(-name) %>% 
        mutate(jfire = as.numeric(as.character(jfire)),
            likX = exp(logLik),
            region = jdata$jdata$FireCentre2[jfire],
            day = jdata$jdata$DayYear[jfire]) %>% 
        group_by(iter2, day) %>% 
        summarise(likX = prod(likX), .groups = "drop") %>% 
        right_join(expand.grid(day = unique(likNi$day), iter2 = unique(likNi$iter2)),
            by = c("iter2", "day"))
    
    likXij$likX[is.na(likXij$likX)] <- 1
    
    full_join(likNi, likXij, by = c("iter2", "day")) %>% 
        mutate(ilik = lik*likX) %>% 
        group_by(day) %>% 
        summarise(lppd = log(mean(ilik)), p = var(log(ilik))) %>% 
        mutate(WAICi = -2*(lppd - p)) %>% 
        pull(WAICi) %>% sum()
}

get_waic_2 <- function(jdf, jdata){
    jdf$iter2 <- 1:nrow(jdf)
    
    likNi <- jdf %>% 
        select(starts_with("NLik"), iter2) %>%
        gather(name, Lik, -iter2) %>% 
        separate(name, sep = "\\_", into = c("name", "iday")) %>% 
        select(-name) %>% 
        mutate(iday = as.numeric(as.character(iday)),
            lik = Lik,
            region = jdata$jdata$FireCentre2N[iday],
            day = jdata$jdata$DayYearN[iday]) %>% head()
        group_by(iter2, day) %>% 
        summarise(lik = prod(lik)) #%>% pull(lik) %>% unique()
    
    likXij <- jdf %>% 
        select(starts_with("XLik"), iter2) %>% 
        gather(name, Lik, -iter2) %>% 
        separate(name, sep = "\\_", into = c("name", "jfire")) %>% 
        select(-name) %>% 
        mutate(jfire = as.numeric(as.character(jfire)),
            likX = Lik,
            region = jdata$jdata$FireCentre2[jfire],
            day = jdata$jdata$DayYear[jfire]) %>% 
        group_by(iter2, day) %>% 
        summarise(likX = prod(likX)) %>% 
        right_join(expand.grid(day = unique(likNi$day), iter2 = unique(likNi$iter2)),
            by = c("iter2", "day"))
    
    likXij$likX[is.na(likXij$likX)] <- 1
    
    full_join(likNi, likXij, by = c("iter2", "day")) %>% 
        mutate(ilik = lik*likX) %>% 
        group_by(day) %>% 
        summarise(lppd = log(mean(ilik, na.rm = TRUE)), 
            p = var(log(ilik), na.rm = TRUE)) %>% 
        mutate(WAICi = -2*(lppd - p)) %>% 
        pull(WAICi) %>% sum(na.rm = TRUE)
}


extract_NX <- function(jdf, jdata){
    OXi <- jdata$idata %>% unlist() - 0.0001
    ONi <- jdata$jdata$Ni %>% as.numeric()
    
    ENi <- jdf %>% 
        select(starts_with("ENi")) %>% 
        apply(2, median)
    EXi <- jdf %>% 
        select(starts_with("EXi")) %>% 
        apply(2, median)
    
    lppdN <- jdf %>% 
        select(starts_with("logLik")) %>% 
        apply(2, function(x) log(mean(exp(x))))
    pWAICN <- jdf %>% 
        select(starts_with("logLik")) %>% 
        apply(2, var)
    
    
    lppdX <- jdf %>% 
        select(starts_with("XLik")) %>% 
        apply(2, function(x) log(mean(exp(x))))
    pWAICX <- jdf %>% 
        select(starts_with("XLik")) %>% 
        apply(2, var)

	realWAIC <- get_waic(jdf = jdf, jdata = jdata)
    
    data.frame(data = deparse(substitute(jdata)),
        Xrmse = sqrt(mean((OXi - EXi)^2)),
        Nrmse = sqrt(mean((ONi - ENi)^2)),
        Xwaic = -2*(sum(lppdX - pWAICX)),
        Nwaic = -2*(sum(lppdN - pWAICN)),
		WAIC = realWAIC)
}



extract_gr <- function(codas){
    noms <- varnames(codas) %>% as.character()
    mainparams <- codas[, startsWith(noms, "beta") | startsWith(noms, "gammai") |
            startsWith(noms, "mux") | startsWith(noms, "sigma") |
            startsWith(noms, "lamb") | startsWith(noms, "pib") | 
            startsWith(noms, "alph") | startsWith(noms, "phi") | 
            startsWith(noms, "hyp")]
    
    # Because Variable selection makes the chol decomp singular
    corrections <- (apply(mainparams[[1]], 2, function(x) length(unique(x)) > 0.5*length(x))) |
        (apply(mainparams[[2]], 2, function(x) length(unique(x)) > 0.5*length(x)))
    corrections <- names(corrections)[corrections]
    
    mainparams <- codas[, corrections]
    
    gr1 <- gelman.diag(mainparams)
    ess1 <- effectiveSize(mainparams)
    
    data.frame(
        coda = deparse(substitute(codas)),
        miness = min(ess1[ess1 > 0]),
        meaness = mean(ess1),
        nparams = nrow(gr1$psrf),
        numOut1.1 = mean(gr1$psrf[,1] > 1.1),
        numOut1.15 = mean(gr1$psrf[,1] > 1.15),
        numCIOut1.1 = mean(gr1$psrf[,2] > 1.1),
        numCIOut1.15 = mean(gr1$psrf[,2] > 1.15),
        mpsrf = gr1$mpsrf
    )
}


# Wide and Rounded Models -----

# Province-wide splines
# Prov-wide re
# Independent re's
# Region-specific beta
# No hyperparameters
wide_jags_base <- {"
model{
    # For the ones trick
    # C chosen so that Lik/c <= 1
    # P(ones == 1) = p, where p is maximized
    # The trick: p is the actual likelihood, we don't care about bernoulli!
    C <- 1000000
    for(d in 1:ndays){
        b[d] ~ dnorm(0, taub)
        
        # Seasonal Affect
        lambda_spline[d] <- inprod(B[d,], lambeta[])
        # Presence spline
        pi_spline[d] <- inprod(B[d,], pibeta[])
        # Dispersion Spline
        alpha_spline[d] <- exp(inprod(B[d,], abeta[]))
    }
    
    for(i in 1:idays){
        # Ones trick
        ones[i] ~ dbern(p[i])
        p[i] <- Lik[i]/C
        Lik[i] <- exp(logLik[i])
        logLik[i] <- (1-z[i])*log(1-w[i]) + 
        #    z[i]*(log(w[i]) + logTPois[i])
            z[i]*(log(w[i]) + LogTruncNB[i])
        
        # For Hurdle
        logit(w[i]) <- pi_spline[DayYearN[i]]
        
        # Truncated Poisson Distribution
        #logTPois[i] <- log(w[i]) + Ni[i]*log(mu[i]) -
        #    mu[i] - log(1 - exp(-mu[i])) - loggam(Ni[i] + 1)
            
        # https://www.cambridge.org/core/services/aop-cambridge-core/content/view/C1FFB63EC7657CD424CE6726F32BC4FE/9781316459515c7_p184-214_CBO.pdf/glms_part_iii_zeroinflated_and_hurdle_models.pdf
        # https://ncss-wpengine.netdna-ssl.com/wp-content/themes/ncss/pdf/Procedures/NCSS/Zero-Inflated_Negative_Binomial_Regression.pdf
        LogTruncNB[i] <-  (1/alpha_spline[DayYearN[i]])*log(u[i])+
            Ni[i]*log(1 - u[i]) + 
            loggam(Ni[i] + 1/alpha_spline[DayYearN[i]]) -
            loggam(1/alpha_spline[DayYearN[i]]) - 
            loggam(Ni[i] + 1) -
            log(1 - (1 + alpha_spline[DayYearN[i]]*mu[i])^(-1/alpha_spline[DayYearN[i]]))
         u[i] <- 1/(1 + alpha_spline[DayYearN[i]]*mu[i])
         
        #https://data.princeton.edu/wws509/notes/countmoments
        ETrunc1[i] <- 1 + alpha_spline[DayYearN[i]]*mu[i]
        ETrunc2[i] <- -1/alpha_spline[DayYearN[i]]
        ETrunc3[i] <- ETrunc1[i]^ETrunc2[i]
        ETruncNB[i] <- mu[i]/(1 - ETrunc3[i])
        #ExTruncNB[i] <- mu[i]/((1-(1 + 
        #    alpha_spline[DayYearN[i]]*mu[i]))^(-1/alpha_spline[DayYearN[i]]))
        ENi[i] <- w[i]*ETruncNB[i]
        
        # mu is the mean of the Poisson
        # lambda[r] is region-specific intercept
        log(mu[i]) <- lambda_spline[DayYearN[i]] + lambda[FireCentre2N[i]] + 
            b[DayYearN[i]] #+ inprod(ncovars, betan)
    }
    
    
    for(j in 1:jfires){
        Xijr_star[j] ~ dlnorm(lp2[j], taux[FireCentre2[j]])
        lp1[j] <- inprod(xcovars[j,], betax[, FireCentre2[j]])
        lp2[j] <- mux[FireCentre2[j]] + 
                lp1[j] + 
                gammai[FireCentre2[j]] * b[DayYear[j]]
            
        Xijr1[j] ~ dinterval(Xijr_star[j], dint[j,])
        
        EXi[j] <- exp(lp2[j] + 1/(2*taux[FireCentre2[j]]))
        XLik[j] <- log(plnorm(dint[j, 2], lp2[j], taux[FireCentre2[j]]) -
            plnorm(dint[j, 1], lp2[j], taux[FireCentre2[j]]))
    }
    
    # Priors
    taub <- pow(sigma_b, -2)
    sigma_b ~ dgamma(8, 2)
    
    # Splines
    lambeta ~ dmnorm(zero_vec, lambeta_posdef)
    pibeta ~ dmnorm(zero_vec, lambeta_posdef)
    abeta ~ dmnorm(zero_vec, lambeta_posdef)
    

    
    for(p in 1:nxcovar){
        for(r in 1:6){
            betax[p, r] ~ dnorm(0, 1/4)
        }
    }
    
    #for(p in 1:nNcovar){
    #    betan[p] ~ dnorm(0, 1/4)
    #}
    
    for(r in 1:6){
        gammai[r] ~ dnorm(0, 1/5)
        mux[r] ~ dnorm(0, 1/5)
        lambda[r] ~ dnorm(0, 1/5)
        taux[r] <- pow(sigma_x[r], -2)
        sigma_x[r] ~ dgamma(8, 2)
        #roundprob[r] ~ dbeta(2,2)
        
    }
}
"}


# Province-wide splines
# Prov-wide re
# Independent re's
# Region-specific beta
# Yes hyperparameters
wide_jags_hyper <- {"
model{
    # For the ones trick
    # C chosen so that Lik/c <= 1
    # P(ones == 1) = p, where p is maximized
    # The trick: p is the actual likelihood, we don't care about bernoulli!
    C <- 1000000
    for(d in 1:ndays){
        b[d] ~ dnorm(0, taub)
        
        # Seasonal Affect
        lambda_spline[d] <- inprod(B[d,], lambeta[])
        # Presence spline
        pi_spline[d] <- inprod(B[d,], pibeta[])
        # Dispersion Spline
        alpha_spline[d] <- exp(inprod(B[d,], abeta[]))
    }
    
    for(i in 1:idays){
        # Ones trick
        ones[i] ~ dbern(p[i])
        p[i] <- Lik[i]/C
        Lik[i] <- exp(logLik[i])
        logLik[i] <- (1-z[i])*log(1-w[i]) + 
        #    z[i]*(log(w[i]) + logTPois[i])
            z[i]*(log(w[i]) + LogTruncNB[i])
        
        # For Hurdle
        logit(w[i]) <- pi_spline[DayYearN[i]]
        
        # Truncated Poisson Distribution
        #logTPois[i] <- log(w[i]) + Ni[i]*log(mu[i]) -
        #    mu[i] - log(1 - exp(-mu[i])) - loggam(Ni[i] + 1)
            
        # https://www.cambridge.org/core/services/aop-cambridge-core/content/view/C1FFB63EC7657CD424CE6726F32BC4FE/9781316459515c7_p184-214_CBO.pdf/glms_part_iii_zeroinflated_and_hurdle_models.pdf
        # https://ncss-wpengine.netdna-ssl.com/wp-content/themes/ncss/pdf/Procedures/NCSS/Zero-Inflated_Negative_Binomial_Regression.pdf
        LogTruncNB[i] <-  (1/alpha_spline[DayYearN[i]])*log(u[i])+
            Ni[i]*log(1 - u[i]) + 
            loggam(Ni[i] + 1/alpha_spline[DayYearN[i]]) -
            loggam(1/alpha_spline[DayYearN[i]]) - 
            loggam(Ni[i] + 1) -
            log(1 - (1 + alpha_spline[DayYearN[i]]*mu[i])^(-1/alpha_spline[DayYearN[i]]))
         u[i] <- 1/(1 + alpha_spline[DayYearN[i]]*mu[i])
         
        #https://data.princeton.edu/wws509/notes/countmoments
        ETrunc1[i] <- 1 + alpha_spline[DayYearN[i]]*mu[i]
        ETrunc2[i] <- -1/alpha_spline[DayYearN[i]]
        ETrunc3[i] <- ETrunc1[i]^ETrunc2[i]
        ETruncNB[i] <- mu[i]/(1 - ETrunc3[i])
        #ExTruncNB[i] <- mu[i]/((1-(1 + 
        #    alpha_spline[DayYearN[i]]*mu[i]))^(-1/alpha_spline[DayYearN[i]]))
        ENi[i] <- w[i]*ETruncNB[i]
        
        # mu is the mean of the Poisson
        # lambda[r] is region-specific intercept
        log(mu[i]) <- lambda_spline[DayYearN[i]] + lambda[FireCentre2N[i]] + 
            b[DayYearN[i]] #+ inprod(ncovars, betan)
    }
    
    
    for(j in 1:jfires){
        Xijr_star[j] ~ dlnorm(lp2[j], taux[FireCentre2[j]])
        lp1[j] <- inprod(xcovars[j,], betax[, FireCentre2[j]])
        lp2[j] <- mux[FireCentre2[j]] + 
                lp1[j] + 
                gammai[FireCentre2[j]] * b[DayYear[j]]
            
        Xijr1[j] ~ dinterval(Xijr_star[j], dint[j,])
        
        EXi[j] <- exp(lp2[j] + 1/(2*taux[FireCentre2[j]]))
        XLik[j] <- log(plnorm(dint[j, 2], lp2[j], taux[FireCentre2[j]]) -
            plnorm(dint[j, 1], lp2[j], taux[FireCentre2[j]]))
    }
    
    # Priors
    taub <- pow(sigma_b, -2)
    sigma_b ~ dgamma(8, 2)
    
    # Splines
    lambeta ~ dmnorm(zero_vec, lambeta_posdef)
    pibeta ~ dmnorm(zero_vec, lambeta_posdef)
    abeta ~ dmnorm(zero_vec, lambeta_posdef)
    

    
    for(p in 1:nxcovar){
        for(r in 1:6){
            betax[p, r] ~ dnorm(0, 1/4)
        }
    }
    
    #for(p in 1:nNcovar){
    #    betan[p] ~ dnorm(0, 1/4)
    #}
    
    for(r in 1:6){
        gammai[r] ~ dnorm(hypgamma, 1/5)
        mux[r] ~ dnorm(hypmux, 1/5)
        lambda[r] ~ dnorm(hyplambda, 1/5)
        taux[r] <- pow(sigma_x[r], -2)
        sigma_x[r] ~ dgamma(8, 2)
        #roundprob[r] ~ dbeta(2,2)
        
    }
	hypgamma ~ dnorm(0, 1/5)
	hypmux ~ dnorm(0, 1/5)
	hyplambda ~ dnorm(0, 1/5)
}
"}


# Province-wide splines
# Prov-wide re
# AR1 re's
# Region-specific beta
# No hyperparameters
wide_jags_AR1 <- {"
model{
    # For the ones trick
    # C chosen so that Lik/c <= 1
    # P(ones == 1) = p, where p is maximized
    # The trick: p is the actual likelihood, we don't care about bernoulli!
    C <- 1000000
    for(d in 1:ndays){
        # Seasonal Affect
        lambda_spline[d] <- inprod(B[d,], lambeta[])
        # Presence spline
        pi_spline[d] <- inprod(B[d,], pibeta[])
        # Dispersion Spline
        alpha_spline[d] <- exp(inprod(B[d,], abeta[]))
    }
    
    for(i in 1:idays){
        # Ones trick
        ones[i] ~ dbern(p[i])
        p[i] <- Lik[i]/C
        Lik[i] <- exp(logLik[i])
        logLik[i] <- (1-z[i])*log(1-w[i]) + 
        #    z[i]*(log(w[i]) + logTPois[i])
            z[i]*(log(w[i]) + LogTruncNB[i])
        
        # For Hurdle
        logit(w[i]) <- pi_spline[DayYearN[i]]
        
        # Truncated Poisson Distribution
        #logTPois[i] <- log(w[i]) + Ni[i]*log(mu[i]) -
        #    mu[i] - log(1 - exp(-mu[i])) - loggam(Ni[i] + 1)
            
        # https://www.cambridge.org/core/services/aop-cambridge-core/content/view/C1FFB63EC7657CD424CE6726F32BC4FE/9781316459515c7_p184-214_CBO.pdf/glms_part_iii_zeroinflated_and_hurdle_models.pdf
        # https://ncss-wpengine.netdna-ssl.com/wp-content/themes/ncss/pdf/Procedures/NCSS/Zero-Inflated_Negative_Binomial_Regression.pdf
        LogTruncNB[i] <-  (1/alpha_spline[DayYearN[i]])*log(u[i])+
            Ni[i]*log(1 - u[i]) + 
            loggam(Ni[i] + 1/alpha_spline[DayYearN[i]]) -
            loggam(1/alpha_spline[DayYearN[i]]) - 
            loggam(Ni[i] + 1) -
            log(1 - (1 + alpha_spline[DayYearN[i]]*mu[i])^(-1/alpha_spline[DayYearN[i]]))
         u[i] <- 1/(1 + alpha_spline[DayYearN[i]]*mu[i])
         
        #https://data.princeton.edu/wws509/notes/countmoments
        ETrunc1[i] <- 1 + alpha_spline[DayYearN[i]]*mu[i]
        ETrunc2[i] <- -1/alpha_spline[DayYearN[i]]
        ETrunc3[i] <- ETrunc1[i]^ETrunc2[i]
        ETruncNB[i] <- mu[i]/(1 - ETrunc3[i])
        #ExTruncNB[i] <- mu[i]/((1-(1 + 
        #    alpha_spline[DayYearN[i]]*mu[i]))^(-1/alpha_spline[DayYearN[i]]))
        ENi[i] <- w[i]*ETruncNB[i]
        
        # mu is the mean of the Poisson
        # lambda[r] is region-specific intercept
        log(mu[i]) <- lambda_spline[DayYearN[i]] + lambda[FireCentre2N[i]] + 
            b[DayYearN[i]] #+ inprod(ncovars, betan)
    }
    
    
    for(j in 1:jfires){
        Xijr_star[j] ~ dlnorm(lp2[j], taux[FireCentre2[j]])
        lp1[j] <- inprod(xcovars[j,], betax[, FireCentre2[j]])
        lp2[j] <- mux[FireCentre2[j]] + 
                lp1[j] + 
                gammai[FireCentre2[j]] * b[DayYear[j]]
            
        Xijr1[j] ~ dinterval(Xijr_star[j], dint[j,])
        
        EXi[j] <- exp(lp2[j] + 1/(2*taux[FireCentre2[j]]))
        XLik[j] <- log(plnorm(dint[j, 2], lp2[j], taux[FireCentre2[j]]) -
            plnorm(dint[j, 1], lp2[j], taux[FireCentre2[j]]))
    }
    
    # Priors
    taub <- pow(sigma_b, -2)
    sigma_b ~ dgamma(8, 2)
    
    # Splines
    lambeta ~ dmnorm(zero_vec, lambeta_posdef)
    pibeta ~ dmnorm(zero_vec, lambeta_posdef)
    abeta ~ dmnorm(zero_vec, lambeta_posdef)
    
    # AR(1) structure for random effects
    phi[1] ~ dnorm(0, 1/5)
    b[1] ~ dnorm(0, taub)
    for(i in 2:ndays){
        b[i] ~ dnorm(phi[1]*b[i-1], taub)
    }

    
    for(p in 1:nxcovar){
        for(r in 1:6){
            betax[p, r] ~ dnorm(0, 1/4)
        }
    }
    
    #for(p in 1:nNcovar){
    #    betan[p] ~ dnorm(0, 1/4)
    #}
    
    for(r in 1:6){
        gammai[r] ~ dnorm(0, 1/5)
        mux[r] ~ dnorm(0, 1/5)
        lambda[r] ~ dnorm(0, 1/5)
        taux[r] <- pow(sigma_x[r], -2)
        sigma_x[r] ~ dgamma(8, 2)
        #roundprob[r] ~ dbeta(2,2)
        
    }
}
"}


# Province-wide splines
# Prov-wide re
# AR2 re's
# Region-specific beta
# No hyperparameters
wide_jags_AR2 <- {"
model{
    # For the ones trick
    # C chosen so that Lik/c <= 1
    # P(ones == 1) = p, where p is maximized
    # The trick: p is the actual likelihood, we don't care about bernoulli!
    C <- 1000000
    for(d in 1:ndays){
        # Seasonal Affect
        lambda_spline[d] <- inprod(B[d,], lambeta[])
        # Presence spline
        pi_spline[d] <- inprod(B[d,], pibeta[])
        # Dispersion Spline
        alpha_spline[d] <- exp(inprod(B[d,], abeta[]))
    }
    
    for(i in 1:idays){
        # Ones trick
        ones[i] ~ dbern(p[i])
        p[i] <- Lik[i]/C
        Lik[i] <- exp(logLik[i])
        logLik[i] <- (1-z[i])*log(1-w[i]) + 
        #    z[i]*(log(w[i]) + logTPois[i])
            z[i]*(log(w[i]) + LogTruncNB[i])
        
        # For Hurdle
        logit(w[i]) <- pi_spline[DayYearN[i]]
        
        # Truncated Poisson Distribution
        #logTPois[i] <- log(w[i]) + Ni[i]*log(mu[i]) -
        #    mu[i] - log(1 - exp(-mu[i])) - loggam(Ni[i] + 1)
            
        # https://www.cambridge.org/core/services/aop-cambridge-core/content/view/C1FFB63EC7657CD424CE6726F32BC4FE/9781316459515c7_p184-214_CBO.pdf/glms_part_iii_zeroinflated_and_hurdle_models.pdf
        # https://ncss-wpengine.netdna-ssl.com/wp-content/themes/ncss/pdf/Procedures/NCSS/Zero-Inflated_Negative_Binomial_Regression.pdf
        LogTruncNB[i] <-  (1/alpha_spline[DayYearN[i]])*log(u[i])+
            Ni[i]*log(1 - u[i]) + 
            loggam(Ni[i] + 1/alpha_spline[DayYearN[i]]) -
            loggam(1/alpha_spline[DayYearN[i]]) - 
            loggam(Ni[i] + 1) -
            log(1 - (1 + alpha_spline[DayYearN[i]]*mu[i])^(-1/alpha_spline[DayYearN[i]]))
         u[i] <- 1/(1 + alpha_spline[DayYearN[i]]*mu[i])
         
        #https://data.princeton.edu/wws509/notes/countmoments
        ETrunc1[i] <- 1 + alpha_spline[DayYearN[i]]*mu[i]
        ETrunc2[i] <- -1/alpha_spline[DayYearN[i]]
        ETrunc3[i] <- ETrunc1[i]^ETrunc2[i]
        ETruncNB[i] <- mu[i]/(1 - ETrunc3[i])
        #ExTruncNB[i] <- mu[i]/((1-(1 + 
        #    alpha_spline[DayYearN[i]]*mu[i]))^(-1/alpha_spline[DayYearN[i]]))
        ENi[i] <- w[i]*ETruncNB[i]
        
        # mu is the mean of the Poisson
        # lambda[r] is region-specific intercept
        log(mu[i]) <- lambda_spline[DayYearN[i]] + lambda[FireCentre2N[i]] + 
            b[DayYearN[i]] #+ inprod(ncovars, betan)
    }
    
    
    for(j in 1:jfires){
        Xijr_star[j] ~ dlnorm(lp2[j], taux[FireCentre2[j]])
        lp1[j] <- inprod(xcovars[j,], betax[, FireCentre2[j]])
        lp2[j] <- mux[FireCentre2[j]] + 
                lp1[j] + 
                gammai[FireCentre2[j]] * b[DayYear[j]]
            
        Xijr1[j] ~ dinterval(Xijr_star[j], dint[j,])
        
        EXi[j] <- exp(lp2[j] + 1/(2*taux[FireCentre2[j]]))
        XLik[j] <- log(plnorm(dint[j, 2], lp2[j], taux[FireCentre2[j]]) -
            plnorm(dint[j, 1], lp2[j], taux[FireCentre2[j]]))
    }
    
    # Priors
    taub <- pow(sigma_b, -2)
    sigma_b ~ dgamma(8, 2)
    
    # Splines
    lambeta ~ dmnorm(zero_vec, lambeta_posdef)
    pibeta ~ dmnorm(zero_vec, lambeta_posdef)
    abeta ~ dmnorm(zero_vec, lambeta_posdef)
    
    # AR(2) structure for random effects
    phi[1] ~ dnorm(0, 1/5)
    phi[2] ~ dnorm(0, 1/5)
    b[1] ~ dnorm(0, taub)
    b[2] ~ dnorm(phi[1]*b[1], taub)
    for(i in 3:idays){
        b[i] ~ dnorm(phi[1]*b[i-1] + phi[2]*b[i-2], taub)
    }

    
    for(p in 1:nxcovar){
        for(r in 1:6){
            betax[p, r] ~ dnorm(0, 1/4)
        }
    }
    
    #for(p in 1:nNcovar){
    #    betan[p] ~ dnorm(0, 1/4)
    #}
    
    for(r in 1:6){
        gammai[r] ~ dnorm(0, 1/5)
        mux[r] ~ dnorm(0, 1/5)
        lambda[r] ~ dnorm(0, 1/5)
        taux[r] <- pow(sigma_x[r], -2)
        sigma_x[r] ~ dgamma(8, 2)
        #roundprob[r] ~ dbeta(2,2)
        
    }
}
"}


# Province-wide splines
	# alpha is not a spline
# Prov-wide re
# Independent re's
# Region-specific beta
# No hyperparameters
wide_jags_daily_alpha <- {"
model{
    # For the ones trick
    # C chosen so that Lik/c <= 1
    # P(ones == 1) = p, where p is maximized
    # The trick: p is the actual likelihood, we don't care about bernoulli!
    C <- 1000000
    for(d in 1:ndays){
        b[d] ~ dnorm(0, taub)
        
        # Seasonal Affect
        lambda_spline[d] <- inprod(B[d,], lambeta[])
        # Presence spline
        pi_spline[d] <- inprod(B[d,], pibeta[])
        # Dispersion Spline
        #alpha_spline[d] <- exp(inprod(B[d,], abeta[]))
		alpha[d] ~ dunif(0, 100)
    }
    
    for(i in 1:idays){
        # Ones trick
        ones[i] ~ dbern(p[i])
        p[i] <- Lik[i]/C
        Lik[i] <- exp(logLik[i])
        logLik[i] <- (1-z[i])*log(1-w[i]) + 
        #    z[i]*(log(w[i]) + logTPois[i])
            z[i]*(log(w[i]) + LogTruncNB[i])
        
        # For Hurdle
        logit(w[i]) <- pi_spline[DayYearN[i]]
        
        # Truncated Poisson Distribution
        #logTPois[i] <- log(w[i]) + Ni[i]*log(mu[i]) -
        #    mu[i] - log(1 - exp(-mu[i])) - loggam(Ni[i] + 1)
            
        # https://www.cambridge.org/core/services/aop-cambridge-core/content/view/C1FFB63EC7657CD424CE6726F32BC4FE/9781316459515c7_p184-214_CBO.pdf/glms_part_iii_zeroinflated_and_hurdle_models.pdf
        # https://ncss-wpengine.netdna-ssl.com/wp-content/themes/ncss/pdf/Procedures/NCSS/Zero-Inflated_Negative_Binomial_Regression.pdf
        LogTruncNB[i] <-  (1/alpha[DayYearN[i]])*log(u[i])+
            Ni[i]*log(1 - u[i]) + 
            loggam(Ni[i] + 1/alpha[DayYearN[i]]) -
            loggam(1/alpha[DayYearN[i]]) - 
            loggam(Ni[i] + 1) -
            log(1 - (1 + alpha[DayYearN[i]]*mu[i])^(-1/alpha[DayYearN[i]]))
         u[i] <- 1/(1 + alpha[DayYearN[i]]*mu[i])
         
        #https://data.princeton.edu/wws509/notes/countmoments
        ETrunc1[i] <- 1 + alpha[DayYearN[i]]*mu[i]
        ETrunc2[i] <- -1/alpha[DayYearN[i]]
        ETrunc3[i] <- ETrunc1[i]^ETrunc2[i]
        ETruncNB[i] <- mu[i]/(1 - ETrunc3[i])
        #ExTruncNB[i] <- mu[i]/((1-(1 + 
        #    alpha[DayYearN[i]]*mu[i]))^(-1/alpha[DayYearN[i]]))
        ENi[i] <- w[i]*ETruncNB[i]
        
        # mu is the mean of the Poisson
        # lambda[r] is region-specific intercept
        log(mu[i]) <- lambda_spline[DayYearN[i]] + lambda[FireCentre2N[i]] + 
            b[DayYearN[i]] #+ inprod(ncovars, betan)
    }
    
    
    for(j in 1:jfires){
        Xijr_star[j] ~ dlnorm(lp2[j], taux[FireCentre2[j]])
        lp1[j] <- inprod(xcovars[j,], betax[, FireCentre2[j]])
        lp2[j] <- mux[FireCentre2[j]] + 
                lp1[j] + 
                gammai[FireCentre2[j]] * b[DayYear[j]]
            
        Xijr1[j] ~ dinterval(Xijr_star[j], dint[j,])
        
        EXi[j] <- exp(lp2[j] + 1/(2*taux[FireCentre2[j]]))
        XLik[j] <- log(plnorm(dint[j, 2], lp2[j], taux[FireCentre2[j]]) -
            plnorm(dint[j, 1], lp2[j], taux[FireCentre2[j]]))
    }
    
    # Priors
    taub <- pow(sigma_b, -2)
    sigma_b ~ dgamma(8, 2)
    
    # Splines
    lambeta ~ dmnorm(zero_vec, lambeta_posdef)
    pibeta ~ dmnorm(zero_vec, lambeta_posdef)
    #abeta ~ dmnorm(zero_vec, lambeta_posdef)
    

    
    for(p in 1:nxcovar){
        for(r in 1:6){
            betax[p, r] ~ dnorm(0, 1/4)
        }
    }
    
    #for(p in 1:nNcovar){
    #    betan[p] ~ dnorm(0, 1/4)
    #}
    
    for(r in 1:6){
        gammai[r] ~ dnorm(0, 1/5)
        mux[r] ~ dnorm(0, 1/5)
        lambda[r] ~ dnorm(0, 1/5)
        taux[r] <- pow(sigma_x[r], -2)
        sigma_x[r] ~ dgamma(8, 2)
        #roundprob[r] ~ dbeta(2,2)
        
    }
}
"}

# Province-wide splines
# Regional re
# Independent re's
# Region-specific beta
# No hyperparameters
wide_jags_regionre <- {"
model{
    # For the ones trick
    # C chosen so that Lik/c <= 1
    # P(ones == 1) = p, where p is maximized
    # The trick: p is the actual likelihood, we don't care about bernoulli!
    C <- 1000000
    for(d in 1:ndays){
		for(r in 1:6){
			b[r, d] ~ dnorm(0, taub)
		}
        
        # Seasonal Affect
        lambda_spline[d] <- inprod(B[d,], lambeta[])
        # Presence spline
        pi_spline[d] <- inprod(B[d,], pibeta[])
        # Dispersion Spline
        alpha_spline[d] <- exp(inprod(B[d,], abeta[]))
    }
    
    for(i in 1:idays){
        # Ones trick
        ones[i] ~ dbern(p[i])
        p[i] <- Lik[i]/C
        Lik[i] <- exp(logLik[i])
        logLik[i] <- (1-z[i])*log(1-w[i]) + 
        #    z[i]*(log(w[i]) + logTPois[i])
            z[i]*(log(w[i]) + LogTruncNB[i])
        
        # For Hurdle
        logit(w[i]) <- pi_spline[DayYearN[i]]
        
        # Truncated Poisson Distribution
        #logTPois[i] <- log(w[i]) + Ni[i]*log(mu[i]) -
        #    mu[i] - log(1 - exp(-mu[i])) - loggam(Ni[i] + 1)
            
        # https://www.cambridge.org/core/services/aop-cambridge-core/content/view/C1FFB63EC7657CD424CE6726F32BC4FE/9781316459515c7_p184-214_CBO.pdf/glms_part_iii_zeroinflated_and_hurdle_models.pdf
        # https://ncss-wpengine.netdna-ssl.com/wp-content/themes/ncss/pdf/Procedures/NCSS/Zero-Inflated_Negative_Binomial_Regression.pdf
        LogTruncNB[i] <-  (1/alpha_spline[DayYearN[i]])*log(u[i])+
            Ni[i]*log(1 - u[i]) + 
            loggam(Ni[i] + 1/alpha_spline[DayYearN[i]]) -
            loggam(1/alpha_spline[DayYearN[i]]) - 
            loggam(Ni[i] + 1) -
            log(1 - (1 + alpha_spline[DayYearN[i]]*mu[i])^(-1/alpha_spline[DayYearN[i]]))
         u[i] <- 1/(1 + alpha_spline[DayYearN[i]]*mu[i])
         
        #https://data.princeton.edu/wws509/notes/countmoments
        ETrunc1[i] <- 1 + alpha_spline[DayYearN[i]]*mu[i]
        ETrunc2[i] <- -1/alpha_spline[DayYearN[i]]
        ETrunc3[i] <- ETrunc1[i]^ETrunc2[i]
        ETruncNB[i] <- mu[i]/(1 - ETrunc3[i])
        #ExTruncNB[i] <- mu[i]/((1-(1 + 
        #    alpha_spline[DayYearN[i]]*mu[i]))^(-1/alpha_spline[DayYearN[i]]))
        ENi[i] <- w[i]*ETruncNB[i]
        
        # mu is the mean of the Poisson
        # lambda[r] is region-specific intercept
        log(mu[i]) <- lambda_spline[DayYearN[i]] + lambda[FireCentre2N[i]] + 
            b[FireCentre2N[i], DayYearN[i]] #+ inprod(ncovars, betan)
    }
    
    
    for(j in 1:jfires){
        Xijr_star[j] ~ dlnorm(lp2[j], taux[FireCentre2[j]])
        lp1[j] <- inprod(xcovars[j,], betax[, FireCentre2[j]])
        lp2[j] <- mux[FireCentre2[j]] + 
                lp1[j] + 
                gammai[FireCentre2[j]] * b[FireCentre2[j], DayYear[j]]
            
        Xijr1[j] ~ dinterval(Xijr_star[j], dint[j,])
        
        EXi[j] <- exp(lp2[j] + 1/(2*taux[FireCentre2[j]]))
        XLik[j] <- log(plnorm(dint[j, 2], lp2[j], taux[FireCentre2[j]]) -
            plnorm(dint[j, 1], lp2[j], taux[FireCentre2[j]]))
    }
    
    # Priors
    taub <- pow(sigma_b, -2)
    sigma_b ~ dgamma(8, 2)
    
    # Splines
    lambeta ~ dmnorm(zero_vec, lambeta_posdef)
    pibeta ~ dmnorm(zero_vec, lambeta_posdef)
    abeta ~ dmnorm(zero_vec, lambeta_posdef)
    

    
    for(p in 1:nxcovar){
        for(r in 1:6){
            betax[p, r] ~ dnorm(0, 1/4)
        }
    }
    
    #for(p in 1:nNcovar){
    #    betan[p] ~ dnorm(0, 1/4)
    #}
    
    for(r in 1:6){
        gammai[r] ~ dnorm(0, 1/5)
        mux[r] ~ dnorm(0, 1/5)
        lambda[r] ~ dnorm(0, 1/5)
        taux[r] <- pow(sigma_x[r], -2)
        sigma_x[r] ~ dgamma(8, 2)
        #roundprob[r] ~ dbeta(2,2)
        
    }
}
"}


# Province-wide splines
# Prov-wide re
# AR1 re's
# Prov-wide beta
# Yes hyperparameters
wide_jags_AR1pbeta <- {"
model{
    # For the ones trick
    # C chosen so that Lik/c <= 1
    # P(ones == 1) = p, where p is maximized
    # The trick: p is the actual likelihood, we don't care about bernoulli!
    C <- 1000000
    for(d in 1:ndays){
        # Seasonal Affect
        lambda_spline[d] <- inprod(B[d,], lambeta[])
        # Presence spline
        pi_spline[d] <- inprod(B[d,], pibeta[])
        # Dispersion Spline
        alpha_spline[d] <- exp(inprod(B[d,], abeta[]))
    }
    
    for(i in 1:idays){
        # Ones trick
        ones[i] ~ dbern(p[i])
        p[i] <- Lik[i]/C
        Lik[i] <- exp(logLik[i])
        logLik[i] <- (1-z[i])*log(1-w[i]) + 
        #    z[i]*(log(w[i]) + logTPois[i])
            z[i]*(log(w[i]) + LogTruncNB[i])
        
        # For Hurdle
        logit(w[i]) <- pi_spline[DayYearN[i]]
        
        # Truncated Poisson Distribution
        #logTPois[i] <- log(w[i]) + Ni[i]*log(mu[i]) -
        #    mu[i] - log(1 - exp(-mu[i])) - loggam(Ni[i] + 1)
            
        # https://www.cambridge.org/core/services/aop-cambridge-core/content/view/C1FFB63EC7657CD424CE6726F32BC4FE/9781316459515c7_p184-214_CBO.pdf/glms_part_iii_zeroinflated_and_hurdle_models.pdf
        # https://ncss-wpengine.netdna-ssl.com/wp-content/themes/ncss/pdf/Procedures/NCSS/Zero-Inflated_Negative_Binomial_Regression.pdf
        LogTruncNB[i] <-  (1/alpha_spline[DayYearN[i]])*log(u[i])+
            Ni[i]*log(1 - u[i]) + 
            loggam(Ni[i] + 1/alpha_spline[DayYearN[i]]) -
            loggam(1/alpha_spline[DayYearN[i]]) - 
            loggam(Ni[i] + 1) -
            log(1 - (1 + alpha_spline[DayYearN[i]]*mu[i])^(-1/alpha_spline[DayYearN[i]]))
         u[i] <- 1/(1 + alpha_spline[DayYearN[i]]*mu[i])
         
        #https://data.princeton.edu/wws509/notes/countmoments
        ETrunc1[i] <- 1 + alpha_spline[DayYearN[i]]*mu[i]
        ETrunc2[i] <- -1/alpha_spline[DayYearN[i]]
        ETrunc3[i] <- ETrunc1[i]^ETrunc2[i]
        ETruncNB[i] <- mu[i]/(1 - ETrunc3[i])
        #ExTruncNB[i] <- mu[i]/((1-(1 + 
        #    alpha_spline[DayYearN[i]]*mu[i]))^(-1/alpha_spline[DayYearN[i]]))
        ENi[i] <- w[i]*ETruncNB[i]
        
        # mu is the mean of the Poisson
        # lambda[r] is region-specific intercept
        log(mu[i]) <- lambda_spline[DayYearN[i]] + lambda[FireCentre2N[i]] + 
            b[DayYearN[i]] #+ inprod(ncovars, betan)
    }
    
    
    for(j in 1:jfires){
        Xijr_star[j] ~ dlnorm(lp2[j], taux[FireCentre2[j]])
        lp1[j] <- inprod(xcovars[j,], betax)
        lp2[j] <- mux[FireCentre2[j]] + 
                lp1[j] + 
                gammai[FireCentre2[j]] * b[DayYear[j]]
            
        Xijr1[j] ~ dinterval(Xijr_star[j], dint[j,])
        
        EXi[j] <- exp(lp2[j] + 1/(2*taux[FireCentre2[j]]))
        XLik[j] <- log(plnorm(dint[j, 2], lp2[j], taux[FireCentre2[j]]) -
            plnorm(dint[j, 1], lp2[j], taux[FireCentre2[j]]))
    }
    
    # Priors
    taub <- pow(sigma_b, -2)
    sigma_b ~ dgamma(8, 2)
    
    # Splines
    lambeta ~ dmnorm(zero_vec, lambeta_posdef)
    pibeta ~ dmnorm(zero_vec, lambeta_posdef)
    abeta ~ dmnorm(zero_vec, lambeta_posdef)
    
    # AR(1) structure for random effects
    phi[1] ~ dnorm(0, 1/5)
    b[1] ~ dnorm(0, taub)
    for(i in 2:ndays){
        b[i] ~ dnorm(phi[1]*b[i-1], taub)
    }

    
    for(p in 1:nxcovar){
		betax[p] ~ dnorm(0, 1/4)
	}
    
    #for(p in 1:nNcovar){
    #    betan[p] ~ dnorm(0, 1/4)
    #}
    
    for(r in 1:6){
        gammai[r] ~ dnorm(hypgamma, 1/5)
        mux[r] ~ dnorm(hypmux, 1/5)
        lambda[r] ~ dnorm(hyplambda, 1/5)
        taux[r] <- pow(sigma_x[r], -2)
        sigma_x[r] ~ dgamma(8, 2)
        roundprob[r] ~ dbeta(2,2)
        
    }
	hypgamma ~ dnorm(0, 1/5)
	hypmux ~ dnorm(0, 1/5)
	hyplambda ~ dnorm(0, 1/5)
}
"}

# gamma is set to exactly 0
wide_jags_AR1pbeta_sep <- {"
model{
    # For the ones trick
    # C chosen so that Lik/c <= 1
    # P(ones == 1) = p, where p is maximized
    # The trick: p is the actual likelihood, we don't care about bernoulli!
    C <- 1000000
    for(d in 1:ndays){
        # Seasonal Affect
        lambda_spline[d] <- inprod(B[d,], lambeta[])
        # Presence spline
        pi_spline[d] <- inprod(B[d,], pibeta[])
        # Dispersion Spline
        alpha_spline[d] <- exp(inprod(B[d,], abeta[]))
    }
    
    for(i in 1:idays){
        # Ones trick
        ones[i] ~ dbern(p[i])
        p[i] <- Lik[i]/C
        Lik[i] <- exp(logLik[i])
        logLik[i] <- (1-z[i])*log(1-w[i]) + 
        #    z[i]*(log(w[i]) + logTPois[i])
            z[i]*(log(w[i]) + LogTruncNB[i])
        
        # For Hurdle
        logit(w[i]) <- pi_spline[DayYearN[i]]
        
        # Truncated Poisson Distribution
        #logTPois[i] <- log(w[i]) + Ni[i]*log(mu[i]) -
        #    mu[i] - log(1 - exp(-mu[i])) - loggam(Ni[i] + 1)
            
        # https://www.cambridge.org/core/services/aop-cambridge-core/content/view/C1FFB63EC7657CD424CE6726F32BC4FE/9781316459515c7_p184-214_CBO.pdf/glms_part_iii_zeroinflated_and_hurdle_models.pdf
        # https://ncss-wpengine.netdna-ssl.com/wp-content/themes/ncss/pdf/Procedures/NCSS/Zero-Inflated_Negative_Binomial_Regression.pdf
        LogTruncNB[i] <-  (1/alpha_spline[DayYearN[i]])*log(u[i])+
            Ni[i]*log(1 - u[i]) + 
            loggam(Ni[i] + 1/alpha_spline[DayYearN[i]]) -
            loggam(1/alpha_spline[DayYearN[i]]) - 
            loggam(Ni[i] + 1) -
            log(1 - (1 + alpha_spline[DayYearN[i]]*mu[i])^(-1/alpha_spline[DayYearN[i]]))
         u[i] <- 1/(1 + alpha_spline[DayYearN[i]]*mu[i])
         
        #https://data.princeton.edu/wws509/notes/countmoments
        ETrunc1[i] <- 1 + alpha_spline[DayYearN[i]]*mu[i]
        ETrunc2[i] <- -1/alpha_spline[DayYearN[i]]
        ETrunc3[i] <- ETrunc1[i]^ETrunc2[i]
        ETruncNB[i] <- mu[i]/(1 - ETrunc3[i])
        #ExTruncNB[i] <- mu[i]/((1-(1 + 
        #    alpha_spline[DayYearN[i]]*mu[i]))^(-1/alpha_spline[DayYearN[i]]))
        ENi[i] <- w[i]*ETruncNB[i]
        
        # mu is the mean of the Poisson
        # lambda[r] is region-specific intercept
        log(mu[i]) <- lambda_spline[DayYearN[i]] + lambda[FireCentre2N[i]] + 
            b[DayYearN[i]] #+ inprod(ncovars, betan)
    }
    
    
    for(j in 1:jfires){
        Xijr_star[j] ~ dlnorm(lp2[j], taux[FireCentre2[j]])
        lp1[j] <- inprod(xcovars[j,], betax)
        lp2[j] <- mux[FireCentre2[j]] + 
                lp1[j] + 
                gammai[FireCentre2[j]] * b[DayYear[j]]
            
        Xijr1[j] ~ dinterval(Xijr_star[j], dint[j,])
        
        EXi[j] <- exp(lp2[j] + 1/(2*taux[FireCentre2[j]]))
        XLik[j] <- log(plnorm(dint[j, 2], lp2[j], taux[FireCentre2[j]]) -
            plnorm(dint[j, 1], lp2[j], taux[FireCentre2[j]]))
    }
    
    # Priors
    taub <- pow(sigma_b, -2)
    sigma_b ~ dgamma(8, 2)
    
    # Splines
    lambeta ~ dmnorm(zero_vec, lambeta_posdef)
    pibeta ~ dmnorm(zero_vec, lambeta_posdef)
    abeta ~ dmnorm(zero_vec, lambeta_posdef)
    
    # AR(1) structure for random effects
    phi[1] ~ dnorm(0, 1/5)
    b[1] ~ dnorm(0, taub)
    for(i in 2:ndays){
        b[i] ~ dnorm(phi[1]*b[i-1], taub)
    }

    
    for(p in 1:nxcovar){
		betax[p] ~ dnorm(0, 1/4)
	}
    
    #for(p in 1:nNcovar){
    #    betan[p] ~ dnorm(0, 1/4)
    #}
    
    for(r in 1:6){
        gammai[r] <- 0 #~ dnorm(hypgamma, 1/5)
        mux[r] ~ dnorm(hypmux, 1/5)
        lambda[r] ~ dnorm(hyplambda, 1/5)
        taux[r] <- pow(sigma_x[r], -2)
        sigma_x[r] ~ dgamma(8, 2)
        roundprob[r] ~ dbeta(2,2)
        
    }
	hypgamma ~ dnorm(0, 1/5)
	hypmux ~ dnorm(0, 1/5)
	hyplambda ~ dnorm(0, 1/5)
}
"}



# Mixture Models -----

# Province-wide splines
# Prov-wide re
# Indep re's
# Region-Specific beta
# No hyperparameters
mix_jags_base <- {"
model{
    # For the ones trick
    # C chosen so that Lik/c <= 1
    # P(ones == 1) = p, where p is maximized
    # The trick: p is the actual likelihood, we don't care about bernoulli!
    C <- 1000000
    for(d in 1:ndays){
        b[d] ~ dnorm(0, taub)
        
        # Seasonal Affect
        lambda_spline[d] <- inprod(B[d,], lambeta[])
        # Presence spline
        pi_spline[d] <- inprod(B[d,], pibeta[])
        # Dispersion Spline
        alpha_spline[d] <- exp(inprod(B[d,], abeta[]))
    }
    
    for(i in 1:idays){
        # Ones trick
        ones[i] ~ dbern(p[i])
        p[i] <- Lik[i]/C
        Lik[i] <- exp(logLik[i])
        logLik[i] <- (1-z[i])*log(1-w[i]) + 
        #    z[i]*(log(w[i]) + logTPois[i])
            z[i]*(log(w[i]) + LogTruncNB[i])
        
        # For Hurdle
        logit(w[i]) <- pi_spline[DayYearN[i]]
        
        # Truncated Poisson Distribution
        #logTPois[i] <- log(w[i]) + Ni[i]*log(mu[i]) -
        #    mu[i] - log(1 - exp(-mu[i])) - loggam(Ni[i] + 1)
            
        # https://www.cambridge.org/core/services/aop-cambridge-core/content/view/C1FFB63EC7657CD424CE6726F32BC4FE/9781316459515c7_p184-214_CBO.pdf/glms_part_iii_zeroinflated_and_hurdle_models.pdf
        # https://ncss-wpengine.netdna-ssl.com/wp-content/themes/ncss/pdf/Procedures/NCSS/Zero-Inflated_Negative_Binomial_Regression.pdf
        LogTruncNB[i] <-  (1/alpha_spline[DayYearN[i]])*log(u[i])+
            Ni[i]*log(1 - u[i]) + 
            loggam(Ni[i] + 1/alpha_spline[DayYearN[i]]) -
            loggam(1/alpha_spline[DayYearN[i]]) - 
            loggam(Ni[i] + 1) -
            log(1 - (1 + alpha_spline[DayYearN[i]]*mu[i])^(-1/alpha_spline[DayYearN[i]]))
         u[i] <- 1/(1 + alpha_spline[DayYearN[i]]*mu[i])
         
        #https://data.princeton.edu/wws509/notes/countmoments
        ETrunc1[i] <- 1 + alpha_spline[DayYearN[i]]*mu[i]
        ETrunc2[i] <- -1/alpha_spline[DayYearN[i]]
        ETrunc3[i] <- ETrunc1[i]^ETrunc2[i]
        ETruncNB[i] <- mu[i]/(1 - ETrunc3[i])
        #ExTruncNB[i] <- mu[i]/((1-(1 + 
        #    alpha_spline[DayYearN[i]]*mu[i]))^(-1/alpha_spline[DayYearN[i]]))
        ENi[i] <- w[i]*ETruncNB[i]
        
        # mu is the mean of the Poisson
        # lambda[r] is region-specific intercept
        log(mu[i]) <- lambda_spline[DayYearN[i]] + lambda[FireCentre2N[i]] + 
            b[DayYearN[i]] #+ inprod(ncovars, betan)
    }
    
    
    for(j in 1:jfires){
        Xijr_star[j] ~ dlnorm(lp2[j], taux[FireCentre2[j]])
        lp1[j] <- inprod(xcovars[j,], betax[, FireCentre2[j]])
        lp2[j] <- mux[FireCentre2[j]] + 
                lp1[j] + 
                gammai[FireCentre2[j]] * b[DayYear[j]]
        
        # px is probabilty of wider regime (digit preference rather than rounding)
        px[j] ~ dbern(roundprob[FireCentre2[j]])
        ind[j, 1] <- dint[j, 2*px[j] + 1]
        ind[j, 2] <- dint[j, 2*px[j] + 2]
            
        Xijr1[j] ~ dinterval(Xijr_star[j], ind[j,])
        
        EXi[j] <- exp(lp2[j] + 1/(2*taux[FireCentre2[j]]))
        XLik[j] <- log(plnorm(ind[j, 2], lp2[j], taux[FireCentre2[j]]) -
            plnorm(ind[j, 1], lp2[j], taux[FireCentre2[j]]))
    }
    
    # Priors
    taub <- pow(sigma_b, -2)
    sigma_b ~ dgamma(8, 2)
    
    # Splines
    lambeta ~ dmnorm(zero_vec, lambeta_posdef)
    pibeta ~ dmnorm(zero_vec, lambeta_posdef)
    abeta ~ dmnorm(zero_vec, lambeta_posdef)
    

    
    for(p in 1:nxcovar){
        for(r in 1:6){
            betax[p, r] ~ dnorm(0, 1/4)
        }
    }
    
    #for(p in 1:nNcovar){
    #    betan[p] ~ dnorm(0, 1/4)
    #}
    
    for(r in 1:6){
        gammai[r] ~ dnorm(0, 1/5)
        mux[r] ~ dnorm(0, 1/5)
        lambda[r] ~ dnorm(0, 1/5)
        taux[r] <- pow(sigma_x[r], -2)
        sigma_x[r] ~ dgamma(8, 2)
        roundprob[r] ~ dbeta(2,2)
        
    }
}
"}

# Province-wide splines
# Prov-wide re
# Indep re's
# Region-Specific beta
# No hyperparameters
# Special: var selection
mix_jags_varsel <- {"
model{
    # For the ones trick
    # C chosen so that Lik/c <= 1
    # P(ones == 1) = p, where p is maximized
    # The trick: p is the actual likelihood, we don't care about bernoulli!
    C <- 1000000
    for(d in 1:ndays){
        b[d] ~ dnorm(0, taub)
        
        # Seasonal Affect
        lambda_spline[d] <- inprod(B[d,], lambeta[])
        # Presence spline
        pi_spline[d] <- inprod(B[d,], pibeta[])
        # Dispersion Spline
        alpha_spline[d] <- exp(inprod(B[d,], abeta[]))
    }
    
    for(i in 1:idays){
        # Ones trick
        ones[i] ~ dbern(p[i])
        p[i] <- Lik[i]/C
        Lik[i] <- exp(logLik[i])
        logLik[i] <- (1-z[i])*log(1-w[i]) + 
        #    z[i]*(log(w[i]) + logTPois[i])
            z[i]*(log(w[i]) + LogTruncNB[i])
        
        # For Hurdle
        logit(w[i]) <- pi_spline[DayYearN[i]]
        
        # Truncated Poisson Distribution
        #logTPois[i] <- log(w[i]) + Ni[i]*log(mu[i]) -
        #    mu[i] - log(1 - exp(-mu[i])) - loggam(Ni[i] + 1)
            
        # https://www.cambridge.org/core/services/aop-cambridge-core/content/view/C1FFB63EC7657CD424CE6726F32BC4FE/9781316459515c7_p184-214_CBO.pdf/glms_part_iii_zeroinflated_and_hurdle_models.pdf
        # https://ncss-wpengine.netdna-ssl.com/wp-content/themes/ncss/pdf/Procedures/NCSS/Zero-Inflated_Negative_Binomial_Regression.pdf
        LogTruncNB[i] <-  (1/alpha_spline[DayYearN[i]])*log(u[i])+
            Ni[i]*log(1 - u[i]) + 
            loggam(Ni[i] + 1/alpha_spline[DayYearN[i]]) -
            loggam(1/alpha_spline[DayYearN[i]]) - 
            loggam(Ni[i] + 1) -
            log(1 - (1 + alpha_spline[DayYearN[i]]*mu[i])^(-1/alpha_spline[DayYearN[i]]))
         u[i] <- 1/(1 + alpha_spline[DayYearN[i]]*mu[i])
         
        #https://data.princeton.edu/wws509/notes/countmoments
        ETrunc1[i] <- 1 + alpha_spline[DayYearN[i]]*mu[i]
        ETrunc2[i] <- -1/alpha_spline[DayYearN[i]]
        ETrunc3[i] <- ETrunc1[i]^ETrunc2[i]
        ETruncNB[i] <- mu[i]/(1 - ETrunc3[i])
        #ExTruncNB[i] <- mu[i]/((1-(1 + 
        #    alpha_spline[DayYearN[i]]*mu[i]))^(-1/alpha_spline[DayYearN[i]]))
        ENi[i] <- w[i]*ETruncNB[i]
        
        # mu is the mean of the Poisson
        # lambda[r] is region-specific intercept
        log(mu[i]) <- lambda_spline[DayYearN[i]] + lambda[FireCentre2N[i]] + 
            b[DayYearN[i]] #+ inprod(ncovars, betan)
    }
    
    
    for(j in 1:jfires){
        Xijr_star[j] ~ dlnorm(lp2[j], taux[FireCentre2[j]])
        lp1[j] <- inprod(xcovars[j,], betax[, FireCentre2[j]])
        lp2[j] <- mux[FireCentre2[j]] + 
                lp1[j] + 
                gammai[FireCentre2[j]] * b[DayYear[j]]
        
        # px is probabilty of wider regime (digit preference rather than rounding)
        px[j] ~ dbern(roundprob[FireCentre2[j]])
        ind[j, 1] <- dint[j, 2*px[j] + 1]
        ind[j, 2] <- dint[j, 2*px[j] + 2]
            
        Xijr1[j] ~ dinterval(Xijr_star[j], ind[j,])
        
        EXi[j] <- exp(lp2[j] + 1/(2*taux[FireCentre2[j]]))
        XLik[j] <- log(plnorm(ind[j, 2], lp2[j], taux[FireCentre2[j]]) -
            plnorm(ind[j, 1], lp2[j], taux[FireCentre2[j]]))
    }
    
    # Priors
    taub <- pow(sigma_b, -2)
    sigma_b ~ dgamma(8, 2)
    
    # Splines
    lambeta ~ dmnorm(zero_vec, lambeta_posdef)
    pibeta ~ dmnorm(zero_vec, lambeta_posdef)
    abeta ~ dmnorm(zero_vec, lambeta_posdef)
    

    pind ~ dbeta(2,8)
    for(p in 1:nxcovar){
		bind[p] ~ dbern(pind)
        for(r in 1:6){
            betaxT[p, r] ~ dnorm(0, 1/4)
            betax[p, r] <- bind[p]*betaxT[p, r]
        }
    }
    
    #for(p in 1:nNcovar){
    #    betan[p] ~ dnorm(0, 1/4)
    #}
    
    for(r in 1:6){
        gammai[r] ~ dnorm(0, 1/5)
        mux[r] ~ dnorm(0, 1/5)
        lambda[r] ~ dnorm(0, 1/5)
        taux[r] <- pow(sigma_x[r], -2)
        sigma_x[r] ~ dgamma(8, 2)
        roundprob[r] ~ dbeta(2,2)
        
    }
}
"}

# Province-wide splines
# Prov-wide re
# Indep re's
# Region-Specific beta
# Yes hyperparameters
mix_jags_hyper <- {"
model{
    # For the ones trick
    # C chosen so that Lik/c <= 1
    # P(ones == 1) = p, where p is maximized
    # The trick: p is the actual likelihood, we don't care about bernoulli!
    C <- 1000000
    for(d in 1:ndays){
        b[d] ~ dnorm(0, taub)
        
        # Seasonal Affect
        lambda_spline[d] <- inprod(B[d,], lambeta[])
        # Presence spline
        pi_spline[d] <- inprod(B[d,], pibeta[])
        # Dispersion Spline
        alpha_spline[d] <- exp(inprod(B[d,], abeta[]))
    }
    
    for(i in 1:idays){
        # Ones trick
        ones[i] ~ dbern(p[i])
        p[i] <- Lik[i]/C
        Lik[i] <- exp(logLik[i])
        logLik[i] <- (1-z[i])*log(1-w[i]) + 
        #    z[i]*(log(w[i]) + logTPois[i])
            z[i]*(log(w[i]) + LogTruncNB[i])
        
        # For Hurdle
        logit(w[i]) <- pi_spline[DayYearN[i]]
        
        # Truncated Poisson Distribution
        #logTPois[i] <- log(w[i]) + Ni[i]*log(mu[i]) -
        #    mu[i] - log(1 - exp(-mu[i])) - loggam(Ni[i] + 1)
            
        # https://www.cambridge.org/core/services/aop-cambridge-core/content/view/C1FFB63EC7657CD424CE6726F32BC4FE/9781316459515c7_p184-214_CBO.pdf/glms_part_iii_zeroinflated_and_hurdle_models.pdf
        # https://ncss-wpengine.netdna-ssl.com/wp-content/themes/ncss/pdf/Procedures/NCSS/Zero-Inflated_Negative_Binomial_Regression.pdf
        LogTruncNB[i] <-  (1/alpha_spline[DayYearN[i]])*log(u[i])+
            Ni[i]*log(1 - u[i]) + 
            loggam(Ni[i] + 1/alpha_spline[DayYearN[i]]) -
            loggam(1/alpha_spline[DayYearN[i]]) - 
            loggam(Ni[i] + 1) -
            log(1 - (1 + alpha_spline[DayYearN[i]]*mu[i])^(-1/alpha_spline[DayYearN[i]]))
         u[i] <- 1/(1 + alpha_spline[DayYearN[i]]*mu[i])
         
        #https://data.princeton.edu/wws509/notes/countmoments
        ETrunc1[i] <- 1 + alpha_spline[DayYearN[i]]*mu[i]
        ETrunc2[i] <- -1/alpha_spline[DayYearN[i]]
        ETrunc3[i] <- ETrunc1[i]^ETrunc2[i]
        ETruncNB[i] <- mu[i]/(1 - ETrunc3[i])
        #ExTruncNB[i] <- mu[i]/((1-(1 + 
        #    alpha_spline[DayYearN[i]]*mu[i]))^(-1/alpha_spline[DayYearN[i]]))
        ENi[i] <- w[i]*ETruncNB[i]
        
        # mu is the mean of the Poisson
        # lambda[r] is region-specific intercept
        log(mu[i]) <- lambda_spline[DayYearN[i]] + lambda[FireCentre2N[i]] + 
            b[DayYearN[i]] #+ inprod(ncovars, betan)
    }
    
    
    for(j in 1:jfires){
        Xijr_star[j] ~ dlnorm(lp2[j], taux[FireCentre2[j]])
        lp1[j] <- inprod(xcovars[j,], betax[, FireCentre2[j]])
        lp2[j] <- mux[FireCentre2[j]] + 
                lp1[j] + 
                gammai[FireCentre2[j]] * b[DayYear[j]]
        
        # px is probabilty of wider regime (digit preference rather than rounding)
        px[j] ~ dbern(roundprob[FireCentre2[j]])
        ind[j, 1] <- dint[j, 2*px[j] + 1]
        ind[j, 2] <- dint[j, 2*px[j] + 2]
            
        Xijr1[j] ~ dinterval(Xijr_star[j], ind[j,])
        
        EXi[j] <- exp(lp2[j] + 1/(2*taux[FireCentre2[j]]))
        XLik[j] <- log(plnorm(ind[j, 2], lp2[j], taux[FireCentre2[j]]) -
            plnorm(ind[j, 1], lp2[j], taux[FireCentre2[j]]))
    }
    
    # Priors
    taub <- pow(sigma_b, -2)
    sigma_b ~ dgamma(8, 2)
    
    # Splines
    lambeta ~ dmnorm(zero_vec, lambeta_posdef)
    pibeta ~ dmnorm(zero_vec, lambeta_posdef)
    abeta ~ dmnorm(zero_vec, lambeta_posdef)
    

    
    for(p in 1:nxcovar){
        for(r in 1:6){
            betax[p, r] ~ dnorm(0, 1/4)
        }
    }
    
    #for(p in 1:nNcovar){
    #    betan[p] ~ dnorm(0, 1/4)
    #}
    
    for(r in 1:6){
        gammai[r] ~ dnorm(hypgamma, 1/5)
        mux[r] ~ dnorm(hypmux, 1/5)
        lambda[r] ~ dnorm(hyplambda, 1/5)
        taux[r] <- pow(sigma_x[r], -2)
        sigma_x[r] ~ dgamma(8, 2)
        roundprob[r] ~ dbeta(2,2)
        
    }
	hypgamma ~ dnorm(0, 1/5)
	hypmux ~ dnorm(0, 1/5)
	hyplambda ~ dnorm(0, 1/5)
}
"}

# Province-wide splines
# Prov-wide re
# AR1 re's
# Region-Specific beta
# Yes hyperparameters
mix_jags_AR1 <- {"
model{
    # For the ones trick
    # C chosen so that Lik/c <= 1
    # P(ones == 1) = p, where p is maximized
    # The trick: p is the actual likelihood, we don't care about bernoulli!
    C <- 1000000
    for(d in 1:ndays){
        # Seasonal Affect
        lambda_spline[d] <- inprod(B[d,], lambeta[])
        # Presence spline
        pi_spline[d] <- inprod(B[d,], pibeta[])
        # Dispersion Spline
        alpha_spline[d] <- exp(inprod(B[d,], abeta[]))
    }
    
    for(i in 1:idays){
        # Ones trick
        ones[i] ~ dbern(p[i])
        p[i] <- Lik[i]/C
        Lik[i] <- exp(logLik[i])
        logLik[i] <- (1-z[i])*log(1-w[i]) + 
        #    z[i]*(log(w[i]) + logTPois[i])
            z[i]*(log(w[i]) + LogTruncNB[i])
        
        # For Hurdle
        logit(w[i]) <- pi_spline[DayYearN[i]]
        
        # Truncated Poisson Distribution
        #logTPois[i] <- log(w[i]) + Ni[i]*log(mu[i]) -
        #    mu[i] - log(1 - exp(-mu[i])) - loggam(Ni[i] + 1)
            
        # https://www.cambridge.org/core/services/aop-cambridge-core/content/view/C1FFB63EC7657CD424CE6726F32BC4FE/9781316459515c7_p184-214_CBO.pdf/glms_part_iii_zeroinflated_and_hurdle_models.pdf
        # https://ncss-wpengine.netdna-ssl.com/wp-content/themes/ncss/pdf/Procedures/NCSS/Zero-Inflated_Negative_Binomial_Regression.pdf
        LogTruncNB[i] <-  (1/alpha_spline[DayYearN[i]])*log(u[i])+
            Ni[i]*log(1 - u[i]) + 
            loggam(Ni[i] + 1/alpha_spline[DayYearN[i]]) -
            loggam(1/alpha_spline[DayYearN[i]]) - 
            loggam(Ni[i] + 1) -
            log(1 - (1 + alpha_spline[DayYearN[i]]*mu[i])^(-1/alpha_spline[DayYearN[i]]))
         u[i] <- 1/(1 + alpha_spline[DayYearN[i]]*mu[i])
         
        #https://data.princeton.edu/wws509/notes/countmoments
        ETrunc1[i] <- 1 + alpha_spline[DayYearN[i]]*mu[i]
        ETrunc2[i] <- -1/alpha_spline[DayYearN[i]]
        ETrunc3[i] <- ETrunc1[i]^ETrunc2[i]
        ETruncNB[i] <- mu[i]/(1 - ETrunc3[i])
        #ExTruncNB[i] <- mu[i]/((1-(1 + 
        #    alpha_spline[DayYearN[i]]*mu[i]))^(-1/alpha_spline[DayYearN[i]]))
        ENi[i] <- w[i]*ETruncNB[i]
        
        # mu is the mean of the Poisson
        # lambda[r] is region-specific intercept
        log(mu[i]) <- lambda_spline[DayYearN[i]] + lambda[FireCentre2N[i]] + 
            b[DayYearN[i]] #+ inprod(ncovars, betan)
    }
    
    
    for(j in 1:jfires){
        Xijr_star[j] ~ dlnorm(lp2[j], taux[FireCentre2[j]])
        lp1[j] <- inprod(xcovars[j,], betax[, FireCentre2[j]])
        lp2[j] <- mux[FireCentre2[j]] + 
                lp1[j] + 
                gammai[FireCentre2[j]] * b[DayYear[j]]
        
        # px is probabilty of wider regime (digit preference rather than rounding)
        px[j] ~ dbern(roundprob[FireCentre2[j]])
        ind[j, 1] <- dint[j, 2*px[j] + 1]
        ind[j, 2] <- dint[j, 2*px[j] + 2]
            
        Xijr1[j] ~ dinterval(Xijr_star[j], ind[j,])
        
        EXi[j] <- exp(lp2[j] + 1/(2*taux[FireCentre2[j]]))
        XLik[j] <- log(plnorm(ind[j, 2], lp2[j], taux[FireCentre2[j]]) -
            plnorm(ind[j, 1], lp2[j], taux[FireCentre2[j]]))
    }
    
    # Priors
    taub <- pow(sigma_b, -2)
    sigma_b ~ dgamma(8, 2)
    
    # Splines
    lambeta ~ dmnorm(zero_vec, lambeta_posdef)
    pibeta ~ dmnorm(zero_vec, lambeta_posdef)
    abeta ~ dmnorm(zero_vec, lambeta_posdef)
    
    # AR(1) structure for random effects
    phi[1] ~ dnorm(0, 1/5)
    b[1] ~ dnorm(0, taub)
    for(i in 2:ndays){
        b[i] ~ dnorm(phi[1]*b[i-1], taub)
    }

    
    for(p in 1:nxcovar){
        for(r in 1:6){
            betax[p, r] ~ dnorm(0, 1/4)
        }
    }
    
    #for(p in 1:nNcovar){
    #    betan[p] ~ dnorm(0, 1/4)
    #}
    
    for(r in 1:6){
        gammai[r] ~ dnorm(hypgamma, 1/5)
        mux[r] ~ dnorm(hypmux, 1/5)
        lambda[r] ~ dnorm(hyplambda, 1/5)
        taux[r] <- pow(sigma_x[r], -2)
        sigma_x[r] ~ dgamma(8, 2)
        roundprob[r] ~ dbeta(2,2)
        
    }
	hypgamma ~ dnorm(0, 1/5)
	hypmux ~ dnorm(0, 1/5)
	hyplambda ~ dnorm(0, 1/5)
}
"}

# Province-wide splines
# Prov-wide re
# AR2 re's
# Region-Specific beta
# Yes hyperparameters
mix_jags_AR2 <- {"
model{
    # For the ones trick
    # C chosen so that Lik/c <= 1
    # P(ones == 1) = p, where p is maximized
    # The trick: p is the actual likelihood, we don't care about bernoulli!
    C <- 1000000
    for(d in 1:ndays){
        # Seasonal Affect
        lambda_spline[d] <- inprod(B[d,], lambeta[])
        # Presence spline
        pi_spline[d] <- inprod(B[d,], pibeta[])
        # Dispersion Spline
        alpha_spline[d] <- exp(inprod(B[d,], abeta[]))
    }
    
    for(i in 1:idays){
        # Ones trick
        ones[i] ~ dbern(p[i])
        p[i] <- Lik[i]/C
        Lik[i] <- exp(logLik[i])
        logLik[i] <- (1-z[i])*log(1-w[i]) + 
        #    z[i]*(log(w[i]) + logTPois[i])
            z[i]*(log(w[i]) + LogTruncNB[i])
        
        # For Hurdle
        logit(w[i]) <- pi_spline[DayYearN[i]]
        
        # Truncated Poisson Distribution
        #logTPois[i] <- log(w[i]) + Ni[i]*log(mu[i]) -
        #    mu[i] - log(1 - exp(-mu[i])) - loggam(Ni[i] + 1)
            
        # https://www.cambridge.org/core/services/aop-cambridge-core/content/view/C1FFB63EC7657CD424CE6726F32BC4FE/9781316459515c7_p184-214_CBO.pdf/glms_part_iii_zeroinflated_and_hurdle_models.pdf
        # https://ncss-wpengine.netdna-ssl.com/wp-content/themes/ncss/pdf/Procedures/NCSS/Zero-Inflated_Negative_Binomial_Regression.pdf
        LogTruncNB[i] <-  (1/alpha_spline[DayYearN[i]])*log(u[i])+
            Ni[i]*log(1 - u[i]) + 
            loggam(Ni[i] + 1/alpha_spline[DayYearN[i]]) -
            loggam(1/alpha_spline[DayYearN[i]]) - 
            loggam(Ni[i] + 1) -
            log(1 - (1 + alpha_spline[DayYearN[i]]*mu[i])^(-1/alpha_spline[DayYearN[i]]))
         u[i] <- 1/(1 + alpha_spline[DayYearN[i]]*mu[i])
         
        #https://data.princeton.edu/wws509/notes/countmoments
        ETrunc1[i] <- 1 + alpha_spline[DayYearN[i]]*mu[i]
        ETrunc2[i] <- -1/alpha_spline[DayYearN[i]]
        ETrunc3[i] <- ETrunc1[i]^ETrunc2[i]
        ETruncNB[i] <- mu[i]/(1 - ETrunc3[i])
        #ExTruncNB[i] <- mu[i]/((1-(1 + 
        #    alpha_spline[DayYearN[i]]*mu[i]))^(-1/alpha_spline[DayYearN[i]]))
        ENi[i] <- w[i]*ETruncNB[i]
        
        # mu is the mean of the Poisson
        # lambda[r] is region-specific intercept
        log(mu[i]) <- lambda_spline[DayYearN[i]] + lambda[FireCentre2N[i]] + 
            b[DayYearN[i]] #+ inprod(ncovars, betan)
    }
    
    
    for(j in 1:jfires){
        Xijr_star[j] ~ dlnorm(lp2[j], taux[FireCentre2[j]])
        lp1[j] <- inprod(xcovars[j,], betax[, FireCentre2[j]])
        lp2[j] <- mux[FireCentre2[j]] + 
                lp1[j] + 
                gammai[FireCentre2[j]] * b[DayYear[j]]
        
        # px is probabilty of wider regime (digit preference rather than rounding)
        px[j] ~ dbern(roundprob[FireCentre2[j]])
        ind[j, 1] <- dint[j, 2*px[j] + 1]
        ind[j, 2] <- dint[j, 2*px[j] + 2]
            
        Xijr1[j] ~ dinterval(Xijr_star[j], ind[j,])
        
        EXi[j] <- exp(lp2[j] + 1/(2*taux[FireCentre2[j]]))
        XLik[j] <- log(plnorm(ind[j, 2], lp2[j], taux[FireCentre2[j]]) -
            plnorm(ind[j, 1], lp2[j], taux[FireCentre2[j]]))
    }
    
    # Priors
    taub <- pow(sigma_b, -2)
    sigma_b ~ dgamma(8, 2)
    
    # Splines
    lambeta ~ dmnorm(zero_vec, lambeta_posdef)
    pibeta ~ dmnorm(zero_vec, lambeta_posdef)
    abeta ~ dmnorm(zero_vec, lambeta_posdef)
    
    # AR(2) structure for random effects
    phi[1] ~ dnorm(0, 1/5)
    phi[2] ~ dnorm(0, 1/5)
    b[1] ~ dnorm(0, taub)
    b[2] ~ dnorm(phi[1]*b[1], taub)
    for(i in 3:idays){
        b[i] ~ dnorm(phi[1]*b[i-1] + phi[2]*b[i-2], taub)
    }

    
    for(p in 1:nxcovar){
        for(r in 1:6){
            betax[p, r] ~ dnorm(0, 1/4)
        }
    }
    
    #for(p in 1:nNcovar){
    #    betan[p] ~ dnorm(0, 1/4)
    #}
    
    for(r in 1:6){
        gammai[r] ~ dnorm(hypgamma, 1/5)
        mux[r] ~ dnorm(hypmux, 1/5)
        lambda[r] ~ dnorm(hyplambda, 1/5)
        taux[r] <- pow(sigma_x[r], -2)
        sigma_x[r] ~ dgamma(8, 2)
        roundprob[r] ~ dbeta(2,2)
        
    }
	hypgamma ~ dnorm(0, 1/5)
	hypmux ~ dnorm(0, 1/5)
	hyplambda ~ dnorm(0, 1/5)
}
"}

# Province-wide splines
	# alpha is not a spline
# Prov-wide re
# AR1 re's
# Region-Specific beta
# No hyperparameters
mix_jags_daily_alpha <- {"
model{
    # For the ones trick
    # C chosen so that Lik/c <= 1
    # P(ones == 1) = p, where p is maximized
    # The trick: p is the actual likelihood, we don't care about bernoulli!
    C <- 1000000
    for(d in 1:ndays){
        b[d] ~ dnorm(0, taub)
        
        # Seasonal Affect
        lambda_spline[d] <- inprod(B[d,], lambeta[])
        # Presence spline
        pi_spline[d] <- inprod(B[d,], pibeta[])
        # Dispersion Spline
        #alpha_spline[d] <- exp(inprod(B[d,], abeta[]))
		alpha[d] ~ dunif(0, 100)
    }
    
    for(i in 1:idays){
        # Ones trick
        ones[i] ~ dbern(p[i])
        p[i] <- Lik[i]/C
        Lik[i] <- exp(logLik[i])
        logLik[i] <- (1-z[i])*log(1-w[i]) + 
        #    z[i]*(log(w[i]) + logTPois[i])
            z[i]*(log(w[i]) + LogTruncNB[i])
        
        # For Hurdle
        logit(w[i]) <- pi_spline[DayYearN[i]]
        
        # Truncated Poisson Distribution
        #logTPois[i] <- log(w[i]) + Ni[i]*log(mu[i]) -
        #    mu[i] - log(1 - exp(-mu[i])) - loggam(Ni[i] + 1)
            
        # https://www.cambridge.org/core/services/aop-cambridge-core/content/view/C1FFB63EC7657CD424CE6726F32BC4FE/9781316459515c7_p184-214_CBO.pdf/glms_part_iii_zeroinflated_and_hurdle_models.pdf
        # https://ncss-wpengine.netdna-ssl.com/wp-content/themes/ncss/pdf/Procedures/NCSS/Zero-Inflated_Negative_Binomial_Regression.pdf
        LogTruncNB[i] <-  (1/alpha[DayYearN[i]])*log(u[i])+
            Ni[i]*log(1 - u[i]) + 
            loggam(Ni[i] + 1/alpha[DayYearN[i]]) -
            loggam(1/alpha[DayYearN[i]]) - 
            loggam(Ni[i] + 1) -
            log(1 - (1 + alpha[DayYearN[i]]*mu[i])^(-1/alpha[DayYearN[i]]))
         u[i] <- 1/(1 + alpha[DayYearN[i]]*mu[i])
         
        #https://data.princeton.edu/wws509/notes/countmoments
        ETrunc1[i] <- 1 + alpha[DayYearN[i]]*mu[i]
        ETrunc2[i] <- -1/alpha[DayYearN[i]]
        ETrunc3[i] <- ETrunc1[i]^ETrunc2[i]
        ETruncNB[i] <- mu[i]/(1 - ETrunc3[i])
        #ExTruncNB[i] <- mu[i]/((1-(1 + 
        #    alpha[DayYearN[i]]*mu[i]))^(-1/alpha[DayYearN[i]]))
        ENi[i] <- w[i]*ETruncNB[i]
        
        # mu is the mean of the Poisson
        # lambda[r] is region-specific intercept
        log(mu[i]) <- lambda_spline[DayYearN[i]] + lambda[FireCentre2N[i]] + 
            b[DayYearN[i]] #+ inprod(ncovars, betan)
    }
    
    
    for(j in 1:jfires){
        Xijr_star[j] ~ dlnorm(lp2[j], taux[FireCentre2[j]])
        lp1[j] <- inprod(xcovars[j,], betax[, FireCentre2[j]])
        lp2[j] <- mux[FireCentre2[j]] + 
                lp1[j] + 
                gammai[FireCentre2[j]] * b[DayYear[j]]
        
        # px is probabilty of wider regime (digit preference rather than rounding)
        px[j] ~ dbern(roundprob[FireCentre2[j]])
        ind[j, 1] <- dint[j, 2*px[j] + 1]
        ind[j, 2] <- dint[j, 2*px[j] + 2]
            
        Xijr1[j] ~ dinterval(Xijr_star[j], ind[j,])
        
        EXi[j] <- exp(lp2[j] + 1/(2*taux[FireCentre2[j]]))
        XLik[j] <- log(plnorm(ind[j, 2], lp2[j], taux[FireCentre2[j]]) -
            plnorm(ind[j, 1], lp2[j], taux[FireCentre2[j]]))
    }
    
    # Priors
    taub <- pow(sigma_b, -2)
    sigma_b ~ dgamma(8, 2)
    
    # Splines
    lambeta ~ dmnorm(zero_vec, lambeta_posdef)
    pibeta ~ dmnorm(zero_vec, lambeta_posdef)
    #abeta ~ dmnorm(zero_vec, lambeta_posdef)
    

    
    for(p in 1:nxcovar){
        for(r in 1:6){
            betax[p, r] ~ dnorm(0, 1/4)
        }
    }
    
    #for(p in 1:nNcovar){
    #    betan[p] ~ dnorm(0, 1/4)
    #}
    
    for(r in 1:6){
        gammai[r] ~ dnorm(0, 1/5)
        mux[r] ~ dnorm(0, 1/5)
        lambda[r] ~ dnorm(0, 1/5)
        taux[r] <- pow(sigma_x[r], -2)
        sigma_x[r] ~ dgamma(8, 2)
        roundprob[r] ~ dbeta(2,2)
        
    }
}
"}

# Province-wide splines
# Region-specific re
# Indep re's
# Region-Specific beta
# No hyperparameters
mix_jags_regionre <- {"
model{
    # For the ones trick
    # C chosen so that Lik/c <= 1
    # P(ones == 1) = p, where p is maximized
    # The trick: p is the actual likelihood, we don't care about bernoulli!
    C <- 1000000
    for(d in 1:ndays){
		for(r in 1:6){
			b[r, d] ~ dnorm(0, taub)
		}
        
        # Seasonal Affect
        lambda_spline[d] <- inprod(B[d,], lambeta[])
        # Presence spline
        pi_spline[d] <- inprod(B[d,], pibeta[])
        # Dispersion Spline
        alpha_spline[d] <- exp(inprod(B[d,], abeta[]))
    }
    
    for(i in 1:idays){
        # Ones trick
        ones[i] ~ dbern(p[i])
        p[i] <- Lik[i]/C
        Lik[i] <- exp(logLik[i])
        logLik[i] <- (1-z[i])*log(1-w[i]) + 
        #    z[i]*(log(w[i]) + logTPois[i])
            z[i]*(log(w[i]) + LogTruncNB[i])
        
        # For Hurdle
        logit(w[i]) <- pi_spline[DayYearN[i]]
        
        # Truncated Poisson Distribution
        #logTPois[i] <- log(w[i]) + Ni[i]*log(mu[i]) -
        #    mu[i] - log(1 - exp(-mu[i])) - loggam(Ni[i] + 1)
            
        # https://www.cambridge.org/core/services/aop-cambridge-core/content/view/C1FFB63EC7657CD424CE6726F32BC4FE/9781316459515c7_p184-214_CBO.pdf/glms_part_iii_zeroinflated_and_hurdle_models.pdf
        # https://ncss-wpengine.netdna-ssl.com/wp-content/themes/ncss/pdf/Procedures/NCSS/Zero-Inflated_Negative_Binomial_Regression.pdf
        LogTruncNB[i] <-  (1/alpha_spline[DayYearN[i]])*log(u[i])+
            Ni[i]*log(1 - u[i]) + 
            loggam(Ni[i] + 1/alpha_spline[DayYearN[i]]) -
            loggam(1/alpha_spline[DayYearN[i]]) - 
            loggam(Ni[i] + 1) -
            log(1 - (1 + alpha_spline[DayYearN[i]]*mu[i])^(-1/alpha_spline[DayYearN[i]]))
         u[i] <- 1/(1 + alpha_spline[DayYearN[i]]*mu[i])
         
        #https://data.princeton.edu/wws509/notes/countmoments
        ETrunc1[i] <- 1 + alpha_spline[DayYearN[i]]*mu[i]
        ETrunc2[i] <- -1/alpha_spline[DayYearN[i]]
        ETrunc3[i] <- ETrunc1[i]^ETrunc2[i]
        ETruncNB[i] <- mu[i]/(1 - ETrunc3[i])
        #ExTruncNB[i] <- mu[i]/((1-(1 + 
        #    alpha_spline[DayYearN[i]]*mu[i]))^(-1/alpha_spline[DayYearN[i]]))
        ENi[i] <- w[i]*ETruncNB[i]
        
        # mu is the mean of the Poisson
        # lambda[r] is region-specific intercept
        log(mu[i]) <- lambda_spline[DayYearN[i]] + lambda[FireCentre2N[i]] + 
            b[FireCentre2N[i], DayYearN[i]] #+ inprod(ncovars, betan)
    }
    
    
    for(j in 1:jfires){
        Xijr_star[j] ~ dlnorm(lp2[j], taux[FireCentre2[j]])
        lp1[j] <- inprod(xcovars[j,], betax[, FireCentre2[j]])
        lp2[j] <- mux[FireCentre2[j]] + 
                lp1[j] + 
                gammai[FireCentre2[j]] * b[FireCentre2[j], DayYear[j]]
        
        # px is probabilty of wider regime (digit preference rather than rounding)
        px[j] ~ dbern(roundprob[FireCentre2[j]])
        ind[j, 1] <- dint[j, 2*px[j] + 1]
        ind[j, 2] <- dint[j, 2*px[j] + 2]
            
        Xijr1[j] ~ dinterval(Xijr_star[j], ind[j,])
        
        EXi[j] <- exp(lp2[j] + 1/(2*taux[FireCentre2[j]]))
        XLik[j] <- log(plnorm(ind[j, 2], lp2[j], taux[FireCentre2[j]]) -
            plnorm(ind[j, 1], lp2[j], taux[FireCentre2[j]]))
    }
    
    # Priors
    taub <- pow(sigma_b, -2)
    sigma_b ~ dgamma(8, 2)
    
    # Splines
    lambeta ~ dmnorm(zero_vec, lambeta_posdef)
    pibeta ~ dmnorm(zero_vec, lambeta_posdef)
    abeta ~ dmnorm(zero_vec, lambeta_posdef)
    

    
    for(p in 1:nxcovar){
        for(r in 1:6){
            betax[p, r] ~ dnorm(0, 1/4)
        }
    }
    
    #for(p in 1:nNcovar){
    #    betan[p] ~ dnorm(0, 1/4)
    #}
    
    for(r in 1:6){
        gammai[r] ~ dnorm(0, 1/5)
        mux[r] ~ dnorm(0, 1/5)
        lambda[r] ~ dnorm(0, 1/5)
        taux[r] <- pow(sigma_x[r], -2)
        sigma_x[r] ~ dgamma(8, 2)
        roundprob[r] ~ dbeta(2,2)
        
    }
}
"}

# Province-wide splines
	# alpha is not a spline
# Prov-wide re
# AR1 re's
# Region-Specific beta
# Yes hyperparameters
mix_jags_AR1dahyp <- {"
model{
    # For the ones trick
    # C chosen so that Lik/c <= 1
    # P(ones == 1) = p, where p is maximized
    # The trick: p is the actual likelihood, we don't care about bernoulli!
    C <- 1000000
    for(d in 1:ndays){
        # Seasonal Affect
        lambda_spline[d] <- inprod(B[d,], lambeta[])
        # Presence spline
        pi_spline[d] <- inprod(B[d,], pibeta[])
        # Dispersion Spline
		alpha[d] ~ dunif(0, 100)
    }
    
    for(i in 1:idays){
        # Ones trick
        ones[i] ~ dbern(p[i])
        p[i] <- Lik[i]/C
        Lik[i] <- exp(logLik[i])
        logLik[i] <- (1-z[i])*log(1-w[i]) + 
        #    z[i]*(log(w[i]) + logTPois[i])
            z[i]*(log(w[i]) + LogTruncNB[i])
        
        # For Hurdle
        logit(w[i]) <- pi_spline[DayYearN[i]]
        
        # Truncated Poisson Distribution
        #logTPois[i] <- log(w[i]) + Ni[i]*log(mu[i]) -
        #    mu[i] - log(1 - exp(-mu[i])) - loggam(Ni[i] + 1)
            
        # https://www.cambridge.org/core/services/aop-cambridge-core/content/view/C1FFB63EC7657CD424CE6726F32BC4FE/9781316459515c7_p184-214_CBO.pdf/glms_part_iii_zeroinflated_and_hurdle_models.pdf
        # https://ncss-wpengine.netdna-ssl.com/wp-content/themes/ncss/pdf/Procedures/NCSS/Zero-Inflated_Negative_Binomial_Regression.pdf
        LogTruncNB[i] <-  (1/alpha[DayYearN[i]])*log(u[i])+
            Ni[i]*log(1 - u[i]) + 
            loggam(Ni[i] + 1/alpha[DayYearN[i]]) -
            loggam(1/alpha[DayYearN[i]]) - 
            loggam(Ni[i] + 1) -
            log(1 - (1 + alpha[DayYearN[i]]*mu[i])^(-1/alpha[DayYearN[i]]))
         u[i] <- 1/(1 + alpha[DayYearN[i]]*mu[i])
         
        #https://data.princeton.edu/wws509/notes/countmoments
        ETrunc1[i] <- 1 + alpha[DayYearN[i]]*mu[i]
        ETrunc2[i] <- -1/alpha[DayYearN[i]]
        ETrunc3[i] <- ETrunc1[i]^ETrunc2[i]
        ETruncNB[i] <- mu[i]/(1 - ETrunc3[i])
        #ExTruncNB[i] <- mu[i]/((1-(1 + 
        #    alpha[DayYearN[i]]*mu[i]))^(-1/alpha[DayYearN[i]]))
        ENi[i] <- w[i]*ETruncNB[i]
        
        # mu is the mean of the Poisson
        # lambda[r] is region-specific intercept
        log(mu[i]) <- lambda_spline[DayYearN[i]] + lambda[FireCentre2N[i]] + 
            b[DayYearN[i]] #+ inprod(ncovars, betan)
    }
    
    
    for(j in 1:jfires){
        Xijr_star[j] ~ dlnorm(lp2[j], taux[FireCentre2[j]])
        lp1[j] <- inprod(xcovars[j,], betax[, FireCentre2[j]])
        lp2[j] <- mux[FireCentre2[j]] + 
                lp1[j] + 
                gammai[FireCentre2[j]] * b[DayYear[j]]
        
        # px is probabilty of wider regime (digit preference rather than rounding)
        px[j] ~ dbern(roundprob[FireCentre2[j]])
        ind[j, 1] <- dint[j, 2*px[j] + 1]
        ind[j, 2] <- dint[j, 2*px[j] + 2]
            
        Xijr1[j] ~ dinterval(Xijr_star[j], ind[j,])
        
        EXi[j] <- exp(lp2[j] + 1/(2*taux[FireCentre2[j]]))
        XLik[j] <- log(plnorm(ind[j, 2], lp2[j], taux[FireCentre2[j]]) -
            plnorm(ind[j, 1], lp2[j], taux[FireCentre2[j]]))
    }
    
    # Priors
    taub <- pow(sigma_b, -2)
    sigma_b ~ dgamma(8, 2)
    
    # Splines
    lambeta ~ dmnorm(zero_vec, lambeta_posdef)
    pibeta ~ dmnorm(zero_vec, lambeta_posdef)
    #abeta ~ dmnorm(zero_vec, lambeta_posdef)
    
    # AR(1) structure for random effects
    phi[1] ~ dnorm(0, 1/5)
    b[1] ~ dnorm(0, taub)
    for(i in 2:ndays){
        b[i] ~ dnorm(phi[1]*b[i-1], taub)
    }

    
    for(p in 1:nxcovar){
        for(r in 1:6){
            betax[p, r] ~ dnorm(0, 1/4)
        }
    }
    
    #for(p in 1:nNcovar){
    #    betan[p] ~ dnorm(0, 1/4)
    #}
    
    for(r in 1:6){
        gammai[r] ~ dnorm(hypgamma, 1/5)
        mux[r] ~ dnorm(hypmux, 1/5)
        lambda[r] ~ dnorm(hyplambda, 1/5)
        taux[r] <- pow(sigma_x[r], -2)
        sigma_x[r] ~ dgamma(8, 2)
        roundprob[r] ~ dbeta(2,2)
        
    }
	hypgamma ~ dnorm(0, 1/5)
	hypmux ~ dnorm(0, 1/5)
	hyplambda ~ dnorm(0, 1/5)
}
"}

# Province-wide splines
# Prov-wide re
# AR1 re's
# Prov-wide beta
# Yes hyperparameters
mix_jags_AR1pbeta <- {"
model{
    # For the ones trick
    # C chosen so that Lik/c <= 1
    # P(ones == 1) = p, where p is maximized
    # The trick: p is the actual likelihood, we don't care about bernoulli!
    C <- 1000000
    for(d in 1:ndays){
        # Seasonal Affect
        lambda_spline[d] <- inprod(B[d,], lambeta[])
        # Presence spline
        pi_spline[d] <- inprod(B[d,], pibeta[])
        # Dispersion Spline
        alpha_spline[d] <- exp(inprod(B[d,], abeta[]))
    }
    
    for(i in 1:idays){
        # Ones trick
        ones[i] ~ dbern(p[i])
        p[i] <- Lik[i]/C
        Lik[i] <- exp(logLik[i])
        logLik[i] <- (1-z[i])*log(1-w[i]) + 
        #    z[i]*(log(w[i]) + logTPois[i])
            z[i]*(log(w[i]) + LogTruncNB[i])
        
        # For Hurdle
        logit(w[i]) <- pi_spline[DayYearN[i]]
        
        # Truncated Poisson Distribution
        #logTPois[i] <- log(w[i]) + Ni[i]*log(mu[i]) -
        #    mu[i] - log(1 - exp(-mu[i])) - loggam(Ni[i] + 1)
            
        # https://www.cambridge.org/core/services/aop-cambridge-core/content/view/C1FFB63EC7657CD424CE6726F32BC4FE/9781316459515c7_p184-214_CBO.pdf/glms_part_iii_zeroinflated_and_hurdle_models.pdf
        # https://ncss-wpengine.netdna-ssl.com/wp-content/themes/ncss/pdf/Procedures/NCSS/Zero-Inflated_Negative_Binomial_Regression.pdf
        LogTruncNB[i] <-  (1/alpha_spline[DayYearN[i]])*log(u[i])+
            Ni[i]*log(1 - u[i]) + 
            loggam(Ni[i] + 1/alpha_spline[DayYearN[i]]) -
            loggam(1/alpha_spline[DayYearN[i]]) - 
            loggam(Ni[i] + 1) -
            log(1 - (1 + alpha_spline[DayYearN[i]]*mu[i])^(-1/alpha_spline[DayYearN[i]]))
         u[i] <- 1/(1 + alpha_spline[DayYearN[i]]*mu[i])
         
        #https://data.princeton.edu/wws509/notes/countmoments
        ETrunc1[i] <- 1 + alpha_spline[DayYearN[i]]*mu[i]
        ETrunc2[i] <- -1/alpha_spline[DayYearN[i]]
        ETrunc3[i] <- ETrunc1[i]^ETrunc2[i]
        ETruncNB[i] <- mu[i]/(1 - ETrunc3[i])
        #ExTruncNB[i] <- mu[i]/((1-(1 + 
        #    alpha_spline[DayYearN[i]]*mu[i]))^(-1/alpha_spline[DayYearN[i]]))
        ENi[i] <- w[i]*ETruncNB[i]
        
        # mu is the mean of the Poisson
        # lambda[r] is region-specific intercept
        log(mu[i]) <- lambda_spline[DayYearN[i]] + lambda[FireCentre2N[i]] + 
            b[DayYearN[i]] #+ inprod(ncovars, betan)
    }
    
    
    for(j in 1:jfires){
        Xijr_star[j] ~ dlnorm(lp2[j], taux[FireCentre2[j]])
        lp1[j] <- inprod(xcovars[j,], betax)
        lp2[j] <- mux[FireCentre2[j]] + 
                lp1[j] + 
                gammai[FireCentre2[j]] * b[DayYear[j]]
        
        # px is probabilty of wider regime (digit preference rather than rounding)
        px[j] ~ dbern(roundprob[FireCentre2[j]])
        ind[j, 1] <- dint[j, 2*px[j] + 1]
        ind[j, 2] <- dint[j, 2*px[j] + 2]
            
        Xijr1[j] ~ dinterval(Xijr_star[j], ind[j,])
        
        EXi[j] <- exp(lp2[j] + 1/(2*taux[FireCentre2[j]]))
        XLik[j] <- log(plnorm(ind[j, 2], lp2[j], taux[FireCentre2[j]]) -
            plnorm(ind[j, 1], lp2[j], taux[FireCentre2[j]]))
    }
    
    # Priors
    taub <- pow(sigma_b, -2)
    sigma_b ~ dgamma(8, 2)
    
    # Splines
    lambeta ~ dmnorm(zero_vec, lambeta_posdef)
    pibeta ~ dmnorm(zero_vec, lambeta_posdef)
    abeta ~ dmnorm(zero_vec, lambeta_posdef)
    
    # AR(1) structure for random effects
    phi[1] ~ dnorm(0, 1/5)
    b[1] ~ dnorm(0, taub)
    for(i in 2:ndays){
        b[i] ~ dnorm(phi[1]*b[i-1], taub)
    }

    
    for(p in 1:nxcovar){
		betax[p] ~ dnorm(0, 1/4)
    }
    
    #for(p in 1:nNcovar){
    #    betan[p] ~ dnorm(0, 1/4)
    #}
    
    for(r in 1:6){
        gammai[r] ~ dnorm(hypgamma, 1/5)
        mux[r] ~ dnorm(hypmux, 1/5)
        lambda[r] ~ dnorm(hyplambda, 1/5)
        taux[r] <- pow(sigma_x[r], -2)
        sigma_x[r] ~ dgamma(8, 2)
        roundprob[r] ~ dbeta(2,2)
        
    }
	hypgamma ~ dnorm(0, 1/5)
	hypmux ~ dnorm(0, 1/5)
	hyplambda ~ dnorm(0, 1/5)
}
"}


# Province-wide splines
# Prov-wide re
# AR1 re's
# Prov-wide beta
# Yes hyperparameters
# gamma set to 0
mix_jags_AR1pbeta_sep <- {"
model{
    # For the ones trick
    # C chosen so that Lik/c <= 1
    # P(ones == 1) = p, where p is maximized
    # The trick: p is the actual likelihood, we don't care about bernoulli!
    C <- 1000000
    for(d in 1:ndays){
        # Seasonal Affect
        lambda_spline[d] <- inprod(B[d,], lambeta[])
        # Presence spline
        pi_spline[d] <- inprod(B[d,], pibeta[])
        # Dispersion Spline
        alpha_spline[d] <- exp(inprod(B[d,], abeta[]))
    }
    
    for(i in 1:idays){
        # Ones trick
        ones[i] ~ dbern(p[i])
        p[i] <- Lik[i]/C
        Lik[i] <- exp(logLik[i])
        logLik[i] <- (1-z[i])*log(1-w[i]) + 
        #    z[i]*(log(w[i]) + logTPois[i])
            z[i]*(log(w[i]) + LogTruncNB[i])
        
        # For Hurdle
        logit(w[i]) <- pi_spline[DayYearN[i]]
        
        # Truncated Poisson Distribution
        #logTPois[i] <- log(w[i]) + Ni[i]*log(mu[i]) -
        #    mu[i] - log(1 - exp(-mu[i])) - loggam(Ni[i] + 1)
            
        # https://www.cambridge.org/core/services/aop-cambridge-core/content/view/C1FFB63EC7657CD424CE6726F32BC4FE/9781316459515c7_p184-214_CBO.pdf/glms_part_iii_zeroinflated_and_hurdle_models.pdf
        # https://ncss-wpengine.netdna-ssl.com/wp-content/themes/ncss/pdf/Procedures/NCSS/Zero-Inflated_Negative_Binomial_Regression.pdf
        LogTruncNB[i] <-  (1/alpha_spline[DayYearN[i]])*log(u[i])+
            Ni[i]*log(1 - u[i]) + 
            loggam(Ni[i] + 1/alpha_spline[DayYearN[i]]) -
            loggam(1/alpha_spline[DayYearN[i]]) - 
            loggam(Ni[i] + 1) -
            log(1 - (1 + alpha_spline[DayYearN[i]]*mu[i])^(-1/alpha_spline[DayYearN[i]]))
         u[i] <- 1/(1 + alpha_spline[DayYearN[i]]*mu[i])
         
        #https://data.princeton.edu/wws509/notes/countmoments
        ETrunc1[i] <- 1 + alpha_spline[DayYearN[i]]*mu[i]
        ETrunc2[i] <- -1/alpha_spline[DayYearN[i]]
        ETrunc3[i] <- ETrunc1[i]^ETrunc2[i]
        ETruncNB[i] <- mu[i]/(1 - ETrunc3[i])
        #ExTruncNB[i] <- mu[i]/((1-(1 + 
        #    alpha_spline[DayYearN[i]]*mu[i]))^(-1/alpha_spline[DayYearN[i]]))
        ENi[i] <- w[i]*ETruncNB[i]
        
        # mu is the mean of the Poisson
        # lambda[r] is region-specific intercept
        log(mu[i]) <- lambda_spline[DayYearN[i]] + lambda[FireCentre2N[i]] + 
            b[DayYearN[i]] #+ inprod(ncovars, betan)
    }
    
    
    for(j in 1:jfires){
        Xijr_star[j] ~ dlnorm(lp2[j], taux[FireCentre2[j]])
        lp1[j] <- inprod(xcovars[j,], betax)
        lp2[j] <- mux[FireCentre2[j]] + 
                lp1[j] + 
                gammai[FireCentre2[j]] * b[DayYear[j]]
        
        # px is probabilty of wider regime (digit preference rather than rounding)
        px[j] ~ dbern(roundprob[FireCentre2[j]])
        ind[j, 1] <- dint[j, 2*px[j] + 1]
        ind[j, 2] <- dint[j, 2*px[j] + 2]
            
        Xijr1[j] ~ dinterval(Xijr_star[j], ind[j,])
        
        EXi[j] <- exp(lp2[j] + 1/(2*taux[FireCentre2[j]]))
        XLik[j] <- log(plnorm(ind[j, 2], lp2[j], taux[FireCentre2[j]]) -
            plnorm(ind[j, 1], lp2[j], taux[FireCentre2[j]]))
    }
    
    # Priors
    taub <- pow(sigma_b, -2)
    sigma_b ~ dgamma(8, 2)
    
    # Splines
    lambeta ~ dmnorm(zero_vec, lambeta_posdef)
    pibeta ~ dmnorm(zero_vec, lambeta_posdef)
    abeta ~ dmnorm(zero_vec, lambeta_posdef)
    
    # AR(1) structure for random effects
    phi[1] ~ dnorm(0, 1/5)
    b[1] ~ dnorm(0, taub)
    for(i in 2:ndays){
        b[i] ~ dnorm(phi[1]*b[i-1], taub)
    }

    
    for(p in 1:nxcovar){
		betax[p] ~ dnorm(0, 1/4)
    }
    
    #for(p in 1:nNcovar){
    #    betan[p] ~ dnorm(0, 1/4)
    #}
    
    for(r in 1:6){
        gammai[r] <- 0#~ dnorm(hypgamma, 1/5)
        mux[r] ~ dnorm(hypmux, 1/5)
        lambda[r] ~ dnorm(hyplambda, 1/5)
        taux[r] <- pow(sigma_x[r], -2)
        sigma_x[r] ~ dgamma(8, 2)
        roundprob[r] ~ dbeta(2,2)
        
    }
	hypgamma ~ dnorm(0, 1/5)
	hypmux ~ dnorm(0, 1/5)
	hyplambda ~ dnorm(0, 1/5)
}
"}

# Province-wide splines
# Prov-wide re
# AR1 re's
# Prov-wide beta
# Yes hyperparameters
# indicator for gamma == 0
mix_jags_AR1pbeta_gamsel <- {"
model{
    # For the ones trick
    # C chosen so that Lik/c <= 1
    # P(ones == 1) = p, where p is maximized
    # The trick: p is the actual likelihood, we don't care about bernoulli!
    C <- 1000000
    for(d in 1:ndays){
        # Seasonal Affect
        lambda_spline[d] <- inprod(B[d,], lambeta[])
        # Presence spline
        pi_spline[d] <- inprod(B[d,], pibeta[])
        # Dispersion Spline
        alpha_spline[d] <- exp(inprod(B[d,], abeta[]))
    }
    
    for(i in 1:idays){
        # Ones trick
        ones[i] ~ dbern(p[i])
        p[i] <- Lik[i]/C
        Lik[i] <- exp(logLik[i])
        logLik[i] <- (1-z[i])*log(1-w[i]) + 
        #    z[i]*(log(w[i]) + logTPois[i])
            z[i]*(log(w[i]) + LogTruncNB[i])
        
        # For Hurdle
        logit(w[i]) <- pi_spline[DayYearN[i]]
        
        # Truncated Poisson Distribution
        #logTPois[i] <- log(w[i]) + Ni[i]*log(mu[i]) -
        #    mu[i] - log(1 - exp(-mu[i])) - loggam(Ni[i] + 1)
            
        # https://www.cambridge.org/core/services/aop-cambridge-core/content/view/C1FFB63EC7657CD424CE6726F32BC4FE/9781316459515c7_p184-214_CBO.pdf/glms_part_iii_zeroinflated_and_hurdle_models.pdf
        # https://ncss-wpengine.netdna-ssl.com/wp-content/themes/ncss/pdf/Procedures/NCSS/Zero-Inflated_Negative_Binomial_Regression.pdf
        LogTruncNB[i] <-  (1/alpha_spline[DayYearN[i]])*log(u[i])+
            Ni[i]*log(1 - u[i]) + 
            loggam(Ni[i] + 1/alpha_spline[DayYearN[i]]) -
            loggam(1/alpha_spline[DayYearN[i]]) - 
            loggam(Ni[i] + 1) -
            log(1 - (1 + alpha_spline[DayYearN[i]]*mu[i])^(-1/alpha_spline[DayYearN[i]]))
         u[i] <- 1/(1 + alpha_spline[DayYearN[i]]*mu[i])
         
        #https://data.princeton.edu/wws509/notes/countmoments
        ETrunc1[i] <- 1 + alpha_spline[DayYearN[i]]*mu[i]
        ETrunc2[i] <- -1/alpha_spline[DayYearN[i]]
        ETrunc3[i] <- ETrunc1[i]^ETrunc2[i]
        ETruncNB[i] <- mu[i]/(1 - ETrunc3[i])
        #ExTruncNB[i] <- mu[i]/((1-(1 + 
        #    alpha_spline[DayYearN[i]]*mu[i]))^(-1/alpha_spline[DayYearN[i]]))
        ENi[i] <- w[i]*ETruncNB[i]
        
        # mu is the mean of the Poisson
        # lambda[r] is region-specific intercept
        log(mu[i]) <- lambda_spline[DayYearN[i]] + lambda[FireCentre2N[i]] + 
            b[DayYearN[i]] #+ inprod(ncovars, betan)
    }
    
    
    for(j in 1:jfires){
        Xijr_star[j] ~ dlnorm(lp2[j], taux[FireCentre2[j]])
        lp1[j] <- inprod(xcovars[j,], betax)
        lp2[j] <- mux[FireCentre2[j]] + 
                lp1[j] + 
                gammai[FireCentre2[j]] * b[DayYear[j]]
        
        # px is probabilty of wider regime (digit preference rather than rounding)
        px[j] ~ dbern(roundprob[FireCentre2[j]])
        ind[j, 1] <- dint[j, 2*px[j] + 1]
        ind[j, 2] <- dint[j, 2*px[j] + 2]
            
        Xijr1[j] ~ dinterval(Xijr_star[j], ind[j,])
        
        EXi[j] <- exp(lp2[j] + 1/(2*taux[FireCentre2[j]]))
        XLik[j] <- log(plnorm(ind[j, 2], lp2[j], taux[FireCentre2[j]]) -
            plnorm(ind[j, 1], lp2[j], taux[FireCentre2[j]]))
    }
    
    # Priors
    taub <- pow(sigma_b, -2)
    sigma_b ~ dgamma(8, 2)
    
    # Splines
    lambeta ~ dmnorm(zero_vec, lambeta_posdef)
    pibeta ~ dmnorm(zero_vec, lambeta_posdef)
    abeta ~ dmnorm(zero_vec, lambeta_posdef)
    
    # AR(1) structure for random effects
    phi[1] ~ dnorm(0, 1/5)
    b[1] ~ dnorm(0, taub)
    for(i in 2:ndays){
        b[i] ~ dnorm(phi[1]*b[i-1], taub)
    }

    
    for(p in 1:nxcovar){
		betax[p] ~ dnorm(0, 1/4)
    }
    
    #for(p in 1:nNcovar){
    #    betan[p] ~ dnorm(0, 1/4)
    #}
    
	gind ~ dbeta(1,1)
    for(r in 1:6){
		gzero[r] ~ dbern(gind)
        gammai_1[r] ~ dnorm(0, 1/5)
		gammai[r] <- gzero[r]*gammai_1[r]
        mux[r] ~ dnorm(0, 1/5)
        lambda[r] ~ dnorm(0, 1/5)
        taux[r] <- pow(sigma_x[r], -2)
        sigma_x[r] ~ dgamma(8, 2)
        roundprob[r] ~ dbeta(2,2)
        
    }
}
"}


# Spline for gamma?

# Province-wide splines
# Prov-wide re
# AR1 re's
# Prov-wide beta
# Yes hyperparameters
mix_jags_AR1pbeta_gspline <- {"
model{
    # For the ones trick
    # C chosen so that Lik/c <= 1
    # P(ones == 1) = p, where p is maximized
    # The trick: p is the actual likelihood, we don't care about bernoulli!
    C <- 1000000
    for(d in 1:ndays){
        # Seasonal Affect
        lambda_spline[d] <- inprod(B[d,], lambeta[])
        # Presence spline
        pi_spline[d] <- inprod(B[d,], pibeta[])
        # Dispersion Spline
        alpha_spline[d] <- exp(inprod(B[d,], abeta[]))
        # Dispersion Spline
		for(r in 1:6){
			gamma_spline[r, d] <- inprod(B[d,], gambeta[r,1:14])
		}
    }
    
    for(i in 1:idays){
        # Ones trick
        ones[i] ~ dbern(p[i])
        p[i] <- Lik[i]/C
        Lik[i] <- exp(logLik[i])
        logLik[i] <- (1-z[i])*log(1-w[i]) + 
        #    z[i]*(log(w[i]) + logTPois[i])
            z[i]*(log(w[i]) + LogTruncNB[i])
        
        # For Hurdle
        logit(w[i]) <- pi_spline[DayYearN[i]]
        
        # Truncated Poisson Distribution
        #logTPois[i] <- log(w[i]) + Ni[i]*log(mu[i]) -
        #    mu[i] - log(1 - exp(-mu[i])) - loggam(Ni[i] + 1)
            
        # https://www.cambridge.org/core/services/aop-cambridge-core/content/view/C1FFB63EC7657CD424CE6726F32BC4FE/9781316459515c7_p184-214_CBO.pdf/glms_part_iii_zeroinflated_and_hurdle_models.pdf
        # https://ncss-wpengine.netdna-ssl.com/wp-content/themes/ncss/pdf/Procedures/NCSS/Zero-Inflated_Negative_Binomial_Regression.pdf
        LogTruncNB[i] <-  (1/alpha_spline[DayYearN[i]])*log(u[i])+
            Ni[i]*log(1 - u[i]) + 
            loggam(Ni[i] + 1/alpha_spline[DayYearN[i]]) -
            loggam(1/alpha_spline[DayYearN[i]]) - 
            loggam(Ni[i] + 1) -
            log(1 - (1 + alpha_spline[DayYearN[i]]*mu[i])^(-1/alpha_spline[DayYearN[i]]))
         u[i] <- 1/(1 + alpha_spline[DayYearN[i]]*mu[i])
         
        #https://data.princeton.edu/wws509/notes/countmoments
        ETrunc1[i] <- 1 + alpha_spline[DayYearN[i]]*mu[i]
        ETrunc2[i] <- -1/alpha_spline[DayYearN[i]]
        ETrunc3[i] <- ETrunc1[i]^ETrunc2[i]
        ETruncNB[i] <- mu[i]/(1 - ETrunc3[i])
        #ExTruncNB[i] <- mu[i]/((1-(1 + 
        #    alpha_spline[DayYearN[i]]*mu[i]))^(-1/alpha_spline[DayYearN[i]]))
        ENi[i] <- w[i]*ETruncNB[i]
        
        # mu is the mean of the Poisson
        # lambda[r] is region-specific intercept
        log(mu[i]) <- lambda_spline[DayYearN[i]] + lambda[FireCentre2N[i]] + 
            b[DayYearN[i]] #+ inprod(ncovars, betan)
    }
    
    
    for(j in 1:jfires){
        Xijr_star[j] ~ dlnorm(lp2[j], taux[FireCentre2[j]])
        lp1[j] <- inprod(xcovars[j,], betax)
        lp2[j] <- mux[FireCentre2[j]] + 
                lp1[j] + 
                gamma_spline[FireCentre2[j], DayYear[j]] * b[DayYear[j]]
        
        # px is probabilty of wider regime (digit preference rather than rounding)
        px[j] ~ dbern(roundprob[FireCentre2[j]])
        ind[j, 1] <- dint[j, 2*px[j] + 1]
        ind[j, 2] <- dint[j, 2*px[j] + 2]
            
        Xijr1[j] ~ dinterval(Xijr_star[j], ind[j,])
        
        EXi[j] <- exp(lp2[j] + 1/(2*taux[FireCentre2[j]]))
        XLik[j] <- log(plnorm(ind[j, 2], lp2[j], taux[FireCentre2[j]]) -
            plnorm(ind[j, 1], lp2[j], taux[FireCentre2[j]]))
    }
    
    # Priors
    taub <- pow(sigma_b, -2)
    sigma_b ~ dgamma(8, 2)
    
    # Splines
    lambeta ~ dmnorm(zero_vec, lambeta_posdef)
    pibeta ~ dmnorm(zero_vec, lambeta_posdef)
    abeta ~ dmnorm(zero_vec, lambeta_posdef)
    
    # AR(1) structure for random effects
    phi[1] ~ dnorm(0, 1/5)
    b[1] ~ dnorm(0, taub)
    for(i in 2:ndays){
        b[i] ~ dnorm(phi[1]*b[i-1], taub)
    }

    
    for(p in 1:nxcovar){
		betax[p] ~ dnorm(0, 1/4)
    }
    
    #for(p in 1:nNcovar){
    #    betan[p] ~ dnorm(0, 1/4)
    #}
    
    for(r in 1:6){
        #gammai[r] ~ dnorm(hypgamma, 1/5)
		gambeta[r, 1:14] ~ dmnorm(zero_vec, lambeta_posdef)
        mux[r] ~ dnorm(hypmux, 1/5)
        lambda[r] ~ dnorm(hyplambda, 1/5)
        taux[r] <- pow(sigma_x[r], -2)
        sigma_x[r] ~ dgamma(8, 2)
        roundprob[r] ~ dbeta(2,2)
        
    }
	hypgamma ~ dnorm(0, 1/5)
	hypmux ~ dnorm(0, 1/5)
	hyplambda ~ dnorm(0, 1/5)
}
"}

# Province-wide splines
# Prov-wide re
# AR1 re's
# Prov-wide beta
# Yes hyperparameters
mix_jags_AR1pbeta_regpi <- {"
model{
    # For the ones trick
    # C chosen so that Lik/c <= 1
    # P(ones == 1) = p, where p is maximized
    # The trick: p is the actual likelihood, we don't care about bernoulli!
    C <- 1000000
    for(d in 1:ndays){
        # Seasonal Affect
        lambda_spline[d] <- inprod(B[d,], lambeta[])
        # Presence spline
		for(r in 1:6){
            pi_spline[r,d] <- inprod(B[d,], pibeta[r,1:J])
		}
        # Dispersion Spline
        alpha_spline[d] <- exp(inprod(B[d,], abeta[]))
    }
    
    for(i in 1:idays){
        # Ones trick
        ones[i] ~ dbern(p[i])
        p[i] <- Lik[i]/C
        Lik[i] <- exp(logLik[i])
        logLik[i] <- (1-z[i])*log(1-w[i]) + 
        #    z[i]*(log(w[i]) + logTPois[i])
            z[i]*(log(w[i]) + LogTruncNB[i])
        
        # For Hurdle
        logit(w[i]) <- pi_spline[FireCentre2N[i], DayYearN[i]]
        
        # Truncated Poisson Distribution
        #logTPois[i] <- log(w[i]) + Ni[i]*log(mu[i]) -
        #    mu[i] - log(1 - exp(-mu[i])) - loggam(Ni[i] + 1)
            
        # https://www.cambridge.org/core/services/aop-cambridge-core/content/view/C1FFB63EC7657CD424CE6726F32BC4FE/9781316459515c7_p184-214_CBO.pdf/glms_part_iii_zeroinflated_and_hurdle_models.pdf
        # https://ncss-wpengine.netdna-ssl.com/wp-content/themes/ncss/pdf/Procedures/NCSS/Zero-Inflated_Negative_Binomial_Regression.pdf
        LogTruncNB[i] <-  (1/alpha_spline[DayYearN[i]])*log(u[i])+
            Ni[i]*log(1 - u[i]) + 
            loggam(Ni[i] + 1/alpha_spline[DayYearN[i]]) -
            loggam(1/alpha_spline[DayYearN[i]]) - 
            loggam(Ni[i] + 1) -
            log(1 - (1 + alpha_spline[DayYearN[i]]*mu[i])^(-1/alpha_spline[DayYearN[i]]))
         u[i] <- 1/(1 + alpha_spline[DayYearN[i]]*mu[i])
         
        #https://data.princeton.edu/wws509/notes/countmoments
        ETrunc1[i] <- 1 + alpha_spline[DayYearN[i]]*mu[i]
        ETrunc2[i] <- -1/alpha_spline[DayYearN[i]]
        ETrunc3[i] <- ETrunc1[i]^ETrunc2[i]
        ETruncNB[i] <- mu[i]/(1 - ETrunc3[i])
        #ExTruncNB[i] <- mu[i]/((1-(1 + 
        #    alpha_spline[DayYearN[i]]*mu[i]))^(-1/alpha_spline[DayYearN[i]]))
        ENi[i] <- w[i]*ETruncNB[i]
        
        # mu is the mean of the Poisson
        # lambda[r] is region-specific intercept
        log(mu[i]) <- lambda_spline[DayYearN[i]] + lambda[FireCentre2N[i]] + 
            b[DayYearN[i]] #+ inprod(ncovars, betan)
    }
    
    
    for(j in 1:jfires){
        Xijr_star[j] ~ dlnorm(lp2[j], taux[FireCentre2[j]])
        lp1[j] <- inprod(xcovars[j,], betax)
        lp2[j] <- mux[FireCentre2[j]] + 
                lp1[j] + 
                gammai[FireCentre2[j]] * b[DayYear[j]]
        
        # px is probabilty of wider regime (digit preference rather than rounding)
        px[j] ~ dbern(roundprob[FireCentre2[j]])
        ind[j, 1] <- dint[j, 2*px[j] + 1]
        ind[j, 2] <- dint[j, 2*px[j] + 2]
            
        Xijr1[j] ~ dinterval(Xijr_star[j], ind[j,])
        
        EXi[j] <- exp(lp2[j] + 1/(2*taux[FireCentre2[j]]))
        XLik[j] <- log(plnorm(ind[j, 2], lp2[j], taux[FireCentre2[j]]) -
            plnorm(ind[j, 1], lp2[j], taux[FireCentre2[j]]))
    }
    
    # Priors
    taub <- pow(sigma_b, -2)
    sigma_b ~ dgamma(8, 2)
    
    # Splines
    lambeta ~ dmnorm(zero_vec, lambeta_posdef)
    for(r in 1:6){
        pibeta[r,1:J] ~ dmnorm(zero_vec, lambeta_posdef)
	}
    abeta ~ dmnorm(zero_vec, lambeta_posdef)
    
    # AR(1) structure for random effects
    phi[1] ~ dnorm(0, 1/5)
    b[1] ~ dnorm(0, taub)
    for(i in 2:ndays){
        b[i] ~ dnorm(phi[1]*b[i-1], taub)
    }

    
    for(p in 1:nxcovar){
		betax[p] ~ dnorm(0, 1/4)
    }
    
    #for(p in 1:nNcovar){
    #    betan[p] ~ dnorm(0, 1/4)
    #}
    
    for(r in 1:6){
        gammai[r] ~ dnorm(hypgamma, 1/5)
        mux[r] ~ dnorm(hypmux, 1/5)
        lambda[r] ~ dnorm(hyplambda, 1/5)
        taux[r] <- pow(sigma_x[r], -2)
        sigma_x[r] ~ dgamma(8, 2)
        roundprob[r] ~ dbeta(2,2)
        
    }
	hypgamma ~ dnorm(0, 1/5)
	hypmux ~ dnorm(0, 1/5)
	hyplambda ~ dnorm(0, 1/5)
}
"}





# Sensitivity analyses

mix_jags_AR1pbeta_gamsel_priorfn <- function(
	phi_prior = "dnorm(0, 1/5)",
	gammai_prior = "dnorm(0, 1/5)",
	mux_prior = "dnorm(0, 1/5)",
	lambdar_prior = "dnorm(0, 1/5)",
	sigmax_prior = "dgamma(8, 2)",
	sigmab_prior = "dgamma(8, 2)",
	betax_prior = "dnorm(0, 1/4)"){paste0("
model{
    # For the ones trick
    # C chosen so that Lik/c <= 1
    # P(ones == 1) = p, where p is maximized
    # The trick: p is the actual likelihood, we don't care about bernoulli!
    C <- 1000000
    for(d in 1:ndays){
        # Seasonal Affect
        lambda_spline[d] <- inprod(B[d,], lambeta[])
        # Presence spline
        pi_spline[d] <- inprod(B[d,], pibeta[])
        # Dispersion Spline
        alpha_spline[d] <- exp(inprod(B[d,], abeta[]))
    }
    
    for(i in 1:idays){
        # Ones trick
        ones[i] ~ dbern(p[i])
        p[i] <- Lik[i]/C
        Lik[i] <- exp(logLik[i])
        logLik[i] <- (1-z[i])*log(1-w[i]) + 
        #    z[i]*(log(w[i]) + logTPois[i])
            z[i]*(log(w[i]) + LogTruncNB[i])
        
        # For Hurdle
        logit(w[i]) <- pi_spline[DayYearN[i]]
        
        # Truncated Poisson Distribution
        #logTPois[i] <- log(w[i]) + Ni[i]*log(mu[i]) -
        #    mu[i] - log(1 - exp(-mu[i])) - loggam(Ni[i] + 1)
            
        # https://www.cambridge.org/core/services/aop-cambridge-core/content/view/C1FFB63EC7657CD424CE6726F32BC4FE/9781316459515c7_p184-214_CBO.pdf/glms_part_iii_zeroinflated_and_hurdle_models.pdf
        # https://ncss-wpengine.netdna-ssl.com/wp-content/themes/ncss/pdf/Procedures/NCSS/Zero-Inflated_Negative_Binomial_Regression.pdf
        LogTruncNB[i] <-  (1/alpha_spline[DayYearN[i]])*log(u[i])+
            Ni[i]*log(1 - u[i]) + 
            loggam(Ni[i] + 1/alpha_spline[DayYearN[i]]) -
            loggam(1/alpha_spline[DayYearN[i]]) - 
            loggam(Ni[i] + 1) -
            log(1 - (1 + alpha_spline[DayYearN[i]]*mu[i])^(-1/alpha_spline[DayYearN[i]]))
         u[i] <- 1/(1 + alpha_spline[DayYearN[i]]*mu[i])
         
        #https://data.princeton.edu/wws509/notes/countmoments
        ETrunc1[i] <- 1 + alpha_spline[DayYearN[i]]*mu[i]
        ETrunc2[i] <- -1/alpha_spline[DayYearN[i]]
        ETrunc3[i] <- ETrunc1[i]^ETrunc2[i]
        ETruncNB[i] <- mu[i]/(1 - ETrunc3[i])
        #ExTruncNB[i] <- mu[i]/((1-(1 + 
        #    alpha_spline[DayYearN[i]]*mu[i]))^(-1/alpha_spline[DayYearN[i]]))
        ENi[i] <- w[i]*ETruncNB[i]
        
        # mu is the mean of the Poisson
        # lambda[r] is region-specific intercept
        log(mu[i]) <- lambda_spline[DayYearN[i]] + lambda[FireCentre2N[i]] + 
            b[DayYearN[i]] #+ inprod(ncovars, betan)
    }
    
    
    for(j in 1:jfires){
        Xijr_star[j] ~ dlnorm(lp2[j], taux[FireCentre2[j]])
        lp1[j] <- inprod(xcovars[j,], betax)
        lp2[j] <- mux[FireCentre2[j]] + 
                lp1[j] + 
                gammai[FireCentre2[j]] * b[DayYear[j]]
        
        # px is probabilty of wider regime (digit preference rather than rounding)
        px[j] ~ dbern(roundprob[FireCentre2[j]])
        ind[j, 1] <- dint[j, 2*px[j] + 1]
        ind[j, 2] <- dint[j, 2*px[j] + 2]
            
        Xijr1[j] ~ dinterval(Xijr_star[j], ind[j,])
        
        EXi[j] <- exp(lp2[j] + 1/(2*taux[FireCentre2[j]]))
        XLik[j] <- log(plnorm(ind[j, 2], lp2[j], taux[FireCentre2[j]]) -
            plnorm(ind[j, 1], lp2[j], taux[FireCentre2[j]]))
    }
    
    # Priors
    taub <- pow(sigma_b, -2)
    sigma_b ~ ", sigmab_prior,"
    
    # Splines
    lambeta ~ dmnorm(zero_vec, lambeta_posdef)
    pibeta ~ dmnorm(zero_vec, lambeta_posdef)
    abeta ~ dmnorm(zero_vec, lambeta_posdef)
    
    # AR(1) structure for random effects
    phi[1] ~ ", phi_prior, "
    b[1] ~ dnorm(0, taub)
    for(i in 2:ndays){
        b[i] ~ dnorm(phi[1]*b[i-1], taub)
    }

    
    for(p in 1:nxcovar){
		betax[p] ~ ", betax_prior, "
    }
    
    #for(p in 1:nNcovar){
    #    betan[p] ~ dnorm(0, 1/4)
    #}
    
	gind ~ dbeta(1,1)
    for(r in 1:6){
		gzero[r] ~ dbern(gind)
        gammai_1[r] ~ ", gammai_prior, "
		gammai[r] <- gzero[r]*gammai_1[r]
        mux[r] ~ ", mux_prior, "
        lambda[r] ~ ", lambdar_prior, "
        taux[r] <- pow(sigma_x[r], -2)
        sigma_x[r] ~ ", sigmax_prior,"
        roundprob[r] ~ dbeta(2,2)
    }
}
")}



mix_jags_AR1pbetahyp_gamsel_priorfn <- function(
	phi_prior = "dnorm(0, 1/5)",
	hypgam_prior = "dnorm(0, 1/5)",
	hypmux_prior = "dnorm(0, 1/5)",
	hyplam_prior = "dnorm(0, 1/5)",
	sigmax_prior = "dgamma(8, 2)",
	sigmab_prior = "dgamma(8, 2)",
	betax_prior = "dnorm(0, 1/4)"){paste0("
model{
    # For the ones trick
    # C chosen so that Lik/c <= 1
    # P(ones == 1) = p, where p is maximized
    # The trick: p is the actual likelihood, we don't care about bernoulli!
    C <- 1000000
    for(d in 1:ndays){
        # Seasonal Affect
        lambda_spline[d] <- inprod(B[d,], lambeta[])
        # Presence spline
        pi_spline[d] <- inprod(B[d,], pibeta[])
        # Dispersion Spline
        alpha_spline[d] <- exp(inprod(B[d,], abeta[]))
    }
    
    for(i in 1:idays){
        # Ones trick
        ones[i] ~ dbern(p[i])
        p[i] <- Lik[i]/C
        Lik[i] <- exp(logLik[i])
        logLik[i] <- (1-z[i])*log(1-w[i]) + 
        #    z[i]*(log(w[i]) + logTPois[i])
            z[i]*(log(w[i]) + LogTruncNB[i])
        
        # For Hurdle
        logit(w[i]) <- pi_spline[DayYearN[i]]
        
        # Truncated Poisson Distribution
        #logTPois[i] <- log(w[i]) + Ni[i]*log(mu[i]) -
        #    mu[i] - log(1 - exp(-mu[i])) - loggam(Ni[i] + 1)
            
        # https://www.cambridge.org/core/services/aop-cambridge-core/content/view/C1FFB63EC7657CD424CE6726F32BC4FE/9781316459515c7_p184-214_CBO.pdf/glms_part_iii_zeroinflated_and_hurdle_models.pdf
        # https://ncss-wpengine.netdna-ssl.com/wp-content/themes/ncss/pdf/Procedures/NCSS/Zero-Inflated_Negative_Binomial_Regression.pdf
        LogTruncNB[i] <-  (1/alpha_spline[DayYearN[i]])*log(u[i])+
            Ni[i]*log(1 - u[i]) + 
            loggam(Ni[i] + 1/alpha_spline[DayYearN[i]]) -
            loggam(1/alpha_spline[DayYearN[i]]) - 
            loggam(Ni[i] + 1) -
            log(1 - (1 + alpha_spline[DayYearN[i]]*mu[i])^(-1/alpha_spline[DayYearN[i]]))
         u[i] <- 1/(1 + alpha_spline[DayYearN[i]]*mu[i])
         
        #https://data.princeton.edu/wws509/notes/countmoments
        ETrunc1[i] <- 1 + alpha_spline[DayYearN[i]]*mu[i]
        ETrunc2[i] <- -1/alpha_spline[DayYearN[i]]
        ETrunc3[i] <- ETrunc1[i]^ETrunc2[i]
        ETruncNB[i] <- mu[i]/(1 - ETrunc3[i])
        #ExTruncNB[i] <- mu[i]/((1-(1 + 
        #    alpha_spline[DayYearN[i]]*mu[i]))^(-1/alpha_spline[DayYearN[i]]))
        ENi[i] <- w[i]*ETruncNB[i]
        
        # mu is the mean of the Poisson
        # lambda[r] is region-specific intercept
        log(mu[i]) <- lambda_spline[DayYearN[i]] + lambda[FireCentre2N[i]] + 
            b[DayYearN[i]] #+ inprod(ncovars, betan)
    }
    
    
    for(j in 1:jfires){
        Xijr_star[j] ~ dlnorm(lp2[j], taux[FireCentre2[j]])
        lp1[j] <- inprod(xcovars[j,], betax)
        lp2[j] <- mux[FireCentre2[j]] + 
                lp1[j] + 
                gammai[FireCentre2[j]] * b[DayYear[j]]
        
        # px is probabilty of wider regime (digit preference rather than rounding)
        px[j] ~ dbern(roundprob[FireCentre2[j]])
        ind[j, 1] <- dint[j, 2*px[j] + 1]
        ind[j, 2] <- dint[j, 2*px[j] + 2]
            
        Xijr1[j] ~ dinterval(Xijr_star[j], ind[j,])
        
        EXi[j] <- exp(lp2[j] + 1/(2*taux[FireCentre2[j]]))
        XLik[j] <- log(plnorm(ind[j, 2], lp2[j], taux[FireCentre2[j]]) -
            plnorm(ind[j, 1], lp2[j], taux[FireCentre2[j]]))
    }
    
    # Priors
    taub <- pow(sigma_b, -2)
    sigma_b ~ ", sigmab_prior,"
    
    # Splines
    lambeta ~ dmnorm(zero_vec, lambeta_posdef)
    pibeta ~ dmnorm(zero_vec, lambeta_posdef)
    abeta ~ dmnorm(zero_vec, lambeta_posdef)
    
    # AR(1) structure for random effects
    phi[1] ~ ", phi_prior, "
    b[1] ~ dnorm(0, taub)
    for(i in 2:ndays){
        b[i] ~ dnorm(phi[1]*b[i-1], taub)
    }

    
    for(p in 1:nxcovar){
		betax[p] ~ ", betax_prior, "
    }
    
    #for(p in 1:nNcovar){
    #    betan[p] ~ dnorm(0, 1/4)
    #}
    
	gind ~ dbeta(1,1)
    for(r in 1:6){
		gzero[r] ~ dbern(gind)
        gammai_1[r] ~ dnorm(hypgam, 1/5)
		gammai[r] <- gzero[r]*gammai_1[r]
        mux[r] ~ dnorm(hypmux, 1/5)
        lambda[r] ~ dnorm(hyplam, 1/5)
        taux[r] <- pow(sigma_x[r], -2)
        sigma_x[r] ~ ", sigmax_prior,"
        roundprob[r] ~ dbeta(2,2)
    }
	
	hypgam ~ ", hypgam_prior,"
	hypmux ~ ", hypmux_prior,"
	hyplam ~ ", hyplam_prior,"
}
")}




# Remember the Poisson Models?


mix_jags_pois_base <- {"
model{
    # For the ones trick
    # C chosen so that Lik/c <= 1
    # P(ones == 1) = p, where p is maximized
    # The trick: p is the actual likelihood, we don't care about bernoulli!
    C <- 1000000
    for(d in 1:ndays){
        b[d] ~ dnorm(0, taub)
        
        # Seasonal Affect
        lambda_spline[d] <- inprod(B[d,], lambeta[])
        # Presence spline
        pi_spline[d] <- inprod(B[d,], pibeta[])
        # Dispersion Spline
        #alpha_spline[d] <- exp(inprod(B[d,], abeta[]))
    }
    
    for(i in 1:idays){
        # Ones trick
        ones[i] ~ dbern(p[i])
        p[i] <- Lik[i]/C
        Lik[i] <- exp(logLik[i])
        logLik[i] <- (1-z[i])*log(1-w[i]) + 
            z[i]*(log(w[i]) + logTPois[i])
        #    z[i]*(log(w[i]) + LogTruncNB[i])
        
        # For Hurdle
        logit(w[i]) <- pi_spline[DayYearN[i]]
        
        # Truncated Poisson Distribution
        logTPois[i] <- log(w[i]) + Ni[i]*log(mu[i]) -
            mu[i] - log(1 - exp(-mu[i])) - loggam(Ni[i] + 1)
            
        # https://www.cambridge.org/core/services/aop-cambridge-core/content/view/C1FFB63EC7657CD424CE6726F32BC4FE/9781316459515c7_p184-214_CBO.pdf/glms_part_iii_zeroinflated_and_hurdle_models.pdf
        # https://ncss-wpengine.netdna-ssl.com/wp-content/themes/ncss/pdf/Procedures/NCSS/Zero-Inflated_Negative_Binomial_Regression.pdf
        #LogTruncNB[i] <-  (1/alpha_spline[DayYearN[i]])*log(u[i])+
        #    Ni[i]*log(1 - u[i]) + 
        #    loggam(Ni[i] + 1/alpha_spline[DayYearN[i]]) -
        #    loggam(1/alpha_spline[DayYearN[i]]) - 
        #    loggam(Ni[i] + 1) -
        #    log(1 - (1 + alpha_spline[DayYearN[i]]*mu[i])^(-1/alpha_spline[DayYearN[i]]))
        # u[i] <- 1/(1 + alpha_spline[DayYearN[i]]*mu[i])
         
        #https://data.princeton.edu/wws509/notes/countmoments
        #ETrunc1[i] <- 1 + alpha_spline[DayYearN[i]]*mu[i]
        #ETrunc2[i] <- -1/alpha_spline[DayYearN[i]]
        #ETrunc3[i] <- ETrunc1[i]^ETrunc2[i]
        #ETruncNB[i] <- mu[i]/(1 - ETrunc3[i])
        #ExTruncNB[i] <- mu[i]/((1-(1 + 
        #    alpha_spline[DayYearN[i]]*mu[i]))^(-1/alpha_spline[DayYearN[i]]))
        #ENi[i] <- w[i]*ETruncNB[i]

		ENi[i] <- mu[i]/(1 + exp(-mu[i]))
        
        # mu is the mean of the Poisson
        # lambda[r] is region-specific intercept
        log(mu[i]) <- lambda_spline[DayYearN[i]] + lambda[FireCentre2N[i]] + 
            b[DayYearN[i]] #+ inprod(ncovars, betan)
    }
    
    
    for(j in 1:jfires){
        Xijr_star[j] ~ dlnorm(lp2[j], taux[FireCentre2[j]])
        lp1[j] <- inprod(xcovars[j,], betax[, FireCentre2[j]])
        lp2[j] <- mux[FireCentre2[j]] + 
                lp1[j] + 
                gammai[FireCentre2[j]] * b[DayYear[j]]
        
        # px is probabilty of wider regime (digit preference rather than rounding)
        px[j] ~ dbern(roundprob[FireCentre2[j]])
        ind[j, 1] <- dint[j, 2*px[j] + 1]
        ind[j, 2] <- dint[j, 2*px[j] + 2]
            
        Xijr1[j] ~ dinterval(Xijr_star[j], ind[j,])
        
        EXi[j] <- exp(lp2[j] + 1/(2*taux[FireCentre2[j]]))
        XLik[j] <- log(plnorm(ind[j, 2], lp2[j], taux[FireCentre2[j]]) -
            plnorm(ind[j, 1], lp2[j], taux[FireCentre2[j]]))
    }
    
    # Priors
    taub <- pow(sigma_b, -2)
    sigma_b ~ dgamma(8, 2)
    
    # Splines
    lambeta ~ dmnorm(zero_vec, lambeta_posdef)
    pibeta ~ dmnorm(zero_vec, lambeta_posdef)
    #abeta ~ dmnorm(zero_vec, lambeta_posdef)
    

    
    for(p in 1:nxcovar){
        for(r in 1:6){
            betax[p, r] ~ dnorm(0, 1/4)
        }
    }
    
    #for(p in 1:nNcovar){
    #    betan[p] ~ dnorm(0, 1/4)
    #}
    
    for(r in 1:6){
        gammai[r] ~ dnorm(0, 1/5)
        mux[r] ~ dnorm(0, 1/5)
        lambda[r] ~ dnorm(0, 1/5)
        taux[r] <- pow(sigma_x[r], -2)
        sigma_x[r] ~ dgamma(8, 2)
        roundprob[r] ~ dbeta(2,2)
        
    }
}
"}



# Covars for Ni

mix_jags_AR1pbeta_NiCov <- {"
model{
    # For the ones trick
    # C chosen so that Lik/c <= 1
    # P(ones == 1) = p, where p is maximized
    # The trick: p is the actual likelihood, we don't care about bernoulli!
    C <- 1000000
    for(d in 1:ndays){
        # Seasonal Affect
        lambda_spline[d] <- inprod(B[d,], lambeta[])
        # Presence spline
        pi_spline[d] <- inprod(B[d,], pibeta[])
        # Dispersion Spline
        alpha_spline[d] <- exp(inprod(B[d,], abeta[]))
    }
    
    for(i in 1:idays){
        # Ones trick
        ones[i] ~ dbern(p[i])
        p[i] <- Lik[i]/C
        Lik[i] <- exp(logLik[i])
        logLik[i] <- (1-z[i])*log(1-w[i]) + 
        #    z[i]*(log(w[i]) + logTPois[i])
            z[i]*(log(w[i]) + LogTruncNB[i])
        
        # For Hurdle
        logit(w[i]) <- pi_spline[DayYearN[i]]
        
        # Truncated Poisson Distribution
        #logTPois[i] <- log(w[i]) + Ni[i]*log(mu[i]) -
        #    mu[i] - log(1 - exp(-mu[i])) - loggam(Ni[i] + 1)
            
        # https://www.cambridge.org/core/services/aop-cambridge-core/content/view/C1FFB63EC7657CD424CE6726F32BC4FE/9781316459515c7_p184-214_CBO.pdf/glms_part_iii_zeroinflated_and_hurdle_models.pdf
        # https://ncss-wpengine.netdna-ssl.com/wp-content/themes/ncss/pdf/Procedures/NCSS/Zero-Inflated_Negative_Binomial_Regression.pdf
        LogTruncNB[i] <-  (1/alpha_spline[DayYearN[i]])*log(u[i])+
            Ni[i]*log(1 - u[i]) + 
            loggam(Ni[i] + 1/alpha_spline[DayYearN[i]]) -
            loggam(1/alpha_spline[DayYearN[i]]) - 
            loggam(Ni[i] + 1) -
            log(1 - (1 + alpha_spline[DayYearN[i]]*mu[i])^(-1/alpha_spline[DayYearN[i]]))
         u[i] <- 1/(1 + alpha_spline[DayYearN[i]]*mu[i])
         
        #https://data.princeton.edu/wws509/notes/countmoments
        ETrunc1[i] <- 1 + alpha_spline[DayYearN[i]]*mu[i]
        ETrunc2[i] <- -1/alpha_spline[DayYearN[i]]
        ETrunc3[i] <- ETrunc1[i]^ETrunc2[i]
        ETruncNB[i] <- mu[i]/(1 - ETrunc3[i])
        #ExTruncNB[i] <- mu[i]/((1-(1 + 
        #    alpha_spline[DayYearN[i]]*mu[i]))^(-1/alpha_spline[DayYearN[i]]))
        ENi[i] <- w[i]*ETruncNB[i]
        
        # mu is the mean of the Poisson
        # lambda[r] is region-specific intercept
        log(mu[i]) <- lambda_spline[DayYearN[i]] + lambda[FireCentre2N[i]] + 
            b[DayYearN[i]] + inprod(ncovars[i,], betan)
    }
    
    
    for(j in 1:jfires){
        Xijr_star[j] ~ dlnorm(lp2[j], taux[FireCentre2[j]])
        lp1[j] <- inprod(xcovars[j,], betax)
        lp2[j] <- mux[FireCentre2[j]] + 
                lp1[j] + 
                gammai[FireCentre2[j]] * b[DayYear[j]]
        
        # px is probabilty of wider regime (digit preference rather than rounding)
        px[j] ~ dbern(roundprob[FireCentre2[j]])
        ind[j, 1] <- dint[j, 2*px[j] + 1]
        ind[j, 2] <- dint[j, 2*px[j] + 2]
            
        Xijr1[j] ~ dinterval(Xijr_star[j], ind[j,])
        
        EXi[j] <- exp(lp2[j] + 1/(2*taux[FireCentre2[j]]))
        XLik[j] <- log(plnorm(ind[j, 2], lp2[j], taux[FireCentre2[j]]) -
            plnorm(ind[j, 1], lp2[j], taux[FireCentre2[j]]))
    }
    
    # Priors
    taub <- pow(sigma_b, -2)
    sigma_b ~ dgamma(8, 2)
    
    # Splines
    lambeta ~ dmnorm(zero_vec, lambeta_posdef)
    pibeta ~ dmnorm(zero_vec, lambeta_posdef)
    abeta ~ dmnorm(zero_vec, lambeta_posdef)
    
    # AR(1) structure for random effects
    phi[1] ~ dnorm(0, 1/5)
    b[1] ~ dnorm(0, taub)
    for(i in 2:ndays){
        b[i] ~ dnorm(phi[1]*b[i-1], taub)
    }

    
    for(p in 1:nxcovar){
        betax[p] ~ dnorm(0, 1/4)
    }
    
    for(p in 1:nNcovar){
        betan[p] ~ dnorm(0, 1/4)
    }
    
    for(r in 1:6){
        gammai[r] ~ dnorm(hypgamma, 1/5)
        mux[r] ~ dnorm(hypmux, 1/5)
        lambda[r] ~ dnorm(hyplambda, 1/5)
        taux[r] <- pow(sigma_x[r], -2)
        sigma_x[r] ~ dgamma(8, 2)
        roundprob[r] ~ dbeta(2,2)
        
    }
	hypgamma ~ dnorm(0, 1/5)
	hypmux ~ dnorm(0, 1/5)
	hyplambda ~ dnorm(0, 1/5)
}
"}


# Possible final?

wide_jags_AR1pbeta_NiCov_regpi <- {"
model{
    # For the ones trick
    # C chosen so that Lik/c <= 1
    # P(ones == 1) = p, where p is maximized
    # The trick: p is the actual likelihood, we don't care about bernoulli!
    C <- 1000000
    for(d in 1:ndays){
        # Seasonal Affect
        lambda_spline[d] <- inprod(B[d,], lambeta[])
        # Presence spline
		for(r in 1:6){
            pi_spline[r,d] <- inprod(B[d,], pibeta[r,1:J])
		}
        # Dispersion Spline
        alpha_spline[d] <- exp(inprod(B[d,], abeta[]))
    }
    
    for(i in 1:idays){
        # Ones trick
        ones[i] ~ dbern(p[i])
        p[i] <- Lik[i]/C
        Lik[i] <- exp(logLik[i])
        logLik[i] <- (1-z[i])*log(1-w[i]) + 
        #    z[i]*(log(w[i]) + logTPois[i])
            z[i]*(log(w[i]) + LogTruncNB[i])
        
        # For Hurdle
        logit(w[i]) <- pi_spline[FireCentre2N[i], DayYearN[i]]
        
        # Truncated Poisson Distribution
        #logTPois[i] <- log(w[i]) + Ni[i]*log(mu[i]) -
        #    mu[i] - log(1 - exp(-mu[i])) - loggam(Ni[i] + 1)
            
        # https://www.cambridge.org/core/services/aop-cambridge-core/content/view/C1FFB63EC7657CD424CE6726F32BC4FE/9781316459515c7_p184-214_CBO.pdf/glms_part_iii_zeroinflated_and_hurdle_models.pdf
        # https://ncss-wpengine.netdna-ssl.com/wp-content/themes/ncss/pdf/Procedures/NCSS/Zero-Inflated_Negative_Binomial_Regression.pdf
        LogTruncNB[i] <-  (1/alpha_spline[DayYearN[i]])*log(u[i])+
            Ni[i]*log(1 - u[i]) + 
            loggam(Ni[i] + 1/alpha_spline[DayYearN[i]]) -
            loggam(1/alpha_spline[DayYearN[i]]) - 
            loggam(Ni[i] + 1) -
            log(1 - (1 + alpha_spline[DayYearN[i]]*mu[i])^(-1/alpha_spline[DayYearN[i]]))
         u[i] <- 1/(1 + alpha_spline[DayYearN[i]]*mu[i])
         
        #https://data.princeton.edu/wws509/notes/countmoments
        ETrunc1[i] <- 1 + alpha_spline[DayYearN[i]]*mu[i]
        ETrunc2[i] <- -1/alpha_spline[DayYearN[i]]
        ETrunc3[i] <- ETrunc1[i]^ETrunc2[i]
        ETruncNB[i] <- mu[i]/(1 - ETrunc3[i])
        #ExTruncNB[i] <- mu[i]/((1-(1 + 
        #    alpha_spline[DayYearN[i]]*mu[i]))^(-1/alpha_spline[DayYearN[i]]))
        ENi[i] <- w[i]*ETruncNB[i]
        
        # mu is the mean of the Poisson
        # lambda[r] is region-specific intercept
        log(mu[i]) <- lambda_spline[DayYearN[i]] + lambda[FireCentre2N[i]] + 
            b[DayYearN[i]] + inprod(ncovars[i,], betan)
    }
    
    
    for(j in 1:jfires){
        Xijr_star[j] ~ dlnorm(lp2[j], taux[FireCentre2[j]])
        lp1[j] <- inprod(xcovars[j,], betax)
        lp2[j] <- mux[FireCentre2[j]] + 
                lp1[j] + 
                gammai[FireCentre2[j]] * b[DayYear[j]]
        
        Xijr1[j] ~ dinterval(Xijr_star[j], dint[j,])
        
        EXi[j] <- exp(lp2[j] + 1/(2*taux[FireCentre2[j]]))
        XLik[j] <- log(plnorm(dint[j, 2], lp2[j], taux[FireCentre2[j]]) -
            plnorm(dint[j, 1], lp2[j], taux[FireCentre2[j]]))
    }
    
    # Priors
    taub <- pow(sigma_b, -2)
    sigma_b ~ dgamma(8, 2)
    
    # Splines
    lambeta ~ dmnorm(zero_vec, lambeta_posdef)
    for(r in 1:6){
        pibeta[r,1:J] ~ dmnorm(zero_vec, lambeta_posdef)
	}
    abeta ~ dmnorm(zero_vec, lambeta_posdef)
    
    # AR(1) structure for random effects
    phi[1] ~ dnorm(0, 1/5)
    b[1] ~ dnorm(0, taub)
    for(i in 2:ndays){
        b[i] ~ dnorm(phi[1]*b[i-1], taub)
    }

    
    for(p in 1:nxcovar){
        betax[p] ~ dnorm(0, 1/4)
    }
    
    for(p in 1:nNcovar){
        betan[p] ~ dnorm(0, 1/4)
    }
    
    for(r in 1:6){
        gammai[r] ~ dnorm(hypgamma, 1/5)
        mux[r] ~ dnorm(hypmux, 1/5)
        lambda[r] ~ dnorm(hyplambda, 1/5)
        taux[r] <- pow(sigma_x[r], -2)
        sigma_x[r] ~ dgamma(8, 2)
        #roundprob[r] ~ dbeta(2,2)
        
    }
	hypgamma ~ dnorm(0, 1/5)
	hypmux ~ dnorm(0, 1/5)
	hyplambda ~ dnorm(0, 1/5)
}
"}

# Possible final?

wide_jags_AR1pbeta_NiCov_regpi_sep <- {"
model{
    # For the ones trick
    # C chosen so that Lik/c <= 1
    # P(ones == 1) = p, where p is maximized
    # The trick: p is the actual likelihood, we don't care about bernoulli!
    C <- 1000000
    for(d in 1:ndays){
        # Seasonal Affect
        lambda_spline[d] <- inprod(B[d,], lambeta[])
        # Presence spline
		for(r in 1:6){
            pi_spline[r,d] <- inprod(B[d,], pibeta[r,1:J])
		}
        # Dispersion Spline
        alpha_spline[d] <- exp(inprod(B[d,], abeta[]))
    }
    
    for(i in 1:idays){
        # Ones trick
        ones[i] ~ dbern(p[i])
        p[i] <- Lik[i]/C
        Lik[i] <- exp(logLik[i])
        logLik[i] <- (1-z[i])*log(1-w[i]) + 
        #    z[i]*(log(w[i]) + logTPois[i])
            z[i]*(log(w[i]) + LogTruncNB[i])
        
        # For Hurdle
        logit(w[i]) <- pi_spline[FireCentre2N[i], DayYearN[i]]
        
        # Truncated Poisson Distribution
        #logTPois[i] <- log(w[i]) + Ni[i]*log(mu[i]) -
        #    mu[i] - log(1 - exp(-mu[i])) - loggam(Ni[i] + 1)
            
        # https://www.cambridge.org/core/services/aop-cambridge-core/content/view/C1FFB63EC7657CD424CE6726F32BC4FE/9781316459515c7_p184-214_CBO.pdf/glms_part_iii_zeroinflated_and_hurdle_models.pdf
        # https://ncss-wpengine.netdna-ssl.com/wp-content/themes/ncss/pdf/Procedures/NCSS/Zero-Inflated_Negative_Binomial_Regression.pdf
        LogTruncNB[i] <-  (1/alpha_spline[DayYearN[i]])*log(u[i])+
            Ni[i]*log(1 - u[i]) + 
            loggam(Ni[i] + 1/alpha_spline[DayYearN[i]]) -
            loggam(1/alpha_spline[DayYearN[i]]) - 
            loggam(Ni[i] + 1) -
            log(1 - (1 + alpha_spline[DayYearN[i]]*mu[i])^(-1/alpha_spline[DayYearN[i]]))
         u[i] <- 1/(1 + alpha_spline[DayYearN[i]]*mu[i])
         
        #https://data.princeton.edu/wws509/notes/countmoments
        ETrunc1[i] <- 1 + alpha_spline[DayYearN[i]]*mu[i]
        ETrunc2[i] <- -1/alpha_spline[DayYearN[i]]
        ETrunc3[i] <- ETrunc1[i]^ETrunc2[i]
        ETruncNB[i] <- mu[i]/(1 - ETrunc3[i])
        #ExTruncNB[i] <- mu[i]/((1-(1 + 
        #    alpha_spline[DayYearN[i]]*mu[i]))^(-1/alpha_spline[DayYearN[i]]))
        ENi[i] <- w[i]*ETruncNB[i]
        
        # mu is the mean of the Poisson
        # lambda[r] is region-specific intercept
        log(mu[i]) <- lambda_spline[DayYearN[i]] + lambda[FireCentre2N[i]] + 
            b[DayYearN[i]] + inprod(ncovars[i,], betan)
    }
    
    
    for(j in 1:jfires){
        Xijr_star[j] ~ dlnorm(lp2[j], taux[FireCentre2[j]])
        lp1[j] <- inprod(xcovars[j,], betax)
        lp2[j] <- mux[FireCentre2[j]] + 
                lp1[j] + 
                gammai[FireCentre2[j]] * b[DayYear[j]]
        
        Xijr1[j] ~ dinterval(Xijr_star[j], dint[j,])
        
        EXi[j] <- exp(lp2[j] + 1/(2*taux[FireCentre2[j]]))
        XLik[j] <- log(plnorm(dint[j, 2], lp2[j], taux[FireCentre2[j]]) -
            plnorm(dint[j, 1], lp2[j], taux[FireCentre2[j]]))
    }
    
    # Priors
    taub <- pow(sigma_b, -2)
    sigma_b ~ dgamma(8, 2)
    
    # Splines
    lambeta ~ dmnorm(zero_vec, lambeta_posdef)
    for(r in 1:6){
        pibeta[r,1:J] ~ dmnorm(zero_vec, lambeta_posdef)
	}
    abeta ~ dmnorm(zero_vec, lambeta_posdef)
    
    # AR(1) structure for random effects
    phi[1] ~ dnorm(0, 1/5)
    b[1] ~ dnorm(0, taub)
    for(i in 2:ndays){
        b[i] ~ dnorm(phi[1]*b[i-1], taub)
    }

    
    for(p in 1:nxcovar){
        betax[p] ~ dnorm(0, 1/4)
    }
    
    for(p in 1:nNcovar){
        betan[p] ~ dnorm(0, 1/4)
    }
    
    for(r in 1:6){
        gammai[r] <- 0 #~ dnorm(hypgamma, 1/5)
        mux[r] ~ dnorm(hypmux, 1/5)
        lambda[r] ~ dnorm(hyplambda, 1/5)
        taux[r] <- pow(sigma_x[r], -2)
        sigma_x[r] ~ dgamma(8, 2)
        #roundprob[r] ~ dbeta(2,2)
        
    }
	#hypgamma ~ dnorm(0, 1/5)
	hypmux ~ dnorm(0, 1/5)
	hyplambda ~ dnorm(0, 1/5)
}
"}


# Possible final?

wide_jags_AR1pbeta_NiCov <- {"
model{
    # For the ones trick
    # C chosen so that Lik/c <= 1
    # P(ones == 1) = p, where p is maximized
    # The trick: p is the actual likelihood, we don't care about bernoulli!
    C <- 1000000
    for(d in 1:ndays){
        # Seasonal Affect
        lambda_spline[d] <- inprod(B[d,], lambeta[])
        # Presence spline
        pi_spline[d] <- inprod(B[d,], pibeta[])
        # Dispersion Spline
        alpha_spline[d] <- exp(inprod(B[d,], abeta[]))
    }
    
    for(i in 1:idays){
        # Ones trick
        ones[i] ~ dbern(p[i])
        p[i] <- Lik[i]/C
        Lik[i] <- exp(logLik[i])
        logLik[i] <- (1-z[i])*log(1-w[i]) + 
        #    z[i]*(log(w[i]) + logTPois[i])
            z[i]*(log(w[i]) + LogTruncNB[i])
        
        # For Hurdle
        logit(w[i]) <- pi_spline[DayYearN[i]]
        
        # Truncated Poisson Distribution
        #logTPois[i] <- log(w[i]) + Ni[i]*log(mu[i]) -
        #    mu[i] - log(1 - exp(-mu[i])) - loggam(Ni[i] + 1)
            
        # https://www.cambridge.org/core/services/aop-cambridge-core/content/view/C1FFB63EC7657CD424CE6726F32BC4FE/9781316459515c7_p184-214_CBO.pdf/glms_part_iii_zeroinflated_and_hurdle_models.pdf
        # https://ncss-wpengine.netdna-ssl.com/wp-content/themes/ncss/pdf/Procedures/NCSS/Zero-Inflated_Negative_Binomial_Regression.pdf
        LogTruncNB[i] <-  (1/alpha_spline[DayYearN[i]])*log(u[i])+
            Ni[i]*log(1 - u[i]) + 
            loggam(Ni[i] + 1/alpha_spline[DayYearN[i]]) -
            loggam(1/alpha_spline[DayYearN[i]]) - 
            loggam(Ni[i] + 1) -
            log(1 - (1 + alpha_spline[DayYearN[i]]*mu[i])^(-1/alpha_spline[DayYearN[i]]))
         u[i] <- 1/(1 + alpha_spline[DayYearN[i]]*mu[i])
         
        #https://data.princeton.edu/wws509/notes/countmoments
        ETrunc1[i] <- 1 + alpha_spline[DayYearN[i]]*mu[i]
        ETrunc2[i] <- -1/alpha_spline[DayYearN[i]]
        ETrunc3[i] <- ETrunc1[i]^ETrunc2[i]
        ETruncNB[i] <- mu[i]/(1 - ETrunc3[i])
        #ExTruncNB[i] <- mu[i]/((1-(1 + 
        #    alpha_spline[DayYearN[i]]*mu[i]))^(-1/alpha_spline[DayYearN[i]]))
        ENi[i] <- w[i]*ETruncNB[i]
        
        # mu is the mean of the Poisson
        # lambda[r] is region-specific intercept
        log(mu[i]) <- lambda_spline[DayYearN[i]] + lambda[FireCentre2N[i]] + 
            b[DayYearN[i]] + inprod(ncovars[i,], betan)
    }
    
    
    for(j in 1:jfires){
        Xijr_star[j] ~ dlnorm(lp2[j], taux[FireCentre2[j]])
        lp1[j] <- inprod(xcovars[j,], betax)
        lp2[j] <- mux[FireCentre2[j]] + 
                lp1[j] + 
                gammai[FireCentre2[j]] * b[DayYear[j]]
        
        Xijr1[j] ~ dinterval(Xijr_star[j], dint[j,])
        
        EXi[j] <- exp(lp2[j] + 1/(2*taux[FireCentre2[j]]))
        XLik[j] <- log(plnorm(dint[j, 2], lp2[j], taux[FireCentre2[j]]) -
            plnorm(dint[j, 1], lp2[j], taux[FireCentre2[j]]))
    }
    
    # Priors
    taub <- pow(sigma_b, -2)
    sigma_b ~ dgamma(8, 2)
    
    # Splines
    lambeta ~ dmnorm(zero_vec, lambeta_posdef)
    pibeta ~ dmnorm(zero_vec, lambeta_posdef)
    abeta ~ dmnorm(zero_vec, lambeta_posdef)
    
    # AR(1) structure for random effects
    phi[1] ~ dnorm(0, 1/5)
    b[1] ~ dnorm(0, taub)
    for(i in 2:ndays){
        b[i] ~ dnorm(phi[1]*b[i-1], taub)
    }

    
    for(p in 1:nxcovar){
        betax[p] ~ dnorm(0, 1/4)
    }
    
    for(p in 1:nNcovar){
        betan[p] ~ dnorm(0, 1/4)
    }
    
    for(r in 1:6){
        gammai[r] ~ dnorm(hypgamma, 1/5)
        mux[r] ~ dnorm(hypmux, 1/5)
        lambda[r] ~ dnorm(hyplambda, 1/5)
        taux[r] <- pow(sigma_x[r], -2)
        sigma_x[r] ~ dgamma(8, 2)
        #roundprob[r] ~ dbeta(2,2)
        
    }
	hypgamma ~ dnorm(0, 1/5)
	hypmux ~ dnorm(0, 1/5)
	hyplambda ~ dnorm(0, 1/5)
}
"}



wide_jags_AR1pbeta_NiCov_gspline <- {"
model{
    # For the ones trick
    # C chosen so that Lik/c <= 1
    # P(ones == 1) = p, where p is maximized
    # The trick: p is the actual likelihood, we don't care about bernoulli!
    C <- 1000000
    for(d in 1:ndays){
        # Seasonal Affect
        lambda_spline[d] <- inprod(B[d,], lambeta[])
        # Presence spline
        pi_spline[d] <- inprod(B[d,], pibeta[])
        # Dispersion Spline
        alpha_spline[d] <- exp(inprod(B[d,], abeta[]))
        # Dispersion Spline
		for(r in 1:6){
			gamma_spline[r, d] <- inprod(B[d,], gambeta[r,1:14])
		}
    }
    
    for(i in 1:idays){
        # Ones trick
        ones[i] ~ dbern(p[i])
        p[i] <- Lik[i]/C
        Lik[i] <- exp(logLik[i])
        logLik[i] <- (1-z[i])*log(1-w[i]) + 
        #    z[i]*(log(w[i]) + logTPois[i])
            z[i]*(log(w[i]) + LogTruncNB[i])
        
        # For Hurdle
        logit(w[i]) <- pi_spline[DayYearN[i]]
        
        # Truncated Poisson Distribution
        #logTPois[i] <- log(w[i]) + Ni[i]*log(mu[i]) -
        #    mu[i] - log(1 - exp(-mu[i])) - loggam(Ni[i] + 1)
            
        # https://www.cambridge.org/core/services/aop-cambridge-core/content/view/C1FFB63EC7657CD424CE6726F32BC4FE/9781316459515c7_p184-214_CBO.pdf/glms_part_iii_zeroinflated_and_hurdle_models.pdf
        # https://ncss-wpengine.netdna-ssl.com/wp-content/themes/ncss/pdf/Procedures/NCSS/Zero-Inflated_Negative_Binomial_Regression.pdf
        LogTruncNB[i] <-  (1/alpha_spline[DayYearN[i]])*log(u[i])+
            Ni[i]*log(1 - u[i]) + 
            loggam(Ni[i] + 1/alpha_spline[DayYearN[i]]) -
            loggam(1/alpha_spline[DayYearN[i]]) - 
            loggam(Ni[i] + 1) -
            log(1 - (1 + alpha_spline[DayYearN[i]]*mu[i])^(-1/alpha_spline[DayYearN[i]]))
         u[i] <- 1/(1 + alpha_spline[DayYearN[i]]*mu[i])
         
        #https://data.princeton.edu/wws509/notes/countmoments
        ETrunc1[i] <- 1 + alpha_spline[DayYearN[i]]*mu[i]
        ETrunc2[i] <- -1/alpha_spline[DayYearN[i]]
        ETrunc3[i] <- ETrunc1[i]^ETrunc2[i]
        ETruncNB[i] <- mu[i]/(1 - ETrunc3[i])
        #ExTruncNB[i] <- mu[i]/((1-(1 + 
        #    alpha_spline[DayYearN[i]]*mu[i]))^(-1/alpha_spline[DayYearN[i]]))
        ENi[i] <- w[i]*ETruncNB[i]
        
        # mu is the mean of the Poisson
        # lambda[r] is region-specific intercept
        log(mu[i]) <- lambda_spline[DayYearN[i]] + lambda[FireCentre2N[i]] + 
            b[DayYearN[i]] + inprod(ncovars[i,], betan)
    }
    
    
    for(j in 1:jfires){
        Xijr_star[j] ~ dlnorm(lp2[j], taux[FireCentre2[j]])
        lp1[j] <- inprod(xcovars[j,], betax)
        lp2[j] <- mux[FireCentre2[j]] + 
                lp1[j] + 
                gamma_spline[FireCentre2[j], DayYear[j]] * b[DayYear[j]]
        
        Xijr1[j] ~ dinterval(Xijr_star[j], dint[j,])
        
        EXi[j] <- exp(lp2[j] + 1/(2*taux[FireCentre2[j]]))
        XLik[j] <- log(plnorm(dint[j, 2], lp2[j], taux[FireCentre2[j]]) -
            plnorm(dint[j, 1], lp2[j], taux[FireCentre2[j]]))
    }
    
    # Priors
    taub <- pow(sigma_b, -2)
    sigma_b ~ dgamma(8, 2)
    
    # Splines
    lambeta ~ dmnorm(zero_vec, lambeta_posdef)
    pibeta ~ dmnorm(zero_vec, lambeta_posdef)
    abeta ~ dmnorm(zero_vec, lambeta_posdef)
    
    # AR(1) structure for random effects
    phi[1] ~ dnorm(0, 1/5)
    b[1] ~ dnorm(0, taub)
    for(i in 2:ndays){
        b[i] ~ dnorm(phi[1]*b[i-1], taub)
    }

    
    for(p in 1:nxcovar){
        betax[p] ~ dnorm(0, 1/4)
    }
    
    for(p in 1:nNcovar){
        betan[p] ~ dnorm(0, 1/4)
    }
    
    for(r in 1:6){
        #gammai[r] ~ dnorm(hypgamma, 1/5)
		gambeta[r, 1:14] ~ dmnorm(zero_vec, lambeta_posdef)
        mux[r] ~ dnorm(hypmux, 1/5)
        lambda[r] ~ dnorm(hyplambda, 1/5)
        taux[r] <- pow(sigma_x[r], -2)
        sigma_x[r] ~ dgamma(8, 2)
        #roundprob[r] ~ dbeta(2,2)
        
    }
	hypgamma ~ dnorm(0, 1/5)
	hypmux ~ dnorm(0, 1/5)
	hyplambda ~ dnorm(0, 1/5)
}
"}



wide_jags_AR1rbeta_NiCov <- {"
model{
    # For the ones trick
    # C chosen so that Lik/c <= 1
    # P(ones == 1) = p, where p is maximized
    # The trick: p is the actual likelihood, we don't care about bernoulli!
    C <- 1000000
    for(d in 1:ndays){
        # Seasonal Affect
        lambda_spline[d] <- inprod(B[d,], lambeta[])
        # Presence spline
        pi_spline[d] <- inprod(B[d,], pibeta[])
        # Dispersion Spline
        alpha_spline[d] <- exp(inprod(B[d,], abeta[]))
    }
    
    for(i in 1:idays){
        # Ones trick
        ones[i] ~ dbern(p[i])
        p[i] <- Lik[i]/C
        Lik[i] <- exp(logLik[i])
        logLik[i] <- (1-z[i])*log(1-w[i]) + 
        #    z[i]*(log(w[i]) + logTPois[i])
            z[i]*(log(w[i]) + LogTruncNB[i])
        
        # For Hurdle
        logit(w[i]) <- pi_spline[DayYearN[i]]
        
        # Truncated Poisson Distribution
        #logTPois[i] <- log(w[i]) + Ni[i]*log(mu[i]) -
        #    mu[i] - log(1 - exp(-mu[i])) - loggam(Ni[i] + 1)
            
        # https://www.cambridge.org/core/services/aop-cambridge-core/content/view/C1FFB63EC7657CD424CE6726F32BC4FE/9781316459515c7_p184-214_CBO.pdf/glms_part_iii_zeroinflated_and_hurdle_models.pdf
        # https://ncss-wpengine.netdna-ssl.com/wp-content/themes/ncss/pdf/Procedures/NCSS/Zero-Inflated_Negative_Binomial_Regression.pdf
        LogTruncNB[i] <-  (1/alpha_spline[DayYearN[i]])*log(u[i])+
            Ni[i]*log(1 - u[i]) + 
            loggam(Ni[i] + 1/alpha_spline[DayYearN[i]]) -
            loggam(1/alpha_spline[DayYearN[i]]) - 
            loggam(Ni[i] + 1) -
            log(1 - (1 + alpha_spline[DayYearN[i]]*mu[i])^(-1/alpha_spline[DayYearN[i]]))
         u[i] <- 1/(1 + alpha_spline[DayYearN[i]]*mu[i])
         
        #https://data.princeton.edu/wws509/notes/countmoments
        ETrunc1[i] <- 1 + alpha_spline[DayYearN[i]]*mu[i]
        ETrunc2[i] <- -1/alpha_spline[DayYearN[i]]
        ETrunc3[i] <- ETrunc1[i]^ETrunc2[i]
        ETruncNB[i] <- mu[i]/(1 - ETrunc3[i])
        #ExTruncNB[i] <- mu[i]/((1-(1 + 
        #    alpha_spline[DayYearN[i]]*mu[i]))^(-1/alpha_spline[DayYearN[i]]))
        ENi[i] <- w[i]*ETruncNB[i]
        
        # mu is the mean of the Poisson
        # lambda[r] is region-specific intercept
        log(mu[i]) <- lambda_spline[DayYearN[i]] + lambda[FireCentre2N[i]] + 
            b[DayYearN[i]] + inprod(ncovars[i,], betan)
    }
    
    
    for(j in 1:jfires){
        Xijr_star[j] ~ dlnorm(lp2[j], taux[FireCentre2[j]])
        lp1[j] <- inprod(xcovars[j,], betax[, FireCentre2[j]])
        lp2[j] <- mux[FireCentre2[j]] + 
                lp1[j] + 
                gammai[FireCentre2[j]] * b[DayYear[j]]
        
        Xijr1[j] ~ dinterval(Xijr_star[j], dint[j,])
        
        EXi[j] <- exp(lp2[j] + 1/(2*taux[FireCentre2[j]]))
        XLik[j] <- log(plnorm(dint[j, 2], lp2[j], taux[FireCentre2[j]]) -
            plnorm(dint[j, 1], lp2[j], taux[FireCentre2[j]]))
    }
    
    # Priors
    taub <- pow(sigma_b, -2)
    sigma_b ~ dgamma(8, 2)
    
    # Splines
    lambeta ~ dmnorm(zero_vec, lambeta_posdef)
    pibeta ~ dmnorm(zero_vec, lambeta_posdef)
    abeta ~ dmnorm(zero_vec, lambeta_posdef)
    
    # AR(1) structure for random effects
    phi[1] ~ dnorm(0, 1/5)
    b[1] ~ dnorm(0, taub)
    for(i in 2:ndays){
        b[i] ~ dnorm(phi[1]*b[i-1], taub)
    }

    
    for(p in 1:nxcovar){
		for(r in 1:6){
			betax[r, p] ~ dnorm(0, 1/4)
		}
    }
    
    for(p in 1:nNcovar){
		betan[p] ~ dnorm(0, 1/4)
    }
    
    for(r in 1:6){
        gammai[r] ~ dnorm(hypgamma, 1/5)
        mux[r] ~ dnorm(hypmux, 1/5)
        lambda[r] ~ dnorm(hyplambda, 1/5)
        taux[r] <- pow(sigma_x[r], -2)
        sigma_x[r] ~ dgamma(8, 2)
        #roundprob[r] ~ dbeta(2,2)
        
    }
	hypgamma ~ dnorm(0, 1/5)
	hypmux ~ dnorm(0, 1/5)
	hyplambda ~ dnorm(0, 1/5)
}
"}




wide_jags_AR1pbeta_NiCov_rega <- {"
model{
    # For the ones trick
    # C chosen so that Lik/c <= 1
    # P(ones == 1) = p, where p is maximized
    # The trick: p is the actual likelihood, we don't care about bernoulli!
    C <- 1000000
    for(d in 1:ndays){
        # Seasonal Affect
        lambda_spline[d] <- inprod(B[d,], lambeta[])
        # Presence spline
        pi_spline[d] <- inprod(B[d,], pibeta[])
        # Dispersion Spline
		for(r in 1:6){
            alpha_spline[r,d] <- 1# inprod(B[d,], abeta[r,1:J])
		}
    }
    
    for(i in 1:idays){
        # Ones trick
        ones[i] ~ dbern(p[i])
        p[i] <- Lik[i]/C
        Lik[i] <- exp(logLik[i])
        logLik[i] <- (1-z[i])*log(1-w[i]) + 
        #    z[i]*(log(w[i]) + logTPois[i])
            z[i]*(log(w[i]) + LogTruncNB[i])
        
        # For Hurdle
        logit(w[i]) <- pi_spline[DayYearN[i]]
        
        # Truncated Poisson Distribution
        #logTPois[i] <- log(w[i]) + Ni[i]*log(mu[i]) -
        #    mu[i] - log(1 - exp(-mu[i])) - loggam(Ni[i] + 1)
            
        # https://www.cambridge.org/core/services/aop-cambridge-core/content/view/C1FFB63EC7657CD424CE6726F32BC4FE/9781316459515c7_p184-214_CBO.pdf/glms_part_iii_zeroinflated_and_hurdle_models.pdf
        # https://ncss-wpengine.netdna-ssl.com/wp-content/themes/ncss/pdf/Procedures/NCSS/Zero-Inflated_Negative_Binomial_Regression.pdf
        LogTruncNB[i] <-  (1/alpha_spline[FireCentre2N[i], DayYearN[i]])*log(u[i])+
            Ni[i]*log(1 - u[i]) + 
            loggam(Ni[i] + 1/alpha_spline[FireCentre2N[i], DayYearN[i]]) -
            loggam(1/alpha_spline[FireCentre2N[i], DayYearN[i]]) - 
            loggam(Ni[i] + 1) -
            log(1 - (1 + alpha_spline[FireCentre2N[i], DayYearN[i]]*mu[i])^(-1/alpha_spline[FireCentre2N[i], DayYearN[i]]))
         u[i] <- 1/(1 + alpha_spline[FireCentre2N[i], DayYearN[i]]*mu[i])
         
        #https://data.princeton.edu/wws509/notes/countmoments
        ETrunc1[i] <- 1 + alpha_spline[FireCentre2N[i], DayYearN[i]]*mu[i]
        ETrunc2[i] <- -1/alpha_spline[FireCentre2N[i], DayYearN[i]]
        ETrunc3[i] <- ETrunc1[i]^ETrunc2[i]
        ETruncNB[i] <- mu[i]/(1 - ETrunc3[i])
        #ExTruncNB[i] <- mu[i]/((1-(1 + 
        #    alpha_spline[FireCentre2N[i], DayYearN[i]]*mu[i]))^(-1/alpha_spline[FireCentre2N[i], DayYearN[i]]))
        ENi[i] <- w[i]*ETruncNB[i]
        
        # mu is the mean of the Poisson
        # lambda[r] is region-specific intercept
        log(mu[i]) <- lambda_spline[DayYearN[i]] + lambda[FireCentre2N[i]] + 
            b[DayYearN[i]] + inprod(ncovars[i,], betan)
    }
    
    
    for(j in 1:jfires){
        Xijr_star[j] ~ dlnorm(lp2[j], taux[FireCentre2[j]])
        lp1[j] <- inprod(xcovars[j,], betax)
        lp2[j] <- mux[FireCentre2[j]] + 
                lp1[j] + 
                gammai[FireCentre2[j]] * b[DayYear[j]]
        
        Xijr1[j] ~ dinterval(Xijr_star[j], dint[j,])
        
        EXi[j] <- exp(lp2[j] + 1/(2*taux[FireCentre2[j]]))
        XLik[j] <- log(plnorm(dint[j, 2], lp2[j], taux[FireCentre2[j]]) -
            plnorm(dint[j, 1], lp2[j], taux[FireCentre2[j]]))
    }
    
    # Priors
    taub <- pow(sigma_b, -2)
    sigma_b ~ dgamma(8, 2)
    
    # Splines
    lambeta ~ dmnorm(zero_vec, lambeta_posdef)
    pibeta ~ dmnorm(zero_vec, lambeta_posdef)
    for(r in 1:6){
        abeta[r,1:J] ~ dmnorm(zero_vec, lambeta_posdef)
	}
    
    # AR(1) structure for random effects
    phi[1] ~ dnorm(0, 1/5)
    b[1] ~ dnorm(0, taub)
    for(i in 2:ndays){
        b[i] ~ dnorm(phi[1]*b[i-1], taub)
    }

    
    for(p in 1:nxcovar){
        betax[p] ~ dnorm(0, 1/4)
    }
    
    for(p in 1:nNcovar){
        betan[p] ~ dnorm(0, 1/4)
    }
    
    for(r in 1:6){
        gammai[r] ~ dnorm(hypgamma, 1/5)
        mux[r] ~ dnorm(hypmux, 1/5)
        lambda[r] ~ dnorm(hyplambda, 1/5)
        taux[r] <- pow(sigma_x[r], -2)
        sigma_x[r] ~ dgamma(8, 2)
        #roundprob[r] ~ dbeta(2,2)
        
    }
	hypgamma ~ dnorm(0, 1/5)
	hypmux ~ dnorm(0, 1/5)
	hyplambda ~ dnorm(0, 1/5)
}
"}



wide_jags_regre_pbeta_NiCov <- {"
model{
    # For the ones trick
    # C chosen so that Lik/c <= 1
    # P(ones == 1) = p, where p is maximized
    # The trick: p is the actual likelihood, we don't care about bernoulli!
    C <- 1000000
    for(d in 1:ndays){
        # Seasonal Affect
        lambda_spline[d] <- inprod(B[d,], lambeta[])
        # Presence spline
        pi_spline[d] <- inprod(B[d,], pibeta[])
        # Dispersion Spline
        alpha_spline[d] <- exp(inprod(B[d,], abeta[]))
    }
    
    for(i in 1:idays){
        # Ones trick
        ones[i] ~ dbern(p[i])
        p[i] <- Lik[i]/C
        Lik[i] <- exp(logLik[i])
        logLik[i] <- (1-z[i])*log(1-w[i]) + 
        #    z[i]*(log(w[i]) + logTPois[i])
            z[i]*(log(w[i]) + LogTruncNB[i])
        
        # For Hurdle
        logit(w[i]) <- pi_spline[DayYearN[i]]
        
        # Truncated Poisson Distribution
        #logTPois[i] <- log(w[i]) + Ni[i]*log(mu[i]) -
        #    mu[i] - log(1 - exp(-mu[i])) - loggam(Ni[i] + 1)
            
        # https://www.cambridge.org/core/services/aop-cambridge-core/content/view/C1FFB63EC7657CD424CE6726F32BC4FE/9781316459515c7_p184-214_CBO.pdf/glms_part_iii_zeroinflated_and_hurdle_models.pdf
        # https://ncss-wpengine.netdna-ssl.com/wp-content/themes/ncss/pdf/Procedures/NCSS/Zero-Inflated_Negative_Binomial_Regression.pdf
        LogTruncNB[i] <-  (1/alpha_spline[DayYearN[i]])*log(u[i])+
            Ni[i]*log(1 - u[i]) + 
            loggam(Ni[i] + 1/alpha_spline[DayYearN[i]]) -
            loggam(1/alpha_spline[DayYearN[i]]) - 
            loggam(Ni[i] + 1) -
            log(1 - (1 + alpha_spline[DayYearN[i]]*mu[i])^(-1/alpha_spline[DayYearN[i]]))
         u[i] <- 1/(1 + alpha_spline[DayYearN[i]]*mu[i])
         
        #https://data.princeton.edu/wws509/notes/countmoments
        ETrunc1[i] <- 1 + alpha_spline[DayYearN[i]]*mu[i]
        ETrunc2[i] <- -1/alpha_spline[DayYearN[i]]
        ETrunc3[i] <- ETrunc1[i]^ETrunc2[i]
        ETruncNB[i] <- mu[i]/(1 - ETrunc3[i])
        #ExTruncNB[i] <- mu[i]/((1-(1 + 
        #    alpha_spline[DayYearN[i]]*mu[i]))^(-1/alpha_spline[DayYearN[i]]))
        ENi[i] <- w[i]*ETruncNB[i]
        
        # mu is the mean of the Poisson
        # lambda[r] is region-specific intercept
        log(mu[i]) <- lambda_spline[DayYearN[i]] + lambda[FireCentre2N[i]] + 
            b[FireCentre2N[i], DayYearN[i]] + inprod(ncovars[i,], betan)
    }
    
    
    for(j in 1:jfires){
        Xijr_star[j] ~ dlnorm(lp2[j], taux[FireCentre2[j]])
        lp1[j] <- inprod(xcovars[j,], betax)
        lp2[j] <- mux[FireCentre2[j]] + 
                lp1[j] + 
                gammai[FireCentre2[j]] * b[FireCentre2[j], DayYear[j]]
        
        Xijr1[j] ~ dinterval(Xijr_star[j], dint[j,])
        
        EXi[j] <- exp(lp2[j] + 1/(2*taux[FireCentre2[j]]))
        XLik[j] <- log(plnorm(dint[j, 2], lp2[j], taux[FireCentre2[j]]) -
            plnorm(dint[j, 1], lp2[j], taux[FireCentre2[j]]))
    }
    
    # Priors
    taub <- pow(sigma_b, -2)
    sigma_b ~ dgamma(8, 2)
    
    # Splines
    lambeta ~ dmnorm(zero_vec, lambeta_posdef)
    pibeta ~ dmnorm(zero_vec, lambeta_posdef)
    abeta ~ dmnorm(zero_vec, lambeta_posdef)
    
    # AR(1) structure for random effects
    for(i in 1:ndays){
        for(r in 1:6){
            b[r, i] ~ dnorm(0, taub)
        }
    }

    
    for(p in 1:nxcovar){
        betax[p] ~ dnorm(0, 1/4)
    }
    
    for(p in 1:nNcovar){
        betan[p] ~ dnorm(0, 1/4)
    }
    
    for(r in 1:6){
        gammai[r] ~ dnorm(hypgamma, 1/5)
        mux[r] ~ dnorm(hypmux, 1/5)
        lambda[r] ~ dnorm(hyplambda, 1/5)
        taux[r] <- pow(sigma_x[r], -2)
        sigma_x[r] ~ dgamma(8, 2)
        #roundprob[r] ~ dbeta(2,2)
        
    }
	hypgamma ~ dnorm(0, 1/5)
	hypmux ~ dnorm(0, 1/5)
	hyplambda ~ dnorm(0, 1/5)
}
"}




wide_jags_AR1pbeta_NiCov_bal <- {"
model{
    # For the ones trick
    # C chosen so that Lik/c <= 1
    # P(ones == 1) = p, where p is maximized
    # The trick: p is the actual likelihood, we don't care about bernoulli!
    C <- 1000000
    for(d in 1:ndays){
        # Seasonal Affect
        lambda_spline[d] <- inprod(B[d,], lambeta[])
        # Presence spline
        pi_spline[d] <- inprod(B[d,], pibeta[])
        # Dispersion Spline
        alpha_spline[d] <- exp(inprod(B[d,], abeta[]))
    }
    
    for(i in 1:idays){
        # Ones trick
        ones[i] ~ dbern(p[i])
        p[i] <- Lik[i]/C
        Lik[i] <- exp(logLik[i])
        logLik[i] <- (1-z[i])*log(1-w[i]) + 
        #    z[i]*(log(w[i]) + logTPois[i])
            z[i]*(log(w[i]) + LogTruncNB[i]*NiStar[i])
        
        # For Hurdle
        logit(w[i]) <- pi_spline[DayYearN[i]]
        
        # Truncated Poisson Distribution
        #logTPois[i] <- log(w[i]) + Ni[i]*log(mu[i]) -
        #    mu[i] - log(1 - exp(-mu[i])) - loggam(Ni[i] + 1)
            
        # https://www.cambridge.org/core/services/aop-cambridge-core/content/view/C1FFB63EC7657CD424CE6726F32BC4FE/9781316459515c7_p184-214_CBO.pdf/glms_part_iii_zeroinflated_and_hurdle_models.pdf
        # https://ncss-wpengine.netdna-ssl.com/wp-content/themes/ncss/pdf/Procedures/NCSS/Zero-Inflated_Negative_Binomial_Regression.pdf
        LogTruncNB[i] <-  (1/alpha_spline[DayYearN[i]])*log(u[i])+
            Ni[i]*log(1 - u[i]) + 
            loggam(Ni[i] + 1/alpha_spline[DayYearN[i]]) -
            loggam(1/alpha_spline[DayYearN[i]]) - 
            loggam(Ni[i] + 1) -
            log(1 - (1 + alpha_spline[DayYearN[i]]*mu[i])^(-1/alpha_spline[DayYearN[i]]))
         u[i] <- 1/(1 + alpha_spline[DayYearN[i]]*mu[i])
         
        #https://data.princeton.edu/wws509/notes/countmoments
        ETrunc1[i] <- 1 + alpha_spline[DayYearN[i]]*mu[i]
        ETrunc2[i] <- -1/alpha_spline[DayYearN[i]]
        ETrunc3[i] <- ETrunc1[i]^ETrunc2[i]
        ETruncNB[i] <- mu[i]/(1 - ETrunc3[i])
        #ExTruncNB[i] <- mu[i]/((1-(1 + 
        #    alpha_spline[DayYearN[i]]*mu[i]))^(-1/alpha_spline[DayYearN[i]]))
        ENi[i] <- w[i]*ETruncNB[i]
        
        # mu is the mean of the Poisson
        # lambda[r] is region-specific intercept
        log(mu[i]) <- lambda_spline[DayYearN[i]] + lambda[FireCentre2N[i]] + 
            b[DayYearN[i]] + inprod(ncovars[i,], betan)
    }
    
    
    for(j in 1:jfires){
        Xijr_star[j] ~ dlnorm(lp2[j], taux[FireCentre2[j]])
        lp1[j] <- inprod(xcovars[j,], betax)
        lp2[j] <- mux[FireCentre2[j]] + 
                lp1[j] + 
                gammai[FireCentre2[j]] * b[DayYear[j]]
        
        Xijr1[j] ~ dinterval(Xijr_star[j], dint[j,])
        
        EXi[j] <- exp(lp2[j] + 1/(2*taux[FireCentre2[j]]))
        XLik[j] <- log(plnorm(dint[j, 2], lp2[j], taux[FireCentre2[j]]) -
            plnorm(dint[j, 1], lp2[j], taux[FireCentre2[j]]))
    }
    
    # Priors
    taub <- pow(sigma_b, -2)
    sigma_b ~ dgamma(8, 2)
    
    # Splines
    lambeta ~ dmnorm(zero_vec, lambeta_posdef)
    pibeta ~ dmnorm(zero_vec, lambeta_posdef)
    abeta ~ dmnorm(zero_vec, lambeta_posdef)
    
    # AR(1) structure for random effects
    phi[1] ~ dnorm(0, 1/5)
    b[1] ~ dnorm(0, taub)
    for(i in 2:ndays){
        b[i] ~ dnorm(phi[1]*b[i-1], taub)
    }

    
    for(p in 1:nxcovar){
        betax[p] ~ dnorm(0, 1/4)
    }
    
    for(p in 1:nNcovar){
        betan[p] ~ dnorm(0, 1/4)
    }
    
    for(r in 1:6){
        gammai[r] ~ dnorm(hypgamma, 1/5)
        mux[r] ~ dnorm(hypmux, 1/5)
        lambda[r] ~ dnorm(hyplambda, 1/5)
        taux[r] <- pow(sigma_x[r], -2)
        sigma_x[r] ~ dgamma(8, 2)
        #roundprob[r] ~ dbeta(2,2)
        
    }
	hypgamma ~ dnorm(0, 1/5)
	hypmux ~ dnorm(0, 1/5)
	hyplambda ~ dnorm(0, 1/5)
}
"}




wide_jags_AR1pbeta_NiCov_regpi_varsel <- {"
model{
    # For the ones trick
    # C chosen so that Lik/c <= 1
    # P(ones == 1) = p, where p is maximized
    # The trick: p is the actual likelihood, we don't care about bernoulli!
    C <- 1000000
    for(d in 1:ndays){
        # Seasonal Affect
        lambda_spline[d] <- inprod(B[d,], lambeta[])
        # Presence spline
		for(r in 1:6){
            pi_spline[r,d] <- inprod(B[d,], pibeta[r,1:J])
		}
        # Dispersion Spline
        alpha_spline[d] <- exp(inprod(B[d,], abeta[]))
    }
    
    for(i in 1:idays){
        # Ones trick
        ones[i] ~ dbern(p[i])
        p[i] <- Lik[i]/C
        Lik[i] <- exp(logLik[i])
        logLik[i] <- (1-z[i])*log(1-w[i]) + 
        #    z[i]*(log(w[i]) + logTPois[i])
            z[i]*(log(w[i]) + LogTruncNB[i])
        
        # For Hurdle
        logit(w[i]) <- pi_spline[FireCentre2N[i], DayYearN[i]]
        
        # Truncated Poisson Distribution
        #logTPois[i] <- log(w[i]) + Ni[i]*log(mu[i]) -
        #    mu[i] - log(1 - exp(-mu[i])) - loggam(Ni[i] + 1)
            
        # https://www.cambridge.org/core/services/aop-cambridge-core/content/view/C1FFB63EC7657CD424CE6726F32BC4FE/9781316459515c7_p184-214_CBO.pdf/glms_part_iii_zeroinflated_and_hurdle_models.pdf
        # https://ncss-wpengine.netdna-ssl.com/wp-content/themes/ncss/pdf/Procedures/NCSS/Zero-Inflated_Negative_Binomial_Regression.pdf
        LogTruncNB[i] <-  (1/alpha_spline[DayYearN[i]])*log(u[i])+
            Ni[i]*log(1 - u[i]) + 
            loggam(Ni[i] + 1/alpha_spline[DayYearN[i]]) -
            loggam(1/alpha_spline[DayYearN[i]]) - 
            loggam(Ni[i] + 1) -
            log(1 - (1 + alpha_spline[DayYearN[i]]*mu[i])^(-1/alpha_spline[DayYearN[i]]))
         u[i] <- 1/(1 + alpha_spline[DayYearN[i]]*mu[i])
         
        #https://data.princeton.edu/wws509/notes/countmoments
        ETrunc1[i] <- 1 + alpha_spline[DayYearN[i]]*mu[i]
        ETrunc2[i] <- -1/alpha_spline[DayYearN[i]]
        ETrunc3[i] <- ETrunc1[i]^ETrunc2[i]
        ETruncNB[i] <- mu[i]/(1 - ETrunc3[i])
        #ExTruncNB[i] <- mu[i]/((1-(1 + 
        #    alpha_spline[DayYearN[i]]*mu[i]))^(-1/alpha_spline[DayYearN[i]]))
        ENi[i] <- w[i]*ETruncNB[i]
        
        # mu is the mean of the Poisson
        # lambda[r] is region-specific intercept
        log(mu[i]) <- lambda_spline[DayYearN[i]] + lambda[FireCentre2N[i]] + 
            b[DayYearN[i]] + inprod(ncovars[i,], betan)
    }
    
    
    for(j in 1:jfires){
        Xijr_star[j] ~ dlnorm(lp2[j], taux[FireCentre2[j]])
        lp1[j] <- inprod(xcovars[j,], betax[FireCentre2[j],])
        lp2[j] <- mux[FireCentre2[j]] + 
                lp1[j] + 
                gammai[FireCentre2[j]] * b[DayYear[j]]
        
        Xijr1[j] ~ dinterval(Xijr_star[j], dint[j,])
        
        EXi[j] <- exp(lp2[j] + 1/(2*taux[FireCentre2[j]]))
        XLik[j] <- log(plnorm(dint[j, 2], lp2[j], taux[FireCentre2[j]]) -
            plnorm(dint[j, 1], lp2[j], taux[FireCentre2[j]]))
    }
    
    # Priors
    taub <- pow(sigma_b, -2)
    sigma_b ~ dgamma(8, 2)
    
    # Splines
    lambeta ~ dmnorm(zero_vec, lambeta_posdef)
    for(r in 1:6){
        pibeta[r,1:J] ~ dmnorm(zero_vec, lambeta_posdef)
	}
    abeta ~ dmnorm(zero_vec, lambeta_posdef)
    
    # AR(1) structure for random effects
    phi[1] ~ dnorm(0, 1/5)
    b[1] ~ dnorm(0, taub)
    for(i in 2:ndays){
        b[i] ~ dnorm(phi[1]*b[i-1], taub)
    }

    
    for(p in 1:nNcovar){
        betan[p] ~ dnorm(0, 1/4)
    }
    
    pind ~ dbeta(2,8)
    for(p in 1:nxcovar){
		bind[p] ~ dbern(pind)
        for(r in 1:6){
            betaxT[p, r] ~ dnorm(0, 1/4)
            betax[p, r] <- bind[p]*betaxT[p, r]
        }
    }
    
    for(r in 1:6){
        gammai[r] ~ dnorm(hypgamma, 1/5)
        mux[r] ~ dnorm(hypmux, 1/5)
        lambda[r] ~ dnorm(hyplambda, 1/5)
        taux[r] <- pow(sigma_x[r], -2)
        sigma_x[r] ~ dgamma(8, 2)
        #roundprob[r] ~ dbeta(2,2)
        
    }
	hypgamma ~ dnorm(0, 1/5)
	hypmux ~ dnorm(0, 1/5)
	hyplambda ~ dnorm(0, 1/5)
}
"}



wide_jags_AR1pbeta_NiCov_regpi_gamsel <- {"
model{
    # For the ones trick
    # C chosen so that Lik/c <= 1
    # P(ones == 1) = p, where p is maximized
    # The trick: p is the actual likelihood, we don't care about bernoulli!
    C <- 1000000
    for(d in 1:ndays){
        # Seasonal Affect
        lambda_spline[d] <- inprod(B[d,], lambeta[])
        # Presence spline
		for(r in 1:6){
            pi_spline[r,d] <- inprod(B[d,], pibeta[r,1:J])
		}
        # Dispersion Spline
        alpha_spline[d] <- exp(inprod(B[d,], abeta[]))
    }
    
    for(i in 1:idays){
        # Ones trick
        ones[i] ~ dbern(p[i])
        p[i] <- Lik[i]/C
        Lik[i] <- exp(logLik[i])
        logLik[i] <- (1-z[i])*log(1-w[i]) + 
        #    z[i]*(log(w[i]) + logTPois[i])
            z[i]*(log(w[i]) + LogTruncNB[i])
        
        # For Hurdle
        logit(w[i]) <- pi_spline[FireCentre2N[i], DayYearN[i]]
        
        # Truncated Poisson Distribution
        #logTPois[i] <- log(w[i]) + Ni[i]*log(mu[i]) -
        #    mu[i] - log(1 - exp(-mu[i])) - loggam(Ni[i] + 1)
            
        # https://www.cambridge.org/core/services/aop-cambridge-core/content/view/C1FFB63EC7657CD424CE6726F32BC4FE/9781316459515c7_p184-214_CBO.pdf/glms_part_iii_zeroinflated_and_hurdle_models.pdf
        # https://ncss-wpengine.netdna-ssl.com/wp-content/themes/ncss/pdf/Procedures/NCSS/Zero-Inflated_Negative_Binomial_Regression.pdf
        LogTruncNB[i] <-  (1/alpha_spline[DayYearN[i]])*log(u[i])+
            Ni[i]*log(1 - u[i]) + 
            loggam(Ni[i] + 1/alpha_spline[DayYearN[i]]) -
            loggam(1/alpha_spline[DayYearN[i]]) - 
            loggam(Ni[i] + 1) -
            log(1 - (1 + alpha_spline[DayYearN[i]]*mu[i])^(-1/alpha_spline[DayYearN[i]]))
         u[i] <- 1/(1 + alpha_spline[DayYearN[i]]*mu[i])
         
        #https://data.princeton.edu/wws509/notes/countmoments
        ETrunc1[i] <- 1 + alpha_spline[DayYearN[i]]*mu[i]
        ETrunc2[i] <- -1/alpha_spline[DayYearN[i]]
        ETrunc3[i] <- ETrunc1[i]^ETrunc2[i]
        ETruncNB[i] <- mu[i]/(1 - ETrunc3[i])
        #ExTruncNB[i] <- mu[i]/((1-(1 + 
        #    alpha_spline[DayYearN[i]]*mu[i]))^(-1/alpha_spline[DayYearN[i]]))
        ENi[i] <- w[i]*ETruncNB[i]
        
        # mu is the mean of the Poisson
        # lambda[r] is region-specific intercept
        log(mu[i]) <- lambda_spline[DayYearN[i]] + lambda[FireCentre2N[i]] + 
            b[DayYearN[i]] + inprod(ncovars[i,], betan)
    }
    
    
    for(j in 1:jfires){
        Xijr_star[j] ~ dlnorm(lp2[j], taux[FireCentre2[j]])
        lp1[j] <- inprod(xcovars[j,], betax[FireCentre2[j],])
        lp2[j] <- mux[FireCentre2[j]] + 
                lp1[j] + 
                gammai[FireCentre2[j]] * b[DayYear[j]]
        
        Xijr1[j] ~ dinterval(Xijr_star[j], dint[j,])
        
        EXi[j] <- exp(lp2[j] + 1/(2*taux[FireCentre2[j]]))
        XLik[j] <- log(plnorm(dint[j, 2], lp2[j], taux[FireCentre2[j]]) -
            plnorm(dint[j, 1], lp2[j], taux[FireCentre2[j]]))
    }
    
    # Priors
    taub <- pow(sigma_b, -2)
    sigma_b ~ dgamma(8, 2)
    
    # Splines
    lambeta ~ dmnorm(zero_vec, lambeta_posdef)
    for(r in 1:6){
        pibeta[r,1:J] ~ dmnorm(zero_vec, lambeta_posdef)
	}
    abeta ~ dmnorm(zero_vec, lambeta_posdef)
    
    # AR(1) structure for random effects
    phi[1] ~ dnorm(0, 1/5)
    b[1] ~ dnorm(0, taub)
    for(i in 2:ndays){
        b[i] ~ dnorm(phi[1]*b[i-1], taub)
    }

    
    for(p in 1:nNcovar){
        betan[p] ~ dnorm(0, 1/4)
    }

    for(p in 1:nxcovar){
        betax[p] ~ dnorm(0, 1/4)
    }
    
	gind ~ dbeta(2,8)
    for(r in 1:6){
        gammaiT[r] ~ dnorm(hypgamma, 1/5)
		gammai[r] <- gind*gammaiT[r]
        mux[r] ~ dnorm(hypmux, 1/5)
        lambda[r] ~ dnorm(hyplambda, 1/5)
        taux[r] <- pow(sigma_x[r], -2)
        sigma_x[r] ~ dgamma(8, 2)
        #roundprob[r] ~ dbeta(2,2)
        
    }
	hypgamma ~ dnorm(0, 1/5)
	hypmux ~ dnorm(0, 1/5)
	hyplambda ~ dnorm(0, 1/5)
}
"}




wide_jags_AR1pbeta_NiCov_gamsel <- {"
model{
    # For the ones trick
    # C chosen so that Lik/c <= 1
    # P(ones == 1) = p, where p is maximized
    # The trick: p is the actual likelihood, we don't care about bernoulli!
    C <- 1000000
    for(d in 1:ndays){
        # Seasonal Affect
        lambda_spline[d] <- inprod(B[d,], lambeta[])
        # Presence spline
        pi_spline[d] <- inprod(B[d,], pibeta[])
        # Dispersion Spline
        alpha_spline[d] <- exp(inprod(B[d,], abeta[]))
    }
    
    for(i in 1:idays){
        # Ones trick
        ones[i] ~ dbern(p[i])
        p[i] <- Lik[i]/C
        Lik[i] <- exp(logLik[i])
        logLik[i] <- (1-z[i])*log(1-w[i]) + 
        #    z[i]*(log(w[i]) + logTPois[i])
            z[i]*(log(w[i]) + LogTruncNB[i])
        
        # For Hurdle
        logit(w[i]) <- pi_spline[DayYearN[i]]
        
        # Truncated Poisson Distribution
        #logTPois[i] <- log(w[i]) + Ni[i]*log(mu[i]) -
        #    mu[i] - log(1 - exp(-mu[i])) - loggam(Ni[i] + 1)
            
        # https://www.cambridge.org/core/services/aop-cambridge-core/content/view/C1FFB63EC7657CD424CE6726F32BC4FE/9781316459515c7_p184-214_CBO.pdf/glms_part_iii_zeroinflated_and_hurdle_models.pdf
        # https://ncss-wpengine.netdna-ssl.com/wp-content/themes/ncss/pdf/Procedures/NCSS/Zero-Inflated_Negative_Binomial_Regression.pdf
        LogTruncNB[i] <-  (1/alpha_spline[DayYearN[i]])*log(u[i])+
            Ni[i]*log(1 - u[i]) + 
            loggam(Ni[i] + 1/alpha_spline[DayYearN[i]]) -
            loggam(1/alpha_spline[DayYearN[i]]) - 
            loggam(Ni[i] + 1) -
            log(1 - (1 + alpha_spline[DayYearN[i]]*mu[i])^(-1/alpha_spline[DayYearN[i]]))
         u[i] <- 1/(1 + alpha_spline[DayYearN[i]]*mu[i])
         
        #https://data.princeton.edu/wws509/notes/countmoments
        ETrunc1[i] <- 1 + alpha_spline[DayYearN[i]]*mu[i]
        ETrunc2[i] <- -1/alpha_spline[DayYearN[i]]
        ETrunc3[i] <- ETrunc1[i]^ETrunc2[i]
        ETruncNB[i] <- mu[i]/(1 - ETrunc3[i])
        #ExTruncNB[i] <- mu[i]/((1-(1 + 
        #    alpha_spline[DayYearN[i]]*mu[i]))^(-1/alpha_spline[DayYearN[i]]))
        ENi[i] <- w[i]*ETruncNB[i]
        
        # mu is the mean of the Poisson
        # lambda[r] is region-specific intercept
        log(mu[i]) <- lambda_spline[DayYearN[i]] + lambda[FireCentre2N[i]] + 
            b[DayYearN[i]] + inprod(ncovars[i,], betan)
    }
    
    
    for(j in 1:jfires){
        Xijr_star[j] ~ dlnorm(lp2[j], taux[FireCentre2[j]])
        lp1[j] <- inprod(xcovars[j,], betax)
        lp2[j] <- mux[FireCentre2[j]] + 
                lp1[j] + 
                gammai[FireCentre2[j]] * b[DayYear[j]]
        
        Xijr1[j] ~ dinterval(Xijr_star[j], dint[j,])
        
        EXi[j] <- exp(lp2[j] + 1/(2*taux[FireCentre2[j]]))
        XLik[j] <- log(plnorm(dint[j, 2], lp2[j], taux[FireCentre2[j]]) -
            plnorm(dint[j, 1], lp2[j], taux[FireCentre2[j]]))
    }
    
    # Priors
    taub <- pow(sigma_b, -2)
    sigma_b ~ dgamma(8, 2)
    
    # Splines
    lambeta ~ dmnorm(zero_vec, lambeta_posdef)
    pibeta ~ dmnorm(zero_vec, lambeta_posdef)
    abeta ~ dmnorm(zero_vec, lambeta_posdef)
    
    # AR(1) structure for random effects
    phi[1] ~ dnorm(0, 1/5)
    b[1] ~ dnorm(0, taub)
    for(i in 2:ndays){
        b[i] ~ dnorm(phi[1]*b[i-1], taub)
    }

    
    for(p in 1:nxcovar){
        betax[p] ~ dnorm(0, 1/4)
    }
    
    for(p in 1:nNcovar){
        betan[p] ~ dnorm(0, 1/4)
    }
    
	gind ~ dbeta(2,8)
    for(r in 1:6){
        gammaiT[r] ~ dnorm(hypgamma, 1/5)
		gammai[r] <- gind*gammaiT[r]
        mux[r] ~ dnorm(hypmux, 1/5)
        lambda[r] ~ dnorm(hyplambda, 1/5)
        taux[r] <- pow(sigma_x[r], -2)
        sigma_x[r] ~ dgamma(8, 2)
        #roundprob[r] ~ dbeta(2,2)
        
    }
	hypgamma ~ dnorm(0, 1/5)
	hypmux ~ dnorm(0, 1/5)
	hyplambda ~ dnorm(0, 1/5)
}
"}





wide_jags_AR1rbeta_NiCov_varsel <- {"
model{
    # For the ones trick
    # C chosen so that Lik/c <= 1
    # P(ones == 1) = p, where p is maximized
    # The trick: p is the actual likelihood, we don't care about bernoulli!
    C <- 1000000
    for(d in 1:ndays){
        # Seasonal Affect
        lambda_spline[d] <- inprod(B[d,], lambeta[])
        # Presence spline
        pi_spline[d] <- inprod(B[d,], pibeta[])
        # Dispersion Spline
        alpha_spline[d] <- exp(inprod(B[d,], abeta[]))
    }
    
    for(i in 1:idays){
        # Ones trick
        ones[i] ~ dbern(p[i])
        p[i] <- Lik[i]/C
        Lik[i] <- exp(logLik[i])
        logLik[i] <- (1-z[i])*log(1-w[i]) + 
        #    z[i]*(log(w[i]) + logTPois[i])
            z[i]*(log(w[i]) + LogTruncNB[i])
        
        # For Hurdle
        logit(w[i]) <- pi_spline[DayYearN[i]]
        
        # Truncated Poisson Distribution
        #logTPois[i] <- log(w[i]) + Ni[i]*log(mu[i]) -
        #    mu[i] - log(1 - exp(-mu[i])) - loggam(Ni[i] + 1)
            
        # https://www.cambridge.org/core/services/aop-cambridge-core/content/view/C1FFB63EC7657CD424CE6726F32BC4FE/9781316459515c7_p184-214_CBO.pdf/glms_part_iii_zeroinflated_and_hurdle_models.pdf
        # https://ncss-wpengine.netdna-ssl.com/wp-content/themes/ncss/pdf/Procedures/NCSS/Zero-Inflated_Negative_Binomial_Regression.pdf
        LogTruncNB[i] <-  (1/alpha_spline[DayYearN[i]])*log(u[i])+
            Ni[i]*log(1 - u[i]) + 
            loggam(Ni[i] + 1/alpha_spline[DayYearN[i]]) -
            loggam(1/alpha_spline[DayYearN[i]]) - 
            loggam(Ni[i] + 1) -
            log(1 - (1 + alpha_spline[DayYearN[i]]*mu[i])^(-1/alpha_spline[DayYearN[i]]))
         u[i] <- 1/(1 + alpha_spline[DayYearN[i]]*mu[i])
         
        #https://data.princeton.edu/wws509/notes/countmoments
        ETrunc1[i] <- 1 + alpha_spline[DayYearN[i]]*mu[i]
        ETrunc2[i] <- -1/alpha_spline[DayYearN[i]]
        ETrunc3[i] <- ETrunc1[i]^ETrunc2[i]
        ETruncNB[i] <- mu[i]/(1 - ETrunc3[i])
        #ExTruncNB[i] <- mu[i]/((1-(1 + 
        #    alpha_spline[DayYearN[i]]*mu[i]))^(-1/alpha_spline[DayYearN[i]]))
        ENi[i] <- w[i]*ETruncNB[i]
        
        # mu is the mean of the Poisson
        # lambda[r] is region-specific intercept
        log(mu[i]) <- lambda_spline[DayYearN[i]] + lambda[FireCentre2N[i]] + 
            b[DayYearN[i]] + inprod(ncovars[i,], betan)
    }
    
    
    for(j in 1:jfires){
        Xijr_star[j] ~ dlnorm(lp2[j], taux[FireCentre2[j]])
        lp1[j] <- inprod(xcovars[j,], betax[, FireCentre2[j]])
        lp2[j] <- mux[FireCentre2[j]] + 
                lp1[j] + 
                gammai[FireCentre2[j]] * b[DayYear[j]]
        
        Xijr1[j] ~ dinterval(Xijr_star[j], dint[j,])
        
        EXi[j] <- exp(lp2[j] + 1/(2*taux[FireCentre2[j]]))
        XLik[j] <- log(plnorm(dint[j, 2], lp2[j], taux[FireCentre2[j]]) -
            plnorm(dint[j, 1], lp2[j], taux[FireCentre2[j]]))
    }
    
    # Priors
    taub <- pow(sigma_b, -2)
    sigma_b ~ dgamma(8, 2)
    
    # Splines
    lambeta ~ dmnorm(zero_vec, lambeta_posdef)
    pibeta ~ dmnorm(zero_vec, lambeta_posdef)
    abeta ~ dmnorm(zero_vec, lambeta_posdef)
    
    # AR(1) structure for random effects
    phi[1] ~ dnorm(0, 1/5)
    b[1] ~ dnorm(0, taub)
    for(i in 2:ndays){
        b[i] ~ dnorm(phi[1]*b[i-1], taub)
    }

    
	pind ~ dbeta(2,8)
    for(p in 1:nxcovar){
		for(r in 1:6){
			bind[r, p] ~ dbern(pind)
			betaxT[r, p] ~ dnorm(0, 1/4)
			betax[r, p] <- betaxT[r, p]*bind[r,p]
		}
    }
    
    for(p in 1:nNcovar){
		betan[p] ~ dnorm(0, 1/4)
    }
    
    for(r in 1:6){
        gammai[r] ~ dnorm(hypgamma, 1/5)
        mux[r] ~ dnorm(hypmux, 1/5)
        lambda[r] ~ dnorm(hyplambda, 1/5)
        taux[r] <- pow(sigma_x[r], -2)
        sigma_x[r] ~ dgamma(8, 2)
        #roundprob[r] ~ dbeta(2,2)
        
    }
	hypgamma ~ dnorm(0, 1/5)
	hypmux ~ dnorm(0, 1/5)
	hyplambda ~ dnorm(0, 1/5)
}
"}



#### P-Splines instead of B-Splines

wide_jags_AR1rbeta_NiCov2 <- {"
model{
    # For the ones trick
    # C chosen so that Lik/c <= 1
    # P(ones == 1) = p, where p is maximized
    # The trick: p is the actual likelihood, we don't care about bernoulli!
    C <- 1000000
    for(d in 1:ndays){
        # Seasonal Affect
        lambda_spline[d] <- inprod(B[d,], lambeta[])
        # Presence spline
        pi_spline[d] <- inprod(B[d,], pibeta[])
        # Dispersion Spline
        alpha_spline[d] <- exp(inprod(B[d,], abeta[]))
    }
    
    for(i in 1:idays){
        # Ones trick
        ones[i] ~ dbern(p[i])
        p[i] <- Lik[i]/C
        Lik[i] <- exp(logLik[i])
        logLik[i] <- (1-z[i])*log(1-w[i]) + 
        #    z[i]*(log(w[i]) + logTPois[i])
            z[i]*(log(w[i]) + LogTruncNB[i])
        
        # For Hurdle
        logit(w[i]) <- pi_spline[DayYearN[i]]
        
        # Truncated Poisson Distribution
        #logTPois[i] <- log(w[i]) + Ni[i]*log(mu[i]) -
        #    mu[i] - log(1 - exp(-mu[i])) - loggam(Ni[i] + 1)
            
        # https://www.cambridge.org/core/services/aop-cambridge-core/content/view/C1FFB63EC7657CD424CE6726F32BC4FE/9781316459515c7_p184-214_CBO.pdf/glms_part_iii_zeroinflated_and_hurdle_models.pdf
        # https://ncss-wpengine.netdna-ssl.com/wp-content/themes/ncss/pdf/Procedures/NCSS/Zero-Inflated_Negative_Binomial_Regression.pdf
        LogTruncNB[i] <-  (1/alpha_spline[DayYearN[i]])*log(u[i])+
            Ni[i]*log(1 - u[i]) + 
            loggam(Ni[i] + 1/alpha_spline[DayYearN[i]]) -
            loggam(1/alpha_spline[DayYearN[i]]) - 
            loggam(Ni[i] + 1) -
            log(1 - (1 + alpha_spline[DayYearN[i]]*mu[i])^(-1/alpha_spline[DayYearN[i]]))
         u[i] <- 1/(1 + alpha_spline[DayYearN[i]]*mu[i])
         
        #https://data.princeton.edu/wws509/notes/countmoments
        ETrunc1[i] <- 1 + alpha_spline[DayYearN[i]]*mu[i]
        ETrunc2[i] <- -1/alpha_spline[DayYearN[i]]
        ETrunc3[i] <- ETrunc1[i]^ETrunc2[i]
        ETruncNB[i] <- mu[i]/mu[i]/max(0.00001, 1 - ETrunc3[i])
        #ExTruncNB[i] <- mu[i]/((1-(1 + 
        #    alpha_spline[DayYearN[i]]*mu[i]))^(-1/alpha_spline[DayYearN[i]]))
        ENi[i] <- w[i]*ETruncNB[i]
        
        # mu is the mean of the Poisson
        # lambda[r] is region-specific intercept
        log(mu[i]) <- lambda_spline[DayYearN[i]] + lambda[FireCentre2N[i]] + 
            b[DayYearN[i]] + inprod(ncovars[i,], betan)
    }
    
    
    for(j in 1:jfires){
        Xijr_star[j] ~ dlnorm(lp2[j], taux[FireCentre2[j]])
        lp1[j] <- inprod(xcovars[j,], betax[, FireCentre2[j]])
        lp2[j] <- mux[FireCentre2[j]] + 
                lp1[j] + 
                gammai[FireCentre2[j]] * b[DayYear[j]]
        
        Xijr1[j] ~ dinterval(Xijr_star[j], dint[j,])
        
        EXi[j] <- exp(lp2[j] + 1/(2*taux[FireCentre2[j]]))
        XLik[j] <- log(plnorm(dint[j, 2], lp2[j], taux[FireCentre2[j]]) -
            plnorm(dint[j, 1], lp2[j], taux[FireCentre2[j]]))
    }
    
    # Priors
    taub <- pow(sigma_b, -2)
    sigma_b ~ dgamma(8, 2)
    
    
    # Splines
    #lambeta ~ dmnorm(zero_vec, lambeta_posdef)
    #pibeta ~ dmnorm(zero_vec, lambeta_posdef)
    #abeta ~ dmnorm(zero_vec, lambeta_posdef)
    lambeta[1] ~ dnorm(0, 1/10)
    pibeta[1] ~ dnorm(0, 1/10)
    abeta[1] ~ dnorm(0, 1/10)
    for(j in 2:J){
        u1[j] ~ dnorm(0, 1/10)
        lambeta[j] <- lambeta[j-1] + u1[j]
        u2[j] ~ dnorm(0, 1/10)
        pibeta[j] <- pibeta[j-1] + u2[j]
        u3[j] ~ dnorm(0, 1/10)
        abeta[j] <- abeta[j-1] + u3[j]
    }
    
    # AR(1) structure for random effects
    phi[1] ~ dnorm(0, 1/5)
    b[1] ~ dnorm(0, taub)
    for(i in 2:ndays){
        b[i] ~ dnorm(phi[1]*b[i-1], taub)
    }

    
    for(p in 1:nxcovar){
		for(r in 1:6){
			betax[r, p] ~ dnorm(0, 1/4)
		}
    }
    
    for(p in 1:nNcovar){
		betan[p] ~ dnorm(0, 1/4)
    }
    
    for(r in 1:6){
        gammai[r] ~ dnorm(hypgamma, 1/5)
        mux[r] ~ dnorm(hypmux, 1/5)
        lambda[r] ~ dnorm(hyplambda, 1/5)
        taux[r] <- pow(sigma_x[r], -2)
        sigma_x[r] ~ dgamma(8, 2)
        #roundprob[r] ~ dbeta(2,2)
        
    }
	hypgamma ~ dnorm(0, 1/5)
	hypmux ~ dnorm(0, 1/5)
	hyplambda ~ dnorm(0, 1/5)
}
"}


wide_jags_AR1pbeta_NiCov_regpi_varsel2 <- {"
model{
    # For the ones trick
    # C chosen so that Lik/c <= 1
    # P(ones == 1) = p, where p is maximized
    # The trick: p is the actual likelihood, we don't care about bernoulli!
    C <- 1000000
    for(d in 1:ndays){
        # Seasonal Affect
        lambda_spline[d] <- inprod(B[d,], lambeta[1:J])
        # Presence spline
		for(r in 1:6){
            pi_spline[r,d] <- inprod(B[d,], pibeta[r,1:J])
		}
        # Dispersion Spline
        alpha_spline[d] <- exp(inprod(B[d,], abeta[1:J]))
    }
    
    for(i in 1:idays){
        # Ones trick
        ones[i] ~ dbern(p[i])
        p[i] <- Lik[i]/C
        Lik[i] <- exp(logLik[i])
        logLik[i] <- (1-z[i])*log(1-w[i]) + 
        #    z[i]*(log(w[i]) + logTPois[i])
            z[i]*(log(w[i]) + LogTruncNB[i])
        
        # For Hurdle
        logit(w[i]) <- pi_spline[FireCentre2N[i], DayYearN[i]]
        
        # Truncated Poisson Distribution
        #logTPois[i] <- log(w[i]) + Ni[i]*log(mu[i]) -
        #    mu[i] - log(1 - exp(-mu[i])) - loggam(Ni[i] + 1)
            
        # https://www.cambridge.org/core/services/aop-cambridge-core/content/view/C1FFB63EC7657CD424CE6726F32BC4FE/9781316459515c7_p184-214_CBO.pdf/glms_part_iii_zeroinflated_and_hurdle_models.pdf
        # https://ncss-wpengine.netdna-ssl.com/wp-content/themes/ncss/pdf/Procedures/NCSS/Zero-Inflated_Negative_Binomial_Regression.pdf
        LogTruncNB[i] <-  (1/alpha_spline[DayYearN[i]])*log(u[i])+
            Ni[i]*log(1 - u[i]) + 
            loggam(Ni[i] + 1/alpha_spline[DayYearN[i]]) -
            loggam(1/alpha_spline[DayYearN[i]]) - 
            loggam(Ni[i] + 1) -
            log(1 - (1 + alpha_spline[DayYearN[i]]*mu[i])^(-1/alpha_spline[DayYearN[i]]))
         u[i] <- 1/(1 + alpha_spline[DayYearN[i]]*mu[i])
         
        #https://data.princeton.edu/wws509/notes/countmoments
        ETrunc1[i] <- 1 + alpha_spline[DayYearN[i]]*mu[i]
        ETrunc2[i] <- -1/alpha_spline[DayYearN[i]]
        ETrunc3[i] <- ETrunc1[i]^ETrunc2[i]
        ETruncNB[i] <- mu[i]/max(0.00001, 1 - ETrunc3[i])
        #ExTruncNB[i] <- mu[i]/((1-(1 + 
        #    alpha_spline[DayYearN[i]]*mu[i]))^(-1/alpha_spline[DayYearN[i]]))
        ENi[i] <- w[i]*ETruncNB[i]
        
        # mu is the mean of the Poisson
        # lambda[r] is region-specific intercept
        log(mu[i]) <- lambda_spline[DayYearN[i]] + lambda[FireCentre2N[i]] + 
            b[DayYearN[i]] + inprod(ncovars[i,], betan)
    }
    
    
    for(j in 1:jfires){
        Xijr_star[j] ~ dlnorm(lp2[j], taux[FireCentre2[j]])
        lp1[j] <- inprod(xcovars[j,], betax[FireCentre2[j],])
        lp2[j] <- mux[FireCentre2[j]] + 
                lp1[j] + 
                gammai[FireCentre2[j]] * b[DayYear[j]]
        
        Xijr1[j] ~ dinterval(Xijr_star[j], dint[j,])
        
        EXi[j] <- exp(lp2[j] + 1/(2*taux[FireCentre2[j]]))
        XLik[j] <- log(plnorm(dint[j, 2], lp2[j], taux[FireCentre2[j]]) -
            plnorm(dint[j, 1], lp2[j], taux[FireCentre2[j]]))
    }
    
    # Priors
    taub <- pow(sigma_b, -2)
    sigma_b ~ dgamma(8, 2)
    
    
    
    # Splines
    u1[1] ~ dnorm(0,1/10)
    u3[1] ~ dnorm(0,1/10)
    lambeta[1] <- u1[1]
    abeta[1] <- u3[1]
    for(r in 1:6){
        u2[r,1] ~ dnorm(0,1/10)
        pibeta[r,1] <- u2[r,1]
    }
    for(j in 2:J){
        u1[j] ~ dnorm(0, 1/10)
        lambeta[j] <- lambeta[j-1] + u1[j]
        u3[j] ~ dnorm(0, 1/10)
        abeta[j] <- abeta[j-1] + u3[j]
        
        for(r in 1:6){
            u2[r, j] ~ dnorm(0, 1/10)
            pibeta[r, j] <- pibeta[r, j-1] + u2[r, j]
        }
    }
    
    # AR(1) structure for random effects
    phi[1] ~ dnorm(0, 1/5)
    b[1] ~ dnorm(0, taub)
    for(i in 2:ndays){
        b[i] ~ dnorm(phi[1]*b[i-1], taub)
    }

    
    for(p in 1:nNcovar){
        betan[p] ~ dnorm(0, 1/4)
    }
    
    pind ~ dbeta(2,8)
    for(p in 1:nxcovar){
		bind[p] ~ dbern(pind)
        for(r in 1:6){
            betaxT[p, r] ~ dnorm(0, 1/4)
            betax[p, r] <- bind[p]*betaxT[p, r]
        }
    }
    
    for(r in 1:6){
        gammai[r] ~ dnorm(hypgamma, 1/5)
        mux[r] ~ dnorm(hypmux, 1/5)
        lambda[r] ~ dnorm(hyplambda, 1/5)
        taux[r] <- pow(sigma_x[r], -2)
        sigma_x[r] ~ dgamma(8, 2)
        #roundprob[r] ~ dbeta(2,2)
        
    }
	hypgamma ~ dnorm(0, 1/5)
	hypmux ~ dnorm(0, 1/5)
	hyplambda ~ dnorm(0, 1/5)
}
"}


















