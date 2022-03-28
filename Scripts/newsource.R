suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(lubridate))
suppressPackageStartupMessages(library(bayestestR))
suppressPackageStartupMessages(library(splines))
suppressPackageStartupMessages(library(rjags))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggmap))
suppressPackageStartupMessages(library(cowplot))
theme_set(theme_bw())
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))

register_google(key = "") # Key has been deleted from my Google API

load.module("glm")
load.module("dic")
# https://www.mrc-bsu.cam.ac.uk/wp-content/uploads/DIC-slides.pdf

# These files were originally made in the commented out 
# sections of JAGSSource.R
bcoutline <- readRDS("Data/Modified/bcoutline2.rds")
FireCentreShape <- readRDS("Data/Modified/shape3.rds")
bc3 <- readRDS("Data/Modified/bc3.rds")

# Plotting maps:
#ggplot() +
#    geom_polygon(aes(x = long, y = lat, colour = id), 
#        data = FireCentreShape, fill = "white") + 
#    geom_point(aes(x = Longitude, y = Latitude, colour = factor(FireCentre2)), 
#        data = filter(bc3, FireYear == 1992))

# Simplify plots for regions
region_vars <- function(jagsdf, varname, plot_it = TRUE, facet = FALSE, vline = NULL,
    newnames = "Region"){
    newdata <- select(jagsdf, starts_with(paste0(varname, "_")))
    if(length(newnames) == ncol(newdata)) {
        names(newdata) <- newnames
        newdata <- gather(newdata, key = "Region", value = "Value")
    } else {
        newdata <- gather(newdata, key = "Region", value = "Value") %>% 
            mutate(Region = substr(Region, nchar(varname) + 2, nchar(varname) + 3))
    }
    
    if(plot_it){
        the_plot <- ggplot(newdata, aes(x = Value, fill = Region)) + geom_density(alpha = 0.3)
        if(facet){
            the_plot <- the_plot + facet_wrap(~ Region, scales = "free")
        } 
        if(!is.null(vline)){
            the_plot <- the_plot + geom_vline(aes(xintercept = vline))
        }
        print(the_plot)
        
    } else {
        return(newdata)
    }
}



rounded_model <- "
model{
    # For the ones trick
    # C chosen so that Lik/c <= 1
    # P(ones == 1) = p, where p is maximized
    # The trick: p is the actual likelihood, we don't care about bernoulli!
    C <- 1000000
    for(i in 1:idays){
        b[i] ~ dnorm(0, taub)
        
        # Seasonal Affect
        lambda_spline[i] <- inprod(B[i,], lambeta[])
        # Presence spline
        pi_spline[i] <- inprod(B[i,], pibeta[])
        
        for(r in 1:6){
            # Ones trick
            ones[i, r] ~ dbern(p[i, r])
            p[i,r] <- Lik[i, r]/C
            Lik[i,r] <- exp(logLik[i, r])
            logLik[i,r] <- (1-z[i,r])*log(1-w[i,r]) + 
                z[i, r]*(log(w[i,r]) + logTPois[i,r])
            
            # For Hurdle
            logit(w[i, r]) <- pi_spline[i]
            
            # Truncated Poisson Distribution
            logTPois[i, r] <- log(w[i,r]) + Ni[i,r]*log(mu[i,r]) -
                mu[i,r] - log(1 - exp(-mu[i,r])) - loggam(Ni[i,r] + 1)
            
            # mu is the mean of the Poisson
            # lambda[r] is region-specific intercept
            log(mu[i,r]) <- lambda_spline[i] + lambda[r] + b[i] #+  betaFN*FWI_N[i, r] + betaWN*WIND_N[i, r]
        }
    }
    
    # Fires must be split up by region... unfortunately.
    # Since JAGS doesn't handle lists and Xjr can't be NA,
    # I can't make a matrix Xj[j, r] (each one has different length)
    # Region 1
    for(j1 in 1:jfires[1]){
        # day is index for 1:ndays, each index repeated Ni times
        Xijr_star1[j1] ~ dlnorm(mux[1] + lp1[j1] + gammai[1]*b[dayr1[j1]], taux[1])
        lp1[j1] <- inprod(xcovar1[j1,], betas)
            
        Xj1[j1] ~ dround(Xijr_star1[j1], 1)
    }
    
    # Region 2
    for(j2 in 1:jfires[2]){
        # day is index for 1:ndays, each index repeated Ni times
        Xijr_star2[j2] ~ dlnorm(mux[2] + lp2[j2] + gammai[2]*b[dayr2[j2]], taux[2])
        lp2[j2] <- inprod(xcovar2[j2,], betas)
            
        Xj2[j2] ~ dround(Xijr_star2[j2], 1)
    }
    
    # Region 3
    for(j3 in 1:jfires[3]){
        # day is index for 1:ndays, each index repeated Ni times
        Xijr_star3[j3] ~ dlnorm(mux[3] + lp3[j3] + gammai[3]*b[dayr3[j3]], taux[3])
        lp3[j3] <- inprod(xcovar3[j3,], betas)
            
        Xj3[j3] ~ dround(Xijr_star3[j3], 1)
    }
    
    # Region 4
    for(j4 in 1:jfires[4]){
        # day is index for 1:ndays, each index repeated Ni times
        Xijr_star4[j4] ~ dlnorm(mux[4] + lp4[j4] + gammai[4]*b[dayr4[j4]], taux[4])
        lp4[j4] <- inprod(xcovar4[j4,], betas)
            
        Xj4[j4] ~ dround(Xijr_star4[j4], 1)
    }
    
    # Region 5
    for(j5 in 1:jfires[5]){
        # day is index for 1:ndays, each index repeated Ni times
        Xijr_star5[j5] ~ dlnorm(mux[5] + lp5[j5] + gammai[5]*b[dayr5[j5]], taux[5])
        lp5[j5] <- inprod(xcovar5[j5,], betas)
            
        Xj5[j5] ~ dround(Xijr_star5[j5], 1)
    }
    
    # Region 6
    for(j6 in 1:jfires[6]){
        # day is index for 1:ndays, each index repeated Ni times
        Xijr_star6[j6] ~ dlnorm(mux[6] + lp6[j6] + gammai[6]*b[dayr6[j6]], taux[6])
        lp6[j6] <- inprod(xcovar6[j6,], betas)
            
        Xj6[j6] ~ dround(Xijr_star6[j6], 1)
    }
    
    taub <- pow(sigma_b, -2)
    sigma_b ~ dgamma(8, 2)
    
    # Splines
    lambeta ~ dmnorm(zero_vec, lambeta_posdef)
    pibeta ~ dmnorm(zero_vec, lambeta_posdef)
    
    
    
    for(p in 1:ncovar){
        betas[p] ~ dnorm(0, 1/5)
    }
    
    for(r in 1:6){
        gammai[r] ~ dnorm(0, 1/5)
        mux[r] ~ dnorm(0, 1/5)
        lambda[r] ~ dnorm(0, 1/5)
        taux[r] <- pow(sigma_x[r], -2)
        sigma_x[r] ~ dgamma(8, 2)
    }
    
    
}"


# Regions
xregion <- function(r){
    regiontext <- "
    # Region ?
    for(j? in 1:jfires[?]){
        # day is index for 1:ndays, each index repeated Ni times
        Xijr_star?[j?] ~ dlnorm(mux[?] + lp?[j?] + gammai[?]*b[dayr?[j?]], taux[?])
        lp?[j?] <- inprod(xcovar?[j?,], betas)
        
        # px? is probabilty of wider regime (digit preference rather than rounding)
        px?[j?] ~ dbern(roundprob[?])
        ind?[j?, 1] <- dint?[j?, 2*px?[j?] + 1]
        ind?[j?, 2] <- dint?[j?, 2*px?[j?] + 2]
        
        Xj1?[j?] ~ dinterval(Xijr_star?[j?], ind?[j?,])
    }
    "
    gsub(pattern = "\\?", replacement = r, x = regiontext)
}

dinterval_mixture <- paste0("
model{
    # For the ones trick
    # C chosen so that Lik/c <= 1
    # P(ones == 1) = p, where p is maximized
    # The trick: p is the actual likelihood, we don't care about bernoulli!
    C <- 1000000
    for(i in 1:idays){
        b[i] ~ dnorm(0, taub)
        
        # Seasonal Affect
        lambda_spline[i] <- inprod(B[i,], lambeta[])
        # Presence spline
        pi_spline[i] <- inprod(B[i,], pibeta[])
        
        for(r in 1:6){
            # Ones trick
            ones[i, r] ~ dbern(p[i, r])
            p[i,r] <- Lik[i, r]/C
            Lik[i,r] <- exp(logLik[i, r])
            logLik[i,r] <- (1-z[i,r])*log(1-w[i,r]) + 
                z[i, r]*(log(w[i,r]) + logTPois[i,r])
            
            # For Hurdle
            logit(w[i, r]) <- pi_spline[i]
            
            # Truncated Poisson Distribution
            logTPois[i, r] <- log(w[i,r]) + Ni[i,r]*log(mu[i,r]) -
                mu[i,r] - log(1 - exp(-mu[i,r])) - loggam(Ni[i,r] + 1)
            
            # mu is the mean of the Poisson
            # lambda[r] is region-specific intercept
            log(mu[i,r]) <- lambda_spline[i] + lambda[r] + b[i] #+  betaFN*FWI_N[i, r] + betaWN*WIND_N[i, r]
        }
    }
    
    # Fires must be split up by region... unfortunately.
    # Since JAGS doesn't handle lists and Xjr can't be NA,
    # I can't make a matrix Xj[j, r] (each one has different length)
    
    
    ", xregion(1), xregion(2), xregion(3), xregion(4), xregion(5), xregion(6), "
    
    taub ~ dgamma(0.01, 0.01)
    sigma_b <- pow(taub, -1/2)
    
    # Splines
    lambeta ~ dmnorm(zero_vec, lambeta_posdef)
    pibeta ~ dmnorm(zero_vec, lambeta_posdef)
    
    
    
    for(p in 1:ncovar){
        betas[p] ~ dnorm(0, 1/5)
    }
    
    for(r in 1:6){
        gammai[r] ~ dnorm(0, 1/5)
        mux[r] ~ dnorm(0, 1/5)
        lambda[r] ~ dnorm(0, 1/5)
        taux[r] ~ dgamma(0.01, 0.01)
        sigma_x[r] <- pow(taux[r], -1/2)
        roundprob[r] ~ dbeta(2,2)
    }
    
}")

# Regions
old_xregion <- function(r){
    regiontext <- "
    # Region ?
    for(j? in 1:jfires[?]){
        # day is index for 1:ndays, each index repeated Ni times
        Xj?[j?] ~ dlnorm(mux[?] + lp?[j?] + gammai[?]*b[dayr?[j?]], taux[?])
        lp?[j?] <- inprod(xcovar?[j?,], betas)
    }
    "
    gsub(pattern = "\\?", replacement = r, x = regiontext)
}

old_model <- paste0("
model{
    # For the ones trick
    # C chosen so that Lik/c <= 1
    # P(ones == 1) = p, where p is maximized
    # The trick: p is the actual likelihood, we don't care about bernoulli!
    C <- 1000000
    for(i in 1:idays){
        b[i] ~ dnorm(0, taub)
        
        # Seasonal Affect
        lambda_spline[i] <- inprod(B[i,], lambeta[])
        # Presence spline
        pi_spline[i] <- inprod(B[i,], pibeta[])
        
        for(r in 1:6){
            # Ones trick
            ones[i, r] ~ dbern(p[i, r])
            p[i,r] <- Lik[i, r]/C
            Lik[i,r] <- exp(logLik[i, r])
            logLik[i,r] <- (1-z[i,r])*log(1-w[i,r]) + 
                z[i, r]*(log(w[i,r]) + logTPois[i,r])
            
            # For Hurdle
            logit(w[i, r]) <- pi_spline[i]
            
            # Truncated Poisson Distribution
            logTPois[i, r] <- log(w[i,r]) + Ni[i,r]*log(mu[i,r]) -
                mu[i,r] - log(1 - exp(-mu[i,r])) - loggam(Ni[i,r] + 1)
            
            # mu is the mean of the Poisson
            # lambda[r] is region-specific intercept
            log(mu[i,r]) <- lambda_spline[i] + lambda[r] + b[i] #+  betaFN*FWI_N[i, r] + betaWN*WIND_N[i, r]
        }
    }
    
    # Fires must be split up by region... unfortunately.
    # Since JAGS doesn't handle lists and Xjr can't be NA,
    # I can't make a matrix Xj[j, r] (each one has different length)
    
    
    ", old_xregion(1), old_xregion(2), old_xregion(3), old_xregion(4), old_xregion(5), old_xregion(6), "
    
    taub <- pow(sigma_b, -2)
    sigma_b ~ dgamma(8, 2)
    
    # Splines
    lambeta ~ dmnorm(zero_vec, lambeta_posdef)
    pibeta ~ dmnorm(zero_vec, lambeta_posdef)
    
    
    
    for(p in 1:ncovar){
        betas[p] ~ dnorm(0, 1/5)
    }
    
    for(r in 1:6){
        gammai[r] ~ dnorm(0, 1/5)
        mux[r] ~ dnorm(0, 1/5)
        lambda[r] ~ dnorm(0, 1/5)
        taux[r] <- pow(sigma_x[r], -2)
        sigma_x[r] ~ dgamma(8, 2)
        #roundprob[r] ~ dbeta(2,2)
    }
    
}")


newnames <- c("betas", "gammai", "mux", "lambda", "sigma_b", "sigma_x",
    "lambeta", "pibeta", "b", "deviance")

newnames_im <- c("betas", "gammai", "mux", "lambda", "sigma_b", "sigma_x",
    "lambeta", "pibeta", "b", "deviance", "roundprob")

round_data_prep <- function(fires, 
    covar_names = c("WIND", "Slope", "Elevation", "ISI", "DC", "DiscoverySize")){
    # howto: pass the year of data, maybe some covar names
    # output is a list with forjags, forinits, and the covar names (might be useful later)
    
    fire.season <<- yday("2003-04-01"):yday("2003-10-01")
    fires <- filter(fires, DayYear %in% fire.season)
    
    Ni <- fires %>% 
        group_by(DayYear, FireCentre2) %>%
        dplyr::select(DayYear, FireCentre2) %>% 
        dplyr::summarise(Count = dplyr::n()) %>% 
        right_join(data.frame(DayYear = 1:365), by = "DayYear") %>% 
        dplyr::select(DayYear, Count, FireCentre2) %>%
        tidyr::spread(data = ., key = FireCentre2, value = Count)
    Ni <- Ni[fire.season, 2:7]
    Ni[is.na(Ni)] <- 0
    
    centreit <- function(x) x - mean(x, na.rm = TRUE)
    
    Xjr1 <- fires %>% 
        filter(!is.na(FinalControlSize), !is.na(FireCentre2)) %>% 
        select(FinalControlSize, FireCentre2, Slope, 
            Elevation, BUI, DC, DMC, DSR, FFMC, FWI, 
            ISI, PCP, RH, TEMP, WIND, DiscoverySize, DayYear) %>% 
        group_by(FireCentre2) %>%
        mutate(Slope = centreit(Slope),
            Elevation = centreit(Elevation),
            BUI = centreit(BUI), DMC = centreit(DMC),
            DC = centreit(DC), DSR = centreit(DSR),
            FFMC = centreit(FFMC), FWI = centreit(FWI),
            ISI = centreit(ISI), PCP = centreit(PCP),
            RH = centreit(RH), TEMP = centreit(TEMP),
            WIND = centreit(WIND),
            DiscoverySize = centreit(DiscoverySize)) %>%
        group_split() %>% 
        lapply(as.list)
    names(Xjr1) <- filter(fires, !is.na(FireCentre2)) %>% 
        group_by(FireCentre2) %>% 
        group_keys() %>% 
        pull() %>% paste("Xj", ., sep = "")
    
    Xjr <- unlist(Xjr1, recursive = FALSE)
    
    # Splines
    B <- bs(fire.season, knots = quantile(fire.season, 0:10/10))
    J <- ncol(B)
    A <- matrix(runif(J^2)*2-1, ncol=J) 
    Sigma <- t(A) %*% A
    
    jagsdata <- list(Ni = Ni, B = B, 
        jfires = sapply(Xjr1, function(x) length(x[[1]])),
        idays = length(fire.season),
        ncovar = length(covar_names),
        ones = matrix(1, ncol = ncol(Ni), nrow = nrow(Ni)),
        z = (Ni > 0)*1,
        lambeta_posdef = Sigma,
        zero_vec = rep(0, J))
    
    
    Xjrs <- unlist(lapply(Xjr1, function(x) x["FinalControlSize"]), recursive = FALSE)
    names(Xjrs) <- substr(names(Xjrs), 1, 3)
    days <- lapply(Xjr1, function(x) x["DayYear"][[1]] - fire.season[1] + 1) 
    names(days) <- paste0("dayr", 1:6)
    
    covars <- lapply(Xjr1, function(x) data.frame(x[covar_names]))
    names(covars) <- paste0("xcovar", 1:6)
    jagsdata <- c(jagsdata, Xjrs, covars, days)
    
    round_inits <- list(Xijr_star1 = jagsdata$Xj1 + 0.001, 
        Xijr_star2 = jagsdata$Xj2 + 0.001,
        Xijr_star3 = jagsdata$Xj3 + 0.001,
        Xijr_star4 = jagsdata$Xj4 + 0.001,
        Xijr_star5 = jagsdata$Xj5 + 0.001,
        Xijr_star6 = jagsdata$Xj6 + 0.001
    )
    
    list(forjags = jagsdata, forinits = round_inits, 
        covar_names = covar_names)
}


ic_heap_ab <- function(x, version = "a"){
    ic_switch <- matrix(ncol = 3, byrow = TRUE,
        data = c(0, 0, 0.05,
            0.1, 0, 0.25,
            0.5, 0.25, 0.75,
            1, 0.75, 1.25,
            1.5, 1.25, 1.75, 
            2, 1.75, 2.25))
    
    if(version == "a"){
        ic_switch <- matrix(ncol = 3, byrow = FALSE,
            data = c(seq(0, 5, 0.1), 
                c(0, seq(0.05, 4.95, 0.1)),
                seq(0.05, 5.05, 0.1))
        )
    }
    
    newx <- matrix(data = rep(x, 2), ncol = 2, byrow = FALSE)
    
    for(i in 1:nrow(newx)){
        if(newx[i, 1] %in% ic_switch[, 1]){
            this <- which(ic_switch[, 1] == newx[i, 1])
            newx[i, ] <- ic_switch[this, 2:3]
        } else {
            newx[i, ] <- newx[i, ] + c(-0.05, 0.05)
        }
    }
    
    newx
}

mixture_data_prep <- function(fires, 
    covar_names = c("WIND", "Slope", "Elevation", "ISI", "DC", "DiscoverySize")){
    
    x <- round_data_prep(fires = fires, covar_names = covar_names)
    
    x$forjags$dint1 <- cbind(ic_heap_ab(x$forjags$Xj1), 
        ic_heap_ab(x$forjags$Xj1, version = "b"))
    x$forjags$dint2 <- cbind(ic_heap_ab(x$forjags$Xj2), 
        ic_heap_ab(x$forjags$Xj2, version = "b"))
    x$forjags$dint3 <- cbind(ic_heap_ab(x$forjags$Xj3), 
        ic_heap_ab(x$forjags$Xj3, version = "b"))
    x$forjags$dint4 <- cbind(ic_heap_ab(x$forjags$Xj4), 
        ic_heap_ab(x$forjags$Xj4, version = "b"))
    x$forjags$dint5 <- cbind(ic_heap_ab(x$forjags$Xj5), 
        ic_heap_ab(x$forjags$Xj5, version = "b"))
    x$forjags$dint6 <- cbind(ic_heap_ab(x$forjags$Xj6), 
        ic_heap_ab(x$forjags$Xj6, version = "b"))
    
    
    
    x$forjags$Xj11 = rep(1, x$forjags$jfires[1])
    x$forjags$Xj12 = rep(1, x$forjags$jfires[2])
    x$forjags$Xj13 = rep(1, x$forjags$jfires[3])
    x$forjags$Xj14 = rep(1, x$forjags$jfires[4])
    x$forjags$Xj15 = rep(1, x$forjags$jfires[5])
    x$forjags$Xj16 = rep(1, x$forjags$jfires[6])
    
    obsize <- list(Xj1 = x$forjags$Xj1, Xj2 = x$forjags$Xj2, Xj3 = x$forjags$Xj3, 
        Xj4 = x$forjags$Xj4, Xj5 = x$forjags$Xj5, Xj6 = x$forjags$Xj6)
    
    x$forjags$Xj1 <- NULL
    x$forjags$Xj2 <- NULL
    x$forjags$Xj3 <- NULL
    x$forjags$Xj4 <- NULL
    x$forjags$Xj5 <- NULL
    x$forjags$Xj6 <- NULL
    
    list(forjags = x$forjags, forinits = x$forinits, 
        covar_names = x$covar_names, obsize = obsize)
}

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

#cat("dinterval_mixture | mixture_data_prep(fires, names) | region_vars() | melt.mcmc.list.rename()")

cat("models: rounded_model, dinterval_mixture, old_model\n")
cat("data: bc3, bcoutline, FireCentreShape, newnames, newnames_im\n    round_data_prep(), mixture_data_prep()\n")
cat("helpers: region_vars(), melt.mcmc.list.rename()\n")
















