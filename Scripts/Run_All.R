# Setting up the models

library(runjags)

re_source <- function() source("Scripts/modelCompendium.R")
re_source()

bc3 <- readRDS("Data/Modified/bc3.rds")
precip <- readRDS("Data/Modified/regionalWeather.rds")

mychains <- 3
myiter <- 500
myadapt <- 200
myburnin <- 300
myupdate <- 50
mythin <- 20
FireYears <- unique(bc3$FireYear)
FireYears <- FireYears[FireYears <= 1995]
#myyear <- 1990
mycause <- "LTG"
mymodelnumber <- 21
mymodelname <- c("mix_jags_AR1pbeta", "mix_jags_AR1pbeta_sep", #2
    "wide_jags_AR1pbeta", "wide_jags_AR1pbeta_sep", #4
    "mix_jags_AR1", "mix_jags_AR2", #6
    "mix_jags_AR1pbeta_gspline", "mix_jags_AR1dahyp", #8
    "mix_jags_hyper", "wide_jags_hyper", #10
    "mix_jags_base", "wide_jags_base", #12
    "mix_jags_AR1pbeta_NiCov", #13
    "wide_jags_AR1pbeta_NiCov",  #14
    "wide_jags_AR1pbeta_NiCov_regpi", #15
    "wide_jags_AR1pbeta_NiCov_regpi_sep", # 16
    "wide_jags_AR1rbeta_NiCov", # 17
    "wide_jags_regre_pbeta_NiCov", # 18
    "wide_jags_AR1pbeta_NiCov_regpi_varsel", # 19
    "wide_jags_AR1rbeta_NiCov_varsel", # 20
    "wide_jags_AR1rbeta_NiCov2", # 21
    "wide_jags_AR1pbeta_NiCov_regpi_varsel2")[mymodelnumber]
mymodel <- eval(as.symbol(mymodelname))

# Will get warnings about failure to set trace. Ignore them.
imodelmonitors <- c("gammai", "mux", "lambda", "sigma_b", 
    "sigma_x", "betax", "betan",
    "lambeta", "pibeta", "abeta", "b", "phi", "hypgamma", "hypmux", "hyplambda",
    "roundprob",
    "logLik", "XLik", "EXi", "ENi",
    "deviance")


modlist <- vector(mode = "list", length = length(FireYears))
grdf <- data.frame(coda = NA, miness = NA, meaness = NA, 
    nparams = NA, numOut1.1 = NA, numOut1.15 = NA, numCIOut1.1 = NA,
    numCIOut1.15 = NA, mpsrf = NA, year = NA)
nxdf <- data.frame(data = NA, Xrmse = NA, Nrmse = NA,
    Xwaic = NA, Nwaic = NA, WAIC = NA, year = NA)


modfilename <- paste0("Big Models/allyears3b_", 
    mymodelname, mycause, "_dflist.rds")
modfilename2 <- paste0("Data/Models/allyears3b_", 
    mymodelname, mycause, "_dflist.rds")
grfilename <- paste0("Data/Models/allyears3b_", 
    mymodelname, mycause, "_gr.rds")
nxfilename <- paste0("Data/Models/allyears3b_", 
    mymodelname, mycause, "_nx.rds")

if(paste0("allyears3_", mymodelname, mycause, "_nx.rds") %in% 
        list.files("Data/Models")){
    modlist <- list(readRDS(modfilename2))
    grdf <- readRDS(grfilename)
    nxdf <- readRDS(nxfilename)
}

t0 <- Sys.time()

thesetimes <- c()

for(i in seq_along(FireYears)[FireYears > max(nxdf$year, na.rm = TRUE)]){
    tryCatch({
        t1 <- Sys.time()
        
        imodeldata <- mix_data_prep_2(ncovar = precip, 
            units = ifelse(FireYears[i] > 1975, 
                "metric", "imperial"), 
            fires = dplyr::filter(bc3, FireYear == FireYears[i], 
                Cause == mycause,
                ISI >= 0, FWI >= 0, Slope > -80, 
                WIND >= 0, TEMP > -50, RH >= 0), 
            mycov = c("ISI", "FWI", "Slope", "WIND", "TEMP", "RH"))
        
        if(substr(mymodelname, 1, 3) == "rou") {
            imodeldata$jdata$dint <- imodeldata$jdata$dint[, 1:2]
        } else if(substr(mymodelname, 1, 3) == "wid") {
            imodeldata$jdata$dint <- imodeldata$jdata$dint[, 3:4]
        }
        
        imodelmodel <- jags.model(file = textConnection(mymodel), 
            data = imodeldata$jdata, inits = imodeldata$idata, 
            n.chains = mychains, n.adapt = myadapt, quiet = TRUE)
        update(imodelmodel, myupdate, progress.bar = "none")
        imodelcoda <- coda.samples(imodelmodel, n.iter = myiter, 
            variable.names = imodelmonitors, thin = mythin,
            progress.bar = "text")

        imodeldf <- melt.mcmc.list.rename(imodelcoda)
        imodeldf$year <- FireYears[i]
        
        modlist[[i]] <- imodeldf %>% select(-starts_with("ENi"), 
            -starts_with("EXi"),
            -starts_with("logLik"), - starts_with("XLik"))
        
        imodel_NX <- extract_NX(jdf = imodeldf, jdata = imodeldata)
        imodel_NX$year <- FireYears[i]
        imodel_NX$data <- mymodelname
        #imodel_gr <- extract_gr(imodelcoda)
        #imodel_gr$year <- FireYears[i]
        #imodel_gr$coda <- mymodelname
        
        #grdf <- bind_rows(grdf, imodel_gr)
        nxdf <- bind_rows(nxdf, imodel_NX)
        
        thesetimes[i] <- difftime(Sys.time(), t1, units = "mins") %>% 
            round(3)
        meantime <- mean(thesetimes, na.rm = TRUE) %>% round(3)
        estleft <- round(meantime*(length(FireYears) - i), 3)
        
        print(paste0(FireYears[i], " (loop ", i, " of ", 
            length(FireYears), ") took ", 
            thesetimes[i], " mins. Approx ",
            estleft, " mins remaining."))
        
        
        allmods <- bind_rows(modlist) %>% select(-starts_with("ENi"), 
            -starts_with("EXi"),
            -starts_with("logLik"), - starts_with("XLik"))
        
        saveRDS(allmods, file = modfilename)
        #saveRDS(grdf, file = grfilename)
        saveRDS(nxdf, file = nxfilename)
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
Sys.time() - t0

saveRDS(allmods, file = modfilename)
#saveRDS(grdf, file = grfilename)
saveRDS(nxdf, file = nxfilename)

#grdf
nxdf

if(FALSE){
    
    allmods %>% 
        select(starts_with("gammai"), year) %>% 
        gather(name, gamma, -year) %>% 
        separate(name, sep = "\\_", into = c("name", "region")) %>% 
        select(-name) %>% 
        mutate(region = as.numeric(as.character(region))) %>% 
        group_by(region, year) %>% 
        summarise(lo = quantile(gamma, 0.055), med = median(gamma), 
            hi = quantile(gamma, 0.945)) %>% 
        mutate(sig = ifelse(lo > 0, "positive", 
            ifelse(hi < 0, "negative", "neutral"))) %>% 
        ggplot(aes(x = year, colour = sig)) + 
        geom_errorbar(aes(ymin = lo, ymax = hi)) + 
        geom_point(aes(y = med)) +
        facet_wrap(~ region) +
        geom_hline(aes(yintercept = 0)) +
        theme_bw() +
        scale_colour_manual(values = c("red", "black", "green"))
    
    allmods %>% 
        select(starts_with("mux"), year) %>% 
        gather(name, mux, -year) %>% 
        separate(name, sep = "\\_", into = c("name", "region")) %>% 
        select(-name) %>% 
        mutate(region = as.numeric(as.character(region))) %>% 
        group_by(region, year) %>% 
        summarise(lo = quantile(mux, 0.055), med = median(mux), 
            hi = quantile(mux, 0.945)) %>% 
        ggplot(aes(x = year)) + 
        geom_errorbar(aes(ymin = lo, ymax = hi)) + 
        geom_point(aes(y = med)) +
        facet_wrap(~ region) +
        theme_bw()
    
    allmods %>% 
        select(starts_with("lambda"), year) %>% 
        gather(name, lambda, -year) %>% 
        separate(name, sep = "\\_", into = c("name", "region")) %>% 
        select(-name) %>% 
        mutate(region = as.numeric(as.character(region))) %>% 
        group_by(region, year) %>% 
        summarise(lo = quantile(lambda, 0.055), med = median(lambda), 
            hi = quantile(lambda, 0.945)) %>% 
        ggplot(aes(x = year)) + 
        geom_errorbar(aes(ymin = lo, ymax = hi)) + 
        geom_point(aes(y = med)) +
        facet_wrap(~ region) +
        theme_bw()
    
    
    allmods %>% 
        select(starts_with("betax"), year) %>% 
        gather(name, betax, -year) %>% 
        separate(name, sep = "\\_", into = c("name", "region")) %>% 
        select(-name) %>% 
        mutate(region = as.numeric(as.character(region))) %>% 
        group_by(region, year) %>% 
        summarise(lo = quantile(betax, 0.055), med = median(betax), 
            hi = quantile(betax, 0.945)) %>% 
        mutate(sig = ifelse(lo > 0, "positive", 
            ifelse(hi < 0, "negative", "neutral"))) %>% 
        ggplot(aes(x = year, colour = sig)) + 
        geom_errorbar(aes(ymin = lo, ymax = hi)) + 
        geom_point(aes(y = med)) +
        facet_wrap(~ region, labeller = labeller(region = imodeldata$mycov)) +
        theme_bw() +
        scale_colour_manual(values = c("red", "black", "green"))
    
    allmods %>% 
        select(starts_with("sigma_x"), year) %>% 
        gather(name, sigma_x, -year) %>% 
        separate(name, sep = "\\_", into = c("name", "name2", "region")) %>% 
        select(-name, -name2) %>% 
        mutate(region = as.numeric(as.character(region))) %>% 
        group_by(region, year) %>% 
        summarise(lo = quantile(sigma_x, 0.055), med = median(sigma_x), 
            hi = quantile(sigma_x, 0.945)) %>% 
        ggplot(aes(x = year)) + 
        geom_errorbar(aes(ymin = lo, ymax = hi)) + 
        geom_point(aes(y = med)) +
        facet_wrap(~ region) +
        theme_bw() 
    
    allmods %>% 
        select(sigma_b, year) %>% 
        group_by(year) %>% 
        summarise(lo = quantile(sigma_b, 0.055), med = median(sigma_b), 
            hi = quantile(sigma_b, 0.945)) %>% 
        ggplot(aes(x = year)) + 
        geom_errorbar(aes(ymin = lo, ymax = hi)) + 
        geom_point(aes(y = med)) +
        theme_bw() 
    
    allmods %>% 
        select(phi, year) %>% 
        group_by(year) %>% 
        summarise(lo = quantile(phi, 0.055), med = median(phi), 
            hi = quantile(phi, 0.945)) %>% 
        ggplot(aes(x = year)) + 
        geom_errorbar(aes(ymin = lo, ymax = hi)) + 
        geom_point(aes(y = med)) +
        theme_bw() 
    
}










































