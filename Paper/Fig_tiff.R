# Figures for Assessing Dependence Between Frequency and Severity Through Shared Random Effects

library(here)
source(here("Scripts", "modelCompendium.R"))

bc3 <- readRDS(here("Data/Modified", "bc3.rds"))

library(ggmap)
library(cowplot)
library(tidyr)
library(stringr)

theme_set(theme_bw())

bcmap <- readRDS("WRAEES Presentation/bcmap.rds")
bcmap2 <- readRDS("WRAEES Presentation/bcmap2.rds")
bcfort <- readRDS("WRAEES Presentation/bcfort.rds")

fivepointtwo <- function() list(
    theme(axis.text = element_text(size = 4), 
        plot.title = element_text(size = 8), 
        axis.title = element_text(size = 6),
        legend.text = element_text(size = 4),
        legend.title = element_text(size = 4), 
        legend.key.size = unit(0.1, units = "in"), 
        strip.text = element_text(size = 5, 
            margin = margin(1/2,1/2,1/2,1/2, "mm")))
)


# Fig1 -------------------------------------------------------------------------

mids_df <- bcfort %>% 
    mutate(id = as.numeric(id) + 1) %>%
    group_by(id) %>% 
    summarise(lat = mean(lat), long = mean(long))

mids_df$lat[c(1, 2, 6)] <- c(50.5, 50.5, 54.5)
mids_df$long <- c(-117, -120, -123, -127, -123, -128)

#mymap2 <- get_map(location = rev(c(50.0024201,-127.0292006)), zoom = 4)
#saveRDS(mymap2, "bcmap2.rds")


mymap1 <- ggplot() + #ggmap(bcmap2) + 
    geom_polygon(data = bcfort, alpha = 0.6,
        mapping = aes(x = long, y = lat, fill = factor(group)), colour = 1) +
    labs(fill = "Fire Centre", title = "Fire Regions",
        y = "Latitude", x = "Longitude") +
    coord_map(xlim = c(-140, -112), ylim = c(47, 60.25)) +
    scale_fill_brewer(type = "qual", palette = "Dark2") + 
    geom_text(data= mids_df, aes(x = long, y = lat, label = id), size = 5) +
    theme_void() + 
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
    fivepointtwo() +
    theme(axis.title.y = element_text(angle = 90))

summaries <- bc3 %>% 
    filter(!is.na(FireCentre2)) %>% 
    group_by(FireYear, FireCentre2) %>% 
    summarise(count = n(), sizes = sum(FinalControlSize, na.rm = TRUE),
        msizes = mean(FinalControlSize, na.rm = TRUE))

freqsev <- plot_grid(
    ggplot(summaries, aes(x = FireYear, y = count, fill = factor(FireCentre2))) + 
        geom_col() +
        scale_fill_brewer(type = "qual", palette = "Dark2") + 
        theme(legend.position = "none") +
        labs(x = "Year", y = "Count", title = "Frequency") + 
        scale_x_continuous(breaks = seq(1945, 2015, 5)) +
        fivepointtwo(),
    ggplot(summaries, aes(x = FireYear, y = sizes, fill = factor(FireCentre2))) + 
        geom_col(na.rm = TRUE) +
        scale_fill_brewer(type = "qual", palette = "Dark2") + 
        theme(legend.position = "none") +
        labs(title = "Size, in Hectares (log 10 scale)", x = "Year", y = "Total Area Burned") + 
        scale_x_continuous(breaks = seq(1945, 2015, 5)) +
        scale_y_log10() +
        fivepointtwo(),
    ncol = 1
)  

gg_1 <- plot_grid(freqsev, nrow = 1, rel_widths = c(1,1.75))

tiff(filename = "Paper/Figures/Fig1.tif", width = 5.2, height = 3, units = "in", res = 600)
freqsev
dev.off()

pdf(file = "Paper/Joint_Count_Files/firezones-1.pdf", width = 5.2, height = 3)
freqsev + theme_bw()
dev.off()








# Fig2 -------------------------------------------------------------------------
#```{r rounding, fig.height = 2.5, fig.cap="\\label{rounding}Evidence of rounding and unit change. Prior to 1985, 75\\% of all fires were recorded as being exactly 0.1ha. The decrease in proportion coincides with the advent of fires recorded as 0ha. In 1975, Canada switched to the metric system. This means that fires were approximated as being 0.5ha rather tham 0.4ha (1 acre). Note that the proportion of 0.1ha fires did not change since 1/10th of a hectare and 1/4th of an acre are both coarsening values."}
tiff(filename = "Paper/Figures/Fig2.tif", width = 5.2, height = 2, units = "in", res = 600)
sml <- bc3 %>% mutate(Size = FinalControlSize) %>%
    filter(Size <= 2, Year < 1996)#, Size > 0)

bigsmall1 <- table(sml$Size)
bigsmall <- bigsmall1[bigsmall1 > 100]
smallsmall <- bigsmall1[bigsmall1 <= 100]

AllDigs <-  sapply(names(bigsmall1), function(x){
    sml %>% 
        group_by(FireYear) %>% 
        summarise(DigitPref = mean(FinalControlSize == x)) %>%
        pull(DigitPref)
}) %>% 
    as.data.frame() %>%
    mutate(FireYear = unique(sml$FireYear)) %>% 
    gather(Digit, Pref, -FireYear) %>% 
    mutate(Digit = ordered(Digit))

#AllDigs %>% 
#    #filter(Digit %in% seq(0,2,0.1)) %>% 
#    filter(Digit %in% c(1,1.5,2)) %>% 
#    ggplot(aes(x = FireYear, y = Pref, colour = Digit)) + 
#        geom_smooth(se = FALSE) +
#        coord_cartesian(ylim = c(0,0.3))

cowplot::plot_grid(
    AllDigs %>% 
        mutate(is0.1 = factor(ifelse(Digit == 0.1, 0.1, "Other"), 
            levels = rev(c(0.1, "Other")), ordered = TRUE)) %>% 
        group_by(is0.1, FireYear) %>% 
        summarise(Pref = sum(Pref)) %>% 
        ggplot(aes(x = FireYear, y = Pref, fill = is0.1)) + 
        theme_bw() + 
        geom_area(colour = "white") +
        scale_fill_manual(values = c("#cf6126", "#2593CF")) +
        labs(fill = "Recorded Size", x = "Year", y = "Proportion",
            title = "Proportion of 0.1 Hectare Fires") +
        annotate("text", x = c(1975, 1975), y = 1-c(0.8, 0.125), 
            label = c("Exactly 0.1ha", "All Other Values"),
            colour = "white", size = 3) +
        theme(legend.position = "none") +
        fivepointtwo(),
    
    AllDigs %>% 
        filter(Digit %in% c(0.3,0.4,0.5)) %>% 
        ggplot(aes(x = FireYear, y = Pref, fill = Digit)) + 
        geom_area(colour = 1) + #, stat = "identity", width = 1) +
        theme_bw() + 
        labs(x = "Year", y = "Preference", fill = "Size (hectares)",
            title = "0.3 to 0.5 Hectare Fires") +
        geom_vline(aes(xintercept = 1975), colour = "red", size = 0.4) +
        theme(legend.position = "bottom") +
        scale_x_continuous(breaks = seq(1945, 2015, 10), 
            minor_breaks = seq(1945, 2020, 1)) +
        fivepointtwo(),# +
    #coord_cartesian(xlim = c(1960,1990)),
    nrow = 1
)

dev.off()





# Fig3 -------------------------------------------------------------------------
#```{r gammaltg, fig.cap="\\label{gammaltg}Posterior median estimates and 90% Credible Intervals for the linking parameter in the model for lightning-caused fires. There are very few years where the credible interval was entirely above 0, but many where it was below. This indicates that the size and count of lightning-caused fires tends to have a negative association. Days with fewer fires had larger fires, and days with more fires had smaller fires."}
ltgest <- readRDS("Big_Models/allyears3b_wide_jags_AR1rbeta_NiCov2LTG_dflist.rds")
gg_ltg <- ltgest %>% 
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
    fivepointtwo() +
    theme_bw() +
    scale_colour_manual(values = c("red", "black", "green")) +
    scale_x_continuous(breaks = seq(1950, 2000, 5), minor_breaks = 1950:2000) +
    theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    labs(y = expression(gamma[r]), x = "Year")
tiff(filename = "Paper/Figures/Fig3.tif", width = 5.2, height = 3, units = "in", res = 600)
gg_ltg
dev.off()

pdf(file = "Paper/Figures/gammaltg.pdf", width = 5.2, height = 3)
gg_ltg
dev.off()







# Fig4 -------------------------------------------------------------------------
#```{r gammaper, fig.cap="\\label{gammaper}Posterior median estimates and 90\\% Credible Intervals for the linking parameter in the model for person-caused fires. In contrast to lightning-caused fires, person-caused fires tend to have positive association."}
perest <- readRDS("Big_Models/allyears3b_wide_jags_AR1pbeta_NiCov_regpi_varsel2PER_dflist.rds")
gg_per <- perest %>% 
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
    fivepointtwo() +
    theme_bw() +
    scale_colour_manual(values = c("red", "black", "green")) +
    scale_x_continuous(breaks = seq(1950, 2000, 5), minor_breaks = 1950:2000) +
    theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    labs(y = expression(gamma[r]), x = "Year")
tiff(filename = "Paper/Figures/Fig4.tif", width = 5.2, height = 3, units = "in", res = 600)
gg_per
dev.off()

pdf(file = "Paper/Figures/gammaper.pdf", width = 5.2, height = 3)
gg_per
dev.off()






# Fig5 -------------------------------------------------------------------------
#```{r hypers, fig.height=4.5, fig.width=11, fig.cap="\\label{hypers}Hyperparameters for the mean of the size and count. By the construction of the model, the regional parameters are deviations from these values, with the hyperparameters representing the mean over all regions. All of these parameters are on the log scale."}
hypmu_ltg <- ltgest %>% select(hypmux, year) %>% 
    group_by(year) %>%
    summarise(median = median(hypmux), 
        qlo = quantile(hypmux, 0.025),
        qhi = quantile(hypmux, 0.975)) %>% 
    ggplot(aes(x = year, y = median)) + 
        geom_errorbar(aes(ymin = qlo, ymax = qhi), size = 0.4) + 
        geom_point(size = 0.4) + #geom_line() +
        labs(x = "Year", y = "Hyperparameter", 
            title = "Hyperparameter for mean fire size - Lightning") +
        fivepointtwo()
hyplam_ltg <- ltgest %>% select(hyplambda, year) %>% 
    group_by(year) %>%
    summarise(median = median(hyplambda), 
        qlo = quantile(hyplambda, 0.025),
        qhi = quantile(hyplambda, 0.975)) %>% 
    ggplot(aes(x = year, y = median)) + 
        geom_errorbar(aes(ymin = qlo, ymax = qhi), size = 0.4) + 
        geom_point(size = 0.4) + #geom_line() +
        labs(x = "Year", y = "Hyperparameter", 
            title = "Hyperparameter mean fire count - Lightning") +
        fivepointtwo()
hypmu_per <- perest %>% select(hypmux, year) %>% 
    group_by(year) %>%
    summarise(median = median(hypmux), 
        qlo = quantile(hypmux, 0.025),
        qhi = quantile(hypmux, 0.975)) %>% 
    ggplot(aes(x = year, y = median)) + 
        geom_errorbar(aes(ymin = qlo, ymax = qhi), size = 0.4) + 
        geom_point(size = 0.4) + #geom_line() +
        labs(x = "Year", y = "Hyperparameter", 
            title = "Hyperparameter for mean fire size - Human") +
        fivepointtwo()
hyplam_per <- perest %>% select(hyplambda, year) %>% 
    group_by(year) %>%
    summarise(median = median(hyplambda), 
        qlo = quantile(hyplambda, 0.025),
        qhi = quantile(hyplambda, 0.975)) %>% 
    ggplot(aes(x = year, y = median)) + 
        geom_errorbar(aes(ymin = qlo, ymax = qhi), size = 0.4) + 
        geom_point(size = 0.4) +# geom_line() +
        labs(x = "Year", y = "Hyperparameter", 
            title = "Hyperparameter for mean fire count - Human") +
        fivepointtwo()

tiff(filename = "Paper/Figures/Fig5.tif", width = 5.2, height = 2.5, units = "in", res = 600)
plot_grid(hypmu_ltg, hyplam_ltg, hypmu_per, hyplam_per, nrow = 2)
dev.off()

pdf(file = "Paper/Joint_Count_Files/hypers-1.pdf", width = 5.2, height = 2.5)
plot_grid(hypmu_ltg, hyplam_ltg, hypmu_per, hyplam_per, nrow = 2)
dev.off()


# Fig6 -------------------------------------------------------------------------
#```{r spatialcor, fig.height = 4.25, fig.cap="\\label{spatialcor}Correlation of the estimated linking parameter across space and time. The linking parameter in one year does not appear to predict the linking parameter in the next year. The linking parameter for lightning-caused fires is correlated for regions 1 and 2 (which are both interior from the coast and with similar forest cover) and regions 4 and 5 (which have no clear similarities). There is province-wide correlation among the linking parameters for person-caused fires."}

gamper <- perest %>% select(starts_with("gammai"), year) %>% 
    group_by(year) %>%
    summarise_all(.fun = median) %>%
    select(-year)
names(gamper) <- paste0("gamma_", 1:6)
gammalabs <- c(bquote(hat(gamma)[1]), bquote(hat(gamma)[2]), bquote(hat(gamma)[3]),
    bquote(hat(gamma)[4]), bquote(hat(gamma)[5]), bquote(hat(gamma)[6]))
gpergg <- GGally::ggcorr(gamper, label = TRUE, digits = 2, legend.size = 6,
    label_size = 2, size = -2, colour = "white", legend.position = "none") +
    labs(title = bquote("Spatial Cor. of "~hat(gamma)[r]))

gammadf <- data.frame(x = 1:6, y = 1:6, labs = as.character(gammalabs))
gpergg <- gpergg + geom_text(aes(x = x, y = x, label = labs), data = gammadf, parse = TRUE)


lamper <- ltgest %>% select(starts_with("gammai"), year) %>% 
    group_by(year) %>%
    summarise_all(.fun = median) %>%
    select(-year)
names(lamper) <- paste0("gamma_", 1:6)
gammalabs <- c(bquote(hat(gamma)[1]), bquote(hat(gamma)[2]), bquote(hat(gamma)[3]),
    bquote(hat(gamma)[4]), bquote(hat(gamma)[5]), bquote(hat(gamma)[6]))
gltggg <- GGally::ggcorr(lamper, label = TRUE, digits = 2, legend.size = 3,
    label_size = 2, size = 0, colour = "white", legend.position = "none") +
    labs(title = bquote("Spatial Cor. of "~hat(gamma)[r]))

gammadf <- data.frame(x = 1:6, y = 1:6, labs = as.character(gammalabs))
gltggg <- gltggg + geom_text(aes(x = x, y = x, label = labs), data = gammadf, parse = TRUE)


ltgacf <- lapply(1:6,function(x) {
    a <- pacf(lamper[,x], plot = FALSE)
    data.frame(ACF = a$acf, Lag = a$lag, gamma = x)
}) %>% 
    bind_rows() %>% 
    ggplot(aes(xend = Lag, x = Lag, yend = 0, y = ACF)) + 
    geom_segment(size = 0.75) + facet_wrap(~ gamma, nrow = 1) +
    geom_hline(yintercept = 0) +
    geom_hline(yintercept = 1.96/sqrt(nrow(lamper)), 
        linetype = "dashed", colour = "blue")+
    geom_hline(yintercept = -1.96/sqrt(nrow(lamper)), 
        linetype = "dashed", colour = "blue") +
    scale_x_continuous(breaks = seq(0, 40, 5), minor_breaks = seq(0,40,1)) +
    labs(title = bquote("Autocorrelations for "~hat(gamma)[r]*" (Lightning)"))
peracf <- lapply(1:6,function(x) {
    a <- pacf(gamper[,x], plot = FALSE)
    data.frame(ACF = a$acf, Lag = a$lag, gamma = x)
}) %>% 
    bind_rows() %>% 
    ggplot(aes(xend = Lag, x = Lag, yend = 0, y = ACF)) + 
    geom_segment(size = 0.75) + facet_wrap(~ gamma, nrow = 1) +
    geom_hline(yintercept = 0) +
    geom_hline(yintercept = 1.96/sqrt(nrow(lamper)), 
        linetype = "dashed", colour = "blue")+
    geom_hline(yintercept = -1.96/sqrt(nrow(lamper)), 
        linetype = "dashed", colour = "blue") +
    scale_x_continuous(breaks = seq(0, 40, 5), minor_breaks = seq(0,40,1)) +
    labs(title = bquote("Autocorrelations for "~hat(gamma)[r]*" (Person)"))

tiff(filename = "Paper/Figures/Fig6.tif", width = 5.2, height = 2.5, units = "in", res = 600)
plot_grid(ltgacf + fivepointtwo(), gltggg + fivepointtwo(), 
    peracf + fivepointtwo(), gpergg + fivepointtwo(), 
    nrow = 2, rel_widths = c(3, 1))
dev.off()

pdf(file = "Paper/Joint_Count_Files/spatialcor-1.pdf", width = 5.2, height = 2.5)
plot_grid(ltgacf + fivepointtwo(), gltggg + fivepointtwo(), 
    peracf + fivepointtwo(), gpergg + fivepointtwo(), 
    nrow = 2, rel_widths = c(3, 1))
dev.off()





# Fig7 -------------------------------------------------------------------------
#```{r sigestl, fig.cap="\\label{sigestl}The variance parameters for lightning-caused fires. there is no apparent trend over time, except perhaps a slight decrease in region 5 and in the random effect variance, with an increase in the AR(1) parameter. The AR(1) parameter is constrained to be in the interval [-1,1]."}
# gamma
gg1 <- ltgest %>% 
    select(year, starts_with("sigma_x")) %>%
    group_by(year) %>%
    summarise_all(.funs = quantile, prob = c(0.025, 0.5, 0.975)) %>%
    mutate(fun = c("qlo", "median", "qhi")) %>%
    tidyr::pivot_longer(cols = 2:7, names_to = "Region") %>%
    mutate(Region = gsub("sigma_x_", "", Region)) %>%
    tidyr::pivot_wider(names_from = "fun", values_from = value) %>% 
    ggplot(aes(x = as.numeric(year), y = median)) + 
        geom_errorbar(mapping = aes(ymin = qlo, ymax = qhi), size = 0.4) +
        geom_point(size = 0.4) + #geom_smooth(method = "loess", se = FALSE) + 
        facet_wrap(~Region, nrow = 2, 
            labeller = labeller(Region = label_both)) +
        theme_bw() +
        labs(x = "year", y = bquote(gamma[r]), 
            title = "LogNormal Variance - Lightning Caused Fires",
            colour = '"Significance"') +
        geom_hline(aes(yintercept = 0), colour = rgb(0,0,0,0)) +
        theme(legend.position = "none", 
            axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
        scale_x_continuous(breaks = seq(1950,2015,10)) +
        fivepointtwo()

gg2 <- ltgest %>% 
    select(year, sigma_b) %>%
    group_by(year) %>%
    summarise(qlo = quantile(sigma_b, 0.025),
        median = median(sigma_b),
        qhi = quantile(sigma_b, 0.975)) %>%
    ggplot(aes(x = year, y = median)) + 
        geom_errorbar(aes(ymin = qlo, ymax = qhi), size = 0.4) +
        geom_point(size = 0.4) +
        labs(x = "Year", y = bquote(sigma[b]), title = "R.E. Variance") +
        fivepointtwo()

gg3 <- ltgest %>% 
    select(year, phi) %>%
    group_by(year) %>%
    summarise(qlo = quantile(phi, 0.025),
        median = median(phi),
        qhi = quantile(phi, 0.975)) %>%
    ggplot(aes(x = year, y = median)) + 
        geom_errorbar(aes(ymin = qlo, ymax = qhi), size = 0.4) +
        geom_point(size = 0.4) +
        labs(x = "Year", y = bquote(phi), title = "R.E. AR(1)") +
        fivepointtwo()

gg23 <- plot_grid(gg2, gg3, ncol = 1)

tiff(filename = "Paper/Figures/Fig7.tif", width = 5.2, height = 2.5, units = "in", res = 600)
plot_grid(gg1, gg23, nrow = 1, rel_widths = c(3,1))
dev.off()

pdf(file = "Paper/Joint_Count_Files/sigestl-1.pdf", width = 5.2, height = 2.5)
plot_grid(gg1, gg23, nrow = 1, rel_widths = c(3,1))
dev.off()









# Fig8 -------------------------------------------------------------------------
#```{r sigestp, fig.cap="\\label{sigestp}Variance parameters for person-caused fires. Again, there is a slight decrease in the estimates for region 5 and $\\sigma_b$ with an increase in the AR(1) parameter."}
# gamma
gg1 <- perest %>% 
    select(year, starts_with("sigma_x")) %>%
    group_by(year) %>%
    summarise_all(.funs = quantile, prob = c(0.025, 0.5, 0.975)) %>%
    mutate(fun = c("qlo", "median", "qhi")) %>%
    tidyr::pivot_longer(cols = 2:7, names_to = "Region") %>%
    mutate(Region = gsub("sigma_x_", "", Region)) %>%
    tidyr::pivot_wider(names_from = "fun", values_from = value) %>% 
    ggplot(aes(x = as.numeric(year), y = median)) + 
        geom_errorbar(mapping = aes(ymin = qlo, ymax = qhi), size = 0.4) +
        geom_point(size = 0.4) + #geom_smooth(method = "loess", se = FALSE) + 
        facet_wrap(~Region, nrow = 2, 
            labeller = labeller(Region = label_both)) +
        theme_bw() +
        labs(x = "year", y = bquote(gamma[r]), 
            title = "LogNormal Variance - Person Caused Fires",
            colour = '"Significance"') +
        geom_hline(aes(yintercept = 0), colour = rgb(0,0,0,0)) +
        theme(legend.position = "none", 
            axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
        scale_x_continuous(breaks = seq(1950,2015,10)) +
        fivepointtwo()

gg2 <- perest %>% 
    select(year, sigma_b) %>%
    group_by(year) %>%
    summarise(qlo = quantile(sigma_b, 0.025),
        median = median(sigma_b),
        qhi = quantile(sigma_b, 0.975)) %>%
    ggplot(aes(x = year, y = median)) + 
        geom_errorbar(aes(ymin = qlo, ymax = qhi), size = 0.4) +
        geom_point(size = 0.4) +
        labs(x = "Year", y = bquote(sigma[b]), title = "R.E. Variance") +
        fivepointtwo()

gg3 <- perest %>% 
    select(year, phi) %>%
    group_by(year) %>%
    summarise(qlo = quantile(phi, 0.025),
        median = median(phi),
        qhi = quantile(phi, 0.975)) %>%
    ggplot(aes(x = year, y = median)) + 
        geom_errorbar(aes(ymin = qlo, ymax = qhi), size = 0.4) +
        geom_point(size = 0.4) +
        labs(x = "Year", y = bquote(phi), title = "R.E. AR(1)") +
        fivepointtwo()

gg23 <- plot_grid(gg2, gg3, ncol = 1)

tiff(filename = "Paper/Figures/Fig8.tif", width = 5.2, height = 2.5, units = "in", res = 600)
plot_grid(gg1, gg23, nrow = 1, rel_widths = c(3,1))
dev.off()

pdf(file = "Paper/Joint_Count_Files/sigestp-1.pdf", width = 5.2, height = 2.5)
plot_grid(gg1, gg23, nrow = 1, rel_widths = c(3,1))
dev.off()








# Fig9 -------------------------------------------------------------------------

#```{r powers, fig.height=8.5, fig.cap = "\\label{powers}Results of the power study. The credible interval test and the Bayes factors that relied on posterior samples performed the best overall. Simply looking at whether the WAIC is smaller has high power at the expense of high type 1 error. The Savage Dickey representation of the BF might have higher power with better prior distributions. Using the variable selection method on $\\gamma$ has very low power, likely because $\\gamma$ is the coefficient of a latent variable rather than an observed covariate. The equivalence test and probability of direction have slightly lower power than simply checking the credible interval."}

pows <- bind_rows(
    readRDS("Data/Models/2019-10-02 SimiItAgain_part1.RDS"),
    readRDS("Data/Models/2019-10-02 SimiItAgain_part2.RDS"),
    readRDS("Data/Models/2019-10-02 SimiItAgain_part3.RDS"))
naive <- readRDS("Data/Models/2019-10-02 SimiItAgainNaive_part1.RDS")

pow1 <- filter(pows, ndays == 50)
pow2 <- filter(pows, ndays != 50)

pow1$NB <- NA # No Naive bayes BFs for ndays == 50
pow2$NB <- naive$BF

pow <- bind_rows(pow1, pow2)

sdb_labs <- c(bquote(sigma[b]*"="*0.25), bquote(sigma[b]*"="*0.5), bquote(sigma[b]*"="*1))
names(sdb_labs) <- c(0.25,0.5,1)

pr_labs <- c(bquote(p[r]*"="*0.25), bquote(p[r]*"="*0.75))
names(pr_labs) <- c(0.25,0.75)

lam_labs <- c(bquote(lambda*"=-1"), bquote(lambda*"=-0.5"), bquote(lambda*"=0"))
names(lam_labs) <- c(-1,-0.5,0)

mux_labs <- c(bquote(mu*"=-2"), bquote(mu*"=-1"))
names(mux_labs) <- c(-2,-1)

#pows$BF <- naive$BF

tencols <- c('#a6cee3','#1f78b4','#b2df8a',
    '#33a02c','#fb9a99', '#e31a1c','#fdbf6f',
    '#ff7f00','#cab2d6','#6a3d9a')



tiff(filename = "Paper/Figures/Fig9-1.tif", width = 5.2, height = 3.5, units = "in", res = 600)
pow %>% 
    filter((!is.na(CI)) | (!is.na(NB))) %>% 
    filter(ndays > 50) %>% 
    group_by(gamma, sdb, ndays, lam, mux, roundprob) %>% 
    summarise(Harmonic_Mean = mean(HM > 3, na.rm = TRUE),
              Savage_Dickey = mean(SD > 3, na.rm = TRUE),
              Credible_Interval = mean(CI > 0, na.rm = TRUE),
              Raw_Mean = mean(RM > 3, na.rm = TRUE),
              #BIC = mean(BI, na.rm = TRUE),
              WAIC = mean(WA, na.rm = TRUE),
              Equiv = mean(EQ == "rejected", na.rm = TRUE),
              pd = mean(PD > 95, na.rm = TRUE),
              GS33 = mean(1/(1-GS) > 3, na.rm = TRUE),
              Naive = mean(NB > 3, na.rm = TRUE)) %>% 
    ungroup() %>% 
    gather(test, power, 
           -gamma, -sdb, -ndays, -lam, -mux, -roundprob) %>%
    filter(sdb == 0.25) %>%
    mutate(sdb = factor(sdb, labels = sdb_labs[1]),
           lam = factor(lam, labels = lam_labs),
           mux = factor(mux, labels = mux_labs),
           roundprob = factor(roundprob, labels = pr_labs)) %>% 
    mutate(test = case_when(
        test == "BIC" ~ "BIC[d] - BIC[i] < 0",
        test == "Equiv" ~ "Equiv. Test",
        test == "Harmonic_Mean" ~ "Harmonic Mean",
        test == "pd" ~ "Prob. of Direction",
        test == "Savage_Dickey" ~ "Savage Dickey",
        test == "Credible_Interval" ~ "Credible Interval",
        test == "GS33" ~ "gamma Selector",
        test == "Naive" ~ "Naive Samples BF",
        test == "Raw_Mean" ~ "Post. Samples BF",
        test == "WAIC" ~ "WAIC[d] - WAIC[i] < 0"
    )) %>% 
    ggplot(aes(x = gamma, colour = test, 
               shape = test, y = power)) + 
    theme_bw() + 
    geom_line(size = 0.3) + geom_point(size = 0.3) + 
    facet_grid(sdb + lam ~ mux + roundprob, 
               labeller = label_parsed) +
    scale_x_continuous(breaks = seq(0,1,0.2), 
                       minor_breaks = seq(0,1,0.1)) +
    scale_y_continuous(breaks = seq(0,1,0.2)) +
    scale_shape_manual(values = 1:10) +
    coord_cartesian(xlim = c(0,1), ylim = c(0,1)) + 
    scale_colour_manual(values = tencols) + 
    theme(legend.position = "bottom") +
    labs(x = bquote(gamma), y = "Power", 
         colour = "Test", shape = "Test") +
    fivepointtwo()

dev.off()

tiff(filename = "Paper/Figures/Fig9-2.tif", width = 5.2, height = 3.5, units = "in", res = 600)
pow %>% 
    filter((!is.na(CI)) | (!is.na(NB))) %>% 
    filter(ndays > 50) %>% 
    group_by(gamma, sdb, ndays, lam, mux, roundprob) %>% 
    summarise(Harmonic_Mean = mean(HM > 3, na.rm = TRUE),
              Savage_Dickey = mean(SD > 3, na.rm = TRUE),
              Credible_Interval = mean(CI > 0, na.rm = TRUE),
              Raw_Mean = mean(RM > 3, na.rm = TRUE),
              #BIC = mean(BI, na.rm = TRUE),
              WAIC = mean(WA, na.rm = TRUE),
              Equiv = mean(EQ == "rejected", na.rm = TRUE),
              pd = mean(PD > 95, na.rm = TRUE),
              GS33 = mean(1/(1-GS) > 3, na.rm = TRUE),
              Naive = mean(NB > 3, na.rm = TRUE)) %>% 
    ungroup() %>% 
    gather(test, power, 
           -gamma, -sdb, -ndays, -lam, -mux, -roundprob) %>%
    filter(sdb == 0.5) %>%
    mutate(sdb = factor(sdb, labels = sdb_labs[2]),
           lam = factor(lam, labels = lam_labs),
           mux = factor(mux, labels = mux_labs),
           roundprob = factor(roundprob, labels = pr_labs)) %>% 
    mutate(test = case_when(
        test == "BIC" ~ "BIC[d] - BIC[i] < 0",
        test == "Equiv" ~ "Equiv. Test",
        test == "Harmonic_Mean" ~ "Harmonic Mean",
        test == "pd" ~ "Prob. of Direction",
        test == "Savage_Dickey" ~ "Savage Dickey",
        test == "Credible_Interval" ~ "Credible Interval",
        test == "GS33" ~ "gamma Selector",
        test == "Naive" ~ "Naive Samples BF",
        test == "Raw_Mean" ~ "Post. Samples BF",
        test == "WAIC" ~ "WAIC[d] - WAIC[i] < 0"
    )) %>% 
    ggplot(aes(x = gamma, colour = test, 
               shape = test, y = power)) + 
    theme_bw() + 
    geom_line(size = 0.3) + geom_point(size = 0.3) + 
    facet_grid(sdb + lam ~ mux + roundprob, 
               labeller = label_parsed) +
    scale_x_continuous(breaks = seq(0,1,0.2), 
                       minor_breaks = seq(0,1,0.1)) +
    scale_y_continuous(breaks = seq(0,1,0.2)) +
    scale_shape_manual(values = 1:10) +
    coord_cartesian(xlim = c(0,1), ylim = c(0,1)) + 
    scale_colour_manual(values = tencols) + 
    theme(legend.position = "bottom") +
    labs(x = bquote(gamma), y = "Power", 
         colour = "Test", shape = "Test") +
    fivepointtwo()
dev.off()

tiff(filename = "Paper/Figures/Fig9-3.tif", width = 5.2, height = 3.5, units = "in", res = 600)
pow %>% 
    filter((!is.na(CI)) | (!is.na(NB))) %>% 
    filter(ndays > 50) %>% 
    group_by(gamma, sdb, ndays, lam, mux, roundprob) %>% 
    summarise(Harmonic_Mean = mean(HM > 3, na.rm = TRUE),
              Savage_Dickey = mean(SD > 3, na.rm = TRUE),
              Credible_Interval = mean(CI > 0, na.rm = TRUE),
              Raw_Mean = mean(RM > 3, na.rm = TRUE),
              #BIC = mean(BI, na.rm = TRUE),
              WAIC = mean(WA, na.rm = TRUE),
              Equiv = mean(EQ == "rejected", na.rm = TRUE),
              pd = mean(PD > 95, na.rm = TRUE),
              GS33 = mean(1/(1-GS) > 3, na.rm = TRUE),
              Naive = mean(NB > 3, na.rm = TRUE)) %>% 
    ungroup() %>% 
    gather(test, power, 
           -gamma, -sdb, -ndays, -lam, -mux, -roundprob) %>%
    filter(sdb == 1) %>%
    mutate(sdb = factor(sdb, labels = sdb_labs[3]),
           lam = factor(lam, labels = lam_labs),
           mux = factor(mux, labels = mux_labs),
           roundprob = factor(roundprob, labels = pr_labs)) %>% 
    mutate(test = case_when(
        test == "BIC" ~ "BIC[d] - BIC[i] < 0",
        test == "Equiv" ~ "Equiv. Test",
        test == "Harmonic_Mean" ~ "Harmonic Mean",
        test == "pd" ~ "Prob. of Direction",
        test == "Savage_Dickey" ~ "Savage Dickey",
        test == "Credible_Interval" ~ "Credible Interval",
        test == "GS33" ~ "gamma Selector",
        test == "Naive" ~ "Naive Samples BF",
        test == "Raw_Mean" ~ "Post. Samples BF",
        test == "WAIC" ~ "WAIC[d] - WAIC[i] < 0"
    )) %>% 
    ggplot(aes(x = gamma, colour = test, 
               shape = test, y = power)) + 
    theme_bw() + 
    geom_line(size = 0.3) + geom_point(size = 0.3) + 
    facet_grid(sdb + lam ~ mux + roundprob, 
               labeller = label_parsed) +
    scale_x_continuous(breaks = seq(0,1,0.2), 
                       minor_breaks = seq(0,1,0.1)) +
    scale_y_continuous(breaks = seq(0,1,0.2)) +
    scale_shape_manual(values = 1:10) +
    coord_cartesian(xlim = c(0,1), ylim = c(0,1)) + 
    scale_colour_manual(values = tencols) + 
    theme(legend.position = "bottom") +
    labs(x = bquote(gamma), y = "Power", 
         colour = "Test", shape = "Test") +
    fivepointtwo()
dev.off()
