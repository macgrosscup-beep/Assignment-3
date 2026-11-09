# Assignment 3, Part 1: model modification and selection
 

library(bbmle)



library(biometrydata)


# 1a: Competition between bean beetles ----------------------------------------

beetles <- beetle_competition_interspecific

beetles <- beetles[beetles$pair == "Cmac_Csub" & !is.na(beetles$progeny_1), ]

plot(beetles$start_females_1 + beetles$start_females_2,
     beetles$progeny_1 / beetles$start_females_1,
     xlab = "females in the dish, both species",
     ylab = "C. maculatus offspring per female")

# Model 1: starting model (Beverton-Holt family with d)
beetle.fit <- mle2(progeny_1 ~ dnbinom(mu = R * start_females_1 /
                                         (1 + (a11 * start_females_1 +
                                                 a12 * start_females_2)^d),
                                       size = k),
                   start = list(R = 20, a11 = 0.1, a12 = 0.1, d = 1, k = 2),
                   lower = c(R = 0.01, a11 = 0, a12 = 0, d = 0.1, k = 0.01),
                   upper = c(R = 1000, a11 = 10, a12 = 10, d = 10, k = 100),
                   method = "L-BFGS-B",
                   data = beetles)

# Model 2: classic Beverton-Holt (d = 1)
beetle.bh <- mle2(progeny_1 ~ dnbinom(mu = R * start_females_1 /
                                        (1 + a11 * start_females_1 +
                                           a12 * start_females_2),
                                      size = k),
                  start = list(R = 20, a11 = 0.1, a12 = 0.1, k = 2),
                  lower = c(R = 0.01, a11 = 0, a12 = 0, k = 0.01),
                  upper = c(R = 1000, a11 = 10, a12 = 10, k = 100),
                  method = "L-BFGS-B",
                  data = beetles)

# Model 3: Beverton-Holt, no effect of C. subinnotatus (a12 = 0)
beetle.no_a12 <- mle2(progeny_1 ~ dnbinom(mu = R * start_females_1 /
                                            (1 + a11 * start_females_1),
                                          size = k),
                      start = list(R = 20, a11 = 0.1, k = 2),
                      lower = c(R = 0.01, a11 = 0, k = 0.01),
                      upper = c(R = 1000, a11 = 10, k = 100),
                      method = "L-BFGS-B",
                      data = beetles)

# Model 4: Beverton-Holt, both species equally harmful (one parameter, a)
beetle.equal <- mle2(progeny_1 ~ dnbinom(mu = R * start_females_1 /
                                           (1 + a * (start_females_1 +
                                                       start_females_2)),
                                         size = k),
                     start = list(R = 20, a = 0.1, k = 2),
                     lower = c(R = 0.01, a = 0, k = 0.01),
                     upper = c(R = 1000, a = 10, k = 100),
                     method = "L-BFGS-B",
                     data = beetles)

# Model 5: Ricker shape
beetle.ricker <- mle2(progeny_1 ~ dnbinom(mu = R * start_females_1 *
                                            exp(-a11 * start_females_1 -
                                                  a12 * start_females_2),
                                          size = k),
                      start = list(R = 20, a11 = 0.01, a12 = 0.01, k = 2),
                      lower = c(R = 0.01, a11 = 0, a12 = 0, k = 0.01),
                      upper = c(R = 1000, a11 = 10, a12 = 10, k = 100),
                      method = "L-BFGS-B",
                      data = beetles)

AICtab(beetle.fit, beetle.bh, beetle.no_a12, beetle.equal, beetle.ricker,
       base = TRUE, weights = TRUE)

coef(beetle.bh)
coef(beetle.bh)["a12"] / coef(beetle.bh)["a11"]

# Conclusion:
#   - d isn't needed. The classic Beverton-Holt (d = 1) has a lower AIC than
#     the starting model (~1 unit), so the extra parameter isn't supported.
#   - C. subinnotatus does affect C. maculatus. Dropping a12 makes the model
#     ~3.3 AIC units worse.
#   - But it isn't as harmful as C. maculatus itself. Forcing one shared
#     competition parameter is ~3.6 units worse, and a12 is only about a
#     third of a11. Competition within species is stronger than competition
#     between species, which is what you'd expect if the two can coexist.
#   - The Ricker isn't nested in the starting model (no value of d turns one
#     into the other), but AIC can still compare them since both use the
#     same data and response. It's ~1 unit better than Beverton-Holt, so the
#     data can't really tell the two shapes apart.
#   - Best model: Beverton-Holt or Ricker, with separate a11 and a12.


# 1b: Wing length in hawks and eagles -----------------------------------------

hawks <- avonet


hawks <- hawks[hawks$Family1 == "Accipitridae", c("Species1", "Wing.Length", "Mass", "Migration")]
hawks <- na.omit(hawks)
hawks$migration <- factor(hawks$Migration, levels = 1:3,
                          labels = c("resident", "partial", "migrant"))

# new column: resident vs. migratory (partial and full migrants together)
hawks$migratory <- factor(hawks$migration == "resident",
                          levels = c(TRUE, FALSE),
                          labels = c("resident", "migratory"))

plot(log(hawks$Mass), log(hawks$Wing.Length), col = hawks$migration,
     xlab = "log body mass", ylab = "log wing length")
legend("topleft", levels(hawks$migration), col = 1:3, pch = 1)

# Model 1: starting model (intercept, slope and scatter all differ)
hawk.fit <- mle2(Wing.Length ~ dlnorm(meanlog = int + b * log(Mass),
                                      sdlog = exp(log_s)),
                 start = list(int = 3, b = 0.3, log_s = -2),
                 parameters = list(int ~ migration,
                                   b ~ migration,
                                   log_s ~ migration),
                 data = hawks)

# Model 2: only the intercept differs among the 3 classes
hawk.int <- mle2(Wing.Length ~ dlnorm(meanlog = int + b * log(Mass),
                                      sdlog = exp(log_s)),
                 start = list(int = 3, b = 0.3, log_s = -2),
                 parameters = list(int ~ migration),
                 data = hawks)

# Model 3: only the intercept differs, resident vs. migratory
hawk.int2 <- mle2(Wing.Length ~ dlnorm(meanlog = int + b * log(Mass),
                                       sdlog = exp(log_s)),
                  start = list(int = 3, b = 0.3, log_s = -2),
                  parameters = list(int ~ migratory),
                  data = hawks)

# Model 4: no migration effect
hawk.null <- mle2(Wing.Length ~ dlnorm(meanlog = int + b * log(Mass),
                                       sdlog = exp(log_s)),
                  start = list(int = 3, b = 0.3, log_s = -2),
                  data = hawks)

AICtab(hawk.fit, hawk.int, hawk.int2, hawk.null,
       base = TRUE, weights = TRUE)

coef(hawk.int2)
exp(coef(hawk.int2)["int.migratorymigratory"])

# Conclusion:
#   - Migrants do have longer wings for their mass. That's the intercept,
#     and the model with no migration effect is ~14 AIC units worse than
#     the best model.
#   - Slope and scatter don't need to differ. The starting model, which lets
#     all three vary, is ~7.5 units worse than letting only the intercept
#     vary. Wing length scales with mass the same way in every group
#     (b ~ 0.33), and no group is more variable.
#   - Partial vs. full migrants isn't needed. The 2-group model is slightly
#     better than the 3-group model (~0.8 units), so the simpler split
#     is enough.
#   - Best model: shared slope and scatter, with migratory species having
#     ~8% longer wings than residents of the same mass. That fits the idea
#     that longer wings make sustained flight cheaper.
#   - Caveat: species aren't fully independent because related species
#     share traits, so the evidence is probably a bit weaker than AIC
#     suggests.


# 1c: Reef fish species richness ----------------------------------------------

reef <- reef_fish_surveys

reef$depth <- reef$DEPTH / 10
reef$visibility <- reef$UNDERWATER_VISIBILITY / 10
reef$year <- reef$YEAR - 2014
reef$region <- factor(reef$REGION_DESCRIPTION,
                      levels = c("Dry Tortugas",
                                 "Florida Keys: from Key West to Miami",
                                 "Miami to North Martin County"),
                      labels = c("Tortugas", "Keys", "North"))
reef$protected <- ifelse(reef$MPA_NAME == "Unprotected", "no", "yes")

plot(reef$DEPTH, reef$species_richness, col = rgb(0, 0, 0, 0.1),
     xlab = "depth (m)", ylab = "species recorded")

# Check the Poisson: variance / mean should be about 1
var(reef$species_richness) / mean(reef$species_richness)

# Model 1: starting model (Poisson)
reef.fit <- mle2(species_richness ~ dpois(lambda = exp(int +
                                                         d1 * depth +
                                                         d2 * depth^2 +
                                                         v * visibility +
                                                         t * year)),
                 start = list(int = 3, d1 = 0, d2 = 0, v = 0, t = 0),
                 parameters = list(int ~ region + protected),
                 data = reef)

# Model 2: same terms, negative binomial
reef.nb <- mle2(species_richness ~ dnbinom(mu = exp(int +
                                                      d1 * depth +
                                                      d2 * depth^2 +
                                                      v * visibility +
                                                      t * year),
                                           size = k),
                start = list(int = 3, d1 = 0, d2 = 0, v = 0, t = 0, k = 5),
                parameters = list(int ~ region + protected),
                data = reef)

# Model 3: negative binomial, no protection
reef.noprot <- mle2(species_richness ~ dnbinom(mu = exp(int +
                                                          d1 * depth +
                                                          d2 * depth^2 +
                                                          v * visibility +
                                                          t * year),
                                               size = k),
                    start = list(int = 3, d1 = 0, d2 = 0, v = 0, t = 0, k = 5),
                    parameters = list(int ~ region),
                    data = reef)

# Model 4: negative binomial, no visibility
reef.novis <- mle2(species_richness ~ dnbinom(mu = exp(int +
                                                         d1 * depth +
                                                         d2 * depth^2 +
                                                         t * year),
                                              size = k),
                   start = list(int = 3, d1 = 0, d2 = 0, t = 0, k = 5),
                   parameters = list(int ~ region + protected),
                   data = reef)

AICtab(reef.fit, reef.nb, reef.noprot, reef.novis,
       base = TRUE, weights = TRUE)


round(exp(coef(reef.nb)), 3)

# Conclusion:
#   - The Poisson is wrong. Variance is ~4.4 times the mean, and the negative
#     binomial lowers AIC by about 13,500. Under the Poisson the standard
#     errors are too small, so every effect looks stronger than it is.
#   - With 8,560 surveys, AIC keeps every term, so the real question is how
#     big each effect is (exp of the coefficients):
#       region: Keys and North have ~35-38% more species than Tortugas
#       depth: hump-shaped, with richness peaking at ~17 m
#       visibility: ~14% more species per extra 10 m of visibility
#       protection: ~7% more species inside MPAs
#       year: under 1% per year
#   - Visibility describes the water on the day, not the reef. Its effect is
#     most likely detection: in clear water, divers can see and record more
#     of the species that are there. So it's worth keeping as a correction
#     for survey conditions, not as a cause of diversity.
#   - Protection still matters once region is in the model (dropping it is
#     ~36 units worse), but it's a modest effect. 
#   - Best model: the negative binomial with all terms. Region and depth
#     drive most of the differences; protection and year are small effects,
#     and visibility mostly reflects detection.


# Assignment 3, Part 2: focal dataset model selection

#
# Data: my coastal plant field data (Masters_Data_Modified.xlsx). 31 plots along a soil
# salinity gradient, and individual plants measured for traits. 
#
# Two response variables:
#   1. Species richness per plot (count). Does richness drop as soil
#      salinity goes up?
#   2. Leaf thickness of individual plants (mm, continuous and positive).
#      Do plants have thicker (more succulent) leaves in saltier soil, and is
#      that just because different species grow there, or do they also
#      change within a species?

library(bbmle)
library(readxl)

plants <- as.data.frame(read_excel("Masters_Data_Modified.xlsx"))

# Fixed some typos; whoops :)
plants$species <- trimws(plants$species)
plants$species <- gsub("Sporobulus", "Sporobolus", plants$species)
plants$species <- gsub("Langularia", "Laguncularia", plants$species)


# Response 1: species richness per plot -----------------------------------------

# One row per plot: salinity and number of species
plots <- aggregate(soil_salinity_ms ~ site + plot, data = plants, FUN = mean)
rich <- aggregate(species ~ site + plot, data = plants,
                  FUN = function(x) length(unique(x)))
names(rich)[3] <- "richness"
plots <- merge(plots, rich, all.x = TRUE)
plots$richness[is.na(plots$richness)] <- 0   # fds plot 2 had no plants

plot(plots$soil_salinity_ms, plots$richness, col = factor(plots$site),
     xlab = "soil salinity (mS)", ylab = "species per plot")
legend("topright", levels(factor(plots$site)), col = 1:5, pch = 1)

var(plots$richness) / mean(plots$richness)

# Model 1: no effect (richness is the same everywhere)
rich.null <- mle2(richness ~ dpois(lambda = exp(int)),
                  start = list(int = 1),
                  data = plots)

# Model 2: richness changes with salinity
rich.sal <- mle2(richness ~ dpois(lambda = exp(int + b * soil_salinity_ms)),
                 start = list(int = 1, b = 0),
                 data = plots)

# Model 3: richness differs by site only
rich.site <- mle2(richness ~ dpois(lambda = exp(int)),
                  start = list(int = 1),
                  parameters = list(int ~ site),
                  data = plots)

# Model 4: salinity effect plus site differences
rich.sal_site <- mle2(richness ~ dpois(lambda = exp(int + b * soil_salinity_ms)),
                      start = list(int = 1, b = 0),
                      parameters = list(int ~ site),
                      data = plots)

AICtab(rich.null, rich.sal, rich.site, rich.sal_site,
       base = TRUE, weights = TRUE)

coef(rich.sal_site)
exp(coef(rich.sal_site)["b"])         # multiplier per +1 mS
exp(4 * coef(rich.sal_site)["b"])     # freshest (~0 mS) vs saltiest (4 mS) plots

# Conclusion:
#   - Salinity and site both matter. The site + salinity model is best by a
#     wide margin over salinity alone, site alone, or neither.
#   - Within a site, each extra 1 mS of soil salinity means ~19% fewer
#     species, so the saltiest plots (4 mS) have about 55% fewer species
#     than the freshest ones. That fits salt stress filtering out species
#     that can't tolerate it.
#   - Site matters mainly because eps is much richer than the others
#     (~2.7x at the same salinity), which makes sense given it's the
#     largest site with the most habitat variety.
#   - Raw richness looks overdispersed, but most of that
#     comes from the site differences, so once site is in the model the
#     Poisson is reasonable.



# Response 2: leaf thickness of individual plants -------------------------------

leaves <- plants[!is.na(plants$leaf_thickness_mm), ]

plot(leaves$soil_salinity_ms, log(leaves$leaf_thickness_mm),
     col = rgb(0, 0, 0, 0.2),
     xlab = "soil salinity (mS)", ylab = "log leaf thickness (mm)")

# Model 1: salinity only (ignores which species the plant is)
leaf.sal <- mle2(leaf_thickness_mm ~ dlnorm(meanlog = int + b * soil_salinity_ms,
                                            sdlog = exp(log_s)),
                 start = list(int = 0, b = 0, log_s = 0),
                 data = leaves)

# Model 2: species only (each species has its own thickness)
leaf.sp <- mle2(leaf_thickness_mm ~ dlnorm(meanlog = int,
                                           sdlog = exp(log_s)),
                start = list(int = 0, log_s = 0),
                parameters = list(int ~ species),
                data = leaves)

# Model 3: species plus a salinity effect within species
leaf.sp_sal <- mle2(leaf_thickness_mm ~ dlnorm(meanlog = int + b * soil_salinity_ms,
                                               sdlog = exp(log_s)),
                    start = list(int = 0, b = 0, log_s = 0),
                    parameters = list(int ~ species),
                    data = leaves)

# Model 4: same as model 2, gamma instead of lognormal
leaf.sp_gam <- mle2(leaf_thickness_mm ~ dgamma(shape = exp(log_sh),
                                               rate = exp(log_sh) / exp(int)),
                    start = list(int = 0, log_sh = 0),
                    parameters = list(int ~ species),
                    data = leaves)

# Model 5: same as model 3, gamma instead of lognormal
leaf.sp_sal_gam <- mle2(leaf_thickness_mm ~ dgamma(shape = exp(log_sh),
                                                   rate = exp(log_sh) /
                                                     exp(int + b * soil_salinity_ms)),
                        start = list(int = 0, b = 0, log_sh = 0),
                        parameters = list(int ~ species),
                        data = leaves)

AICtab(leaf.sal, leaf.sp, leaf.sp_sal, leaf.sp_gam, leaf.sp_sal_gam,
       base = TRUE, weights = TRUE)

exp(coef(leaf.sal)["b"])                    # per +1 mS, ignoring species
exp(coef(leaf.sp_sal_gam)["b"])             # per +1 mS, within species
exp(4 * coef(leaf.sp_sal_gam)["b"])         # 0 vs 4 mS, within species

# Conclusion:
#   - Species identity explains most of the variation in leaf thickness.
#     Salinity alone is by far the worst model.
#   - Salinity still matters once species is accounted for. Adding it
#     improves the species model under both distributions.
#   - Ignoring species, leaves look ~20% thicker per 1 mS, but within a
#     species the effect is only ~6% per 1 mS (~26% thicker from 0 to
#     4 mS). So most of the overall pattern comes from turnover in which
#     species grow where (succulents like Sesuvium and Suaeda in salty
#     plots), with a smaller real within-species shift toward thicker,
#     more succulent leaves in saltier soil.
#   - Distribution: the gamma fits much better than the lognormal with the
#     same mean structure, and the salinity result is similar under both.
#   - Best model: gamma, species-specific thickness plus a shared salinity
#     effect.