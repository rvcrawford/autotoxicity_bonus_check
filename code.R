library(data.table)
library(tidyverse)

dat <- fread("./input_data/selectomatic AA results(seedlot validation test).csv")

n_entries <- dat[,.N, by = entry]
in_all <- n_entries[N>20]

in_all
dat

# just look at rightside up
has_data <- dat[!is.na(rightsideup)]
has_data[,.N, by = entry]

dat$tray |> unique()

dat

library(ordinal)

dat$score <- ordered(dat$score)
dat$entry <- factor(dat$entry)  # Treat as categorical
dat$tray <- factor(dat$tray, levels = c("weakest", "weak", "strong", "strongest"))

# Entry effects with tray as block
model1 <- clm(score ~ entry + tray, data = dat)
summary(model1)

# Get entry comparisons
library(emmeans)
emmeans(model1, pairwise ~ entry)


# set reference level as entry 25
dat$entry <- relevel(factor(dat$entry), ref = "25")
dat2 <- dat[tray=="weakest"|tray=="strongest",]
model_check <- clm(score ~ entry + tray, data = dat2)
summary(model_check)

library(emmeans)

# Get estimated marginal means for all entries
em <- emmeans(model_check, ~ entry)



em |> 
  data.frame() |> 
  ggplot() + 
  aes(x = fct_reorder(entry, emmean, .desc = T), emmean) + 
  geom_point() + 
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL))+
  theme_bw()


contrast(em, "trt.vs.ctrl", ref = c(1,2,3))


has_data

for_rightsideup <- dat |> filter(as.numeric(score)>2) |> 
  mutate(rightsideup = replace_na(rightsideup,1))

library(brms)

for_rightsideup$entry <- factor(for_rightsideup$entry)

bayes_model <- brm(
  rightsideup ~ entry + tray,
  data = for_rightsideup,
  family = bernoulli(link = "logit"),
  prior = c(
    prior(normal(0, 2), class = "b"),
    prior(normal(0, 5), class = "Intercept")
  ),
  chains = 4,
  iter = 2000,
  warmup = 1000,
  cores = 4
)

# check priors
prior_check <- brm(
  rightsideup ~ entry + tray,
  data = dat_binary,
  family = bernoulli(link = "logit"),
  prior = c(
    prior(normal(0, 2), class = "b"),
    prior(normal(0, 5), class = "Intercept")
  ),
  sample_prior = "only",  # Only sample from prior, ignore data
  chains = 4,
  iter = 2000
)

# Look at predicted probabilities from the prior
pp_check(prior_check, type = "bars", ndraws = 50)

post_samples <- as_draws_df(bayes_model) |> select(b_Intercept,contains("entry"))

inv_logit <- function(x) exp(x) / (1 + exp(x))

answers <- apply(post_samples[,2:79],2,
    function(x)quantile(post_samples$b_Intercept+x, probs = c(0.025, 0.5, 0.975)) |> inv_logit()) |> as_tibble() |> 
  mutate(prob = c("pct2.5","pct50", "pct97.5"), .before = 1) |> 
  pivot_longer(-1) |> 
  pivot_wider(names_from = prob, values_from = value) |> 
  rename(variable = name) |> 
  mutate(variable = parse_number(variable) |> as.character())
  
answers |> 
  ggplot() + 
  aes(x = fct_reorder(variable, pct50, .desc = T), y = pct50) + 
  geom_point() + 
  geom_errorbar(aes(ymin = pct2.5, ymax = pct97.5)) + 
  theme_bw()
