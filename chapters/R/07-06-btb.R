# Chapter 7, Section 5
# Spatio-temporal analysis of Bovine Tubercolosis in Cornwall
source("00-prelim.R")

## ---- data.btb ----
load("data/BTBppp.RData")
plot.df <- as.data.frame(pppdata$window)
W <- pppdata$window
simpW <- simplify.owin(W, 1000)  # This reduces resolution of plot
plot.df <- as.data.frame(W)

# Get data points (location, spoligotypes and year)
btb <- read.table("data/BTB_spoligotype_data.txt", header = FALSE, sep = " ",
                  skip = 1)
colnames(btb) <- c("x", "y", "year", "sp")
btb$sp <- factor(btb$sp)

# Keep the four largest classes
btb %>%
  group_by(sp) %>%
  tally() %>%
  arrange(desc(n)) -> btb.n
(classes.to.keep <- as.character(btb.n$sp[1:4]))
levels(btb$sp) <- c("Sp9",    "Others", "Others", "Sp12" ,  "Sp15",
                    "Others", "Sp20",   "Others", "Others", "Others")
btb$sp <- factor(btb$sp, levels = c("Sp9", "Sp12", "Sp15", "Sp20", "Others"))

# A similar table to Diggle et al. (2005)
btb.summary <- btb %>%
  group_by(year, sp) %>%
  arrange(year) %>%
  tally() %>%
  complete(sp, year, fill = list(n = 0)) %>%
  dcast(year ~ sp)

# Group the years into 4 categories: 1) < 1997; 2) 1997-1998; 3) 1999-2000; 4)
# 2001-02
btb_year_fn <- function(x) {
  res <- x
  res[x < 1997] <- "< 1997"
  res[x == 1997 | x == 1998] <- "1997-1998"
  res[x == 1999 | x == 2000] <- "1999-2000"
  res[x == 2001 | x == 2002] <- "2001-2002"
  res
}
as.tibble(btb) %>%
  mutate(period = factor(btb_year_fn(year))) -> btb

## ---- plot.cow ----
as.tibble(btb.summary) %>%
  melt(id.vars = "year") %>%
  ggplot(aes(year, value, col = variable)) +
  geom_segment(aes(xend = year, yend = 0), size = 12) +
  scale_x_continuous(breaks = seq(1989, 2002, by = 1),
                     minor_breaks = seq(1989, 2002, by = 1)) +
  labs(x = "Year", y = "Count") +
  scale_color_discrete(name = "Spoligotype") +
  theme_bw() +
  guides(col = guide_legend(override.aes = list(size = 6)))

## ---- plot.cornwall ----
ggplot() +
  geom_polygon(data = plot.df, aes(x, y, group = id), fill = NA, col = "grey25") +
  geom_point(data = btb, aes(x, y, col = sp)) +
  labs(x = "Eastings (1,000 km)", y = "Northings (1,000 km)",
       col = "Spoligotype") +
  scale_x_continuous(labels = function(x) x / 1000) +
  scale_y_continuous(labels = function(x) x / 1000) +
  theme_bw() +
  theme(legend.position = c(0.98, 0.005), legend.justification = c(1, 0))

## ---- mod.btb ----
Spoligotype <- btb$sp
X <- scale(btb[, 1:2])
mu.x <- attr(X, "scaled:center")
sd.x <- attr(X, "scaled:scale")
year <- scale(btb$year, scale = FALSE)
mu.year <- attr(year, "scaled:center")
period <- btb$period

# # Spatial model only
# mod1 <- iprobit(y = Spoligotype, X, kernel = "FBM",
#                 control = list(restarts = 8, maxit = 500,
#                                restart.method = "error",
#                                common.RKHS.scale = TRUE,
#                                common.intercept = FALSE))
#
# mod1a <- iprobit(y = Spoligotype, X, kernel = "FBM",
#                  control = list(restarts = 8, maxit = 500,
#                                 restart.method = "error",
#                                common.RKHS.scale = FALSE,
#                                common.intercept = FALSE))
#
# # Temporal only
# mod2 <- iprobit(y = Spoligotype, year, kernel = "FBM",
#                 control = list(restarts = 8, maxit = 500,
#                                restart.method = "error",
#                                common.RKHS.scale = TRUE,
#                                common.intercept = FALSE))
#
# mod2a <- iprobit(y = Spoligotype, year, kernel = "FBM",
#                  control = list(restarts = 8, maxit = 500,
#                                 restart.method = "error",
#                                 common.RKHS.scale = FALSE,
#                                 common.intercept = FALSE))
#
# mod2.1 <- iprobit(y = Spoligotype, period, kernel = "FBM",
#                   control = list(restarts = 8, maxit = 500,
#                                  restart.method = "error",
#                                  common.RKHS.scale = TRUE,
#                                  common.intercept = FALSE))
#
# mod2.1a <- iprobit(y = Spoligotype, period, kernel = "FBM",
#                    control = list(restarts = 8, maxit = 500,
#                                   restart.method = "error",
#                                   common.RKHS.scale = FALSE,
#                                   common.intercept = FALSE))
#
# # Spatial and temporal
# mod3 <- iprobit(y = Spoligotype, X, year, kernel = "FBM",
#                 control = list(restarts = 8, maxit = 500,
#                                restart.method = "error",
#                                common.RKHS.scale = TRUE,
#                                common.intercept = FALSE))
#
# mod3a <- iprobit(y = Spoligotype, X, year, kernel = "FBM",
#                  control = list(restarts = 8, maxit = 500,
#                                 restart.method = "error",
#                                 common.RKHS.scale = FALSE,
#                                 common.intercept = FALSE))
#
# mod3.1 <- iprobit(y = Spoligotype, X, period, kernel = "FBM",
#                   control = list(restarts = 8, maxit = 500,
#                                  restart.method = "error",
#                                  common.RKHS.scale = TRUE,
#                                  common.intercept = FALSE))
#
# mod3.1a <- iprobit(y = Spoligotype, X, period, kernel = "FBM",
#                    control = list(restarts = 8, maxit = 500,
#                                   restart.method = "error",
#                                   common.RKHS.scale = FALSE,
#                                   common.intercept = FALSE))
#
# # Save models
# save(mod1, file = "data/mod1-btb")
# save(mod1a, file = "data/mod1a-btb")
# save(mod2, file = "data/mod2-btb")
# save(mod2a, file = "data/mod2a-btb")
# save(mod2.1, file = "data/mod21-btb")
# save(mod2.1a, file = "data/mod21a-btb")
# save(mod3, file = "data/mod3-btb")
# save(mod3a, file = "data/mod3a-btb")
# save(mod3.1, file = "data/mod31-btb")
# save(mod3.1a, file = "data/mod31a-btb")
load("data/mod1-btb")
load("data/mod1a-btb")
# load("data/mod2-btb")
# load("data/mod2a-btb")
# load("data/mod21-btb")
# load("data/mod21a-btb")
load("data/mod3-btb")
load("data/mod3a-btb")
load("data/mod31-btb")
load("data/mod31a-btb")

# For reference, comparison of the models
# tab <- rbind(
#   c("spatial", logLik(mod1)),
#   c("spatial (multi)", logLik(mod1a)),
#   c("temporal", logLik(mod2)),
#   c("temporal (multi)", logLik(mod2a)),
#   c("temporal (period)", logLik(mod2.1)),
#   c("temporal (period, multi)", logLik(mod2.1a)),
#   c("spatial + temporal", logLik(mod3)),
#   c("spatial + temporal (multi)", logLik(mod3a)),
#   c("spatial + temporal (period)", logLik(mod3.1)),
#   c("spatial + temporal (period, multi)", logLik(mod3.1a))
# )
# tab <- as.tibble(tab) %>%
#   mutate(model = V1, lb = as.numeric(V2)) %>%
#   select(model, lb)

## ---- caterpillar.plot.btb ----
res1 <- summary(mod1)$tab
res1 <- cbind(param = rownames(res1), res1, model = "Spatial model")
res1[6, 1] <- "lambda[spatial]"

res2 <- summary(mod2)$tab
res2 <- cbind(param = rownames(res2), res2, model = "Temporal model")
res2[6, 1] <- "lambda[temporal]"

res3 <- summary(mod3)$tab
res3 <- cbind(param = rownames(res3), res3, model = "Spatial and temporal model")
res3[6, 1] <- "lambda[spatial]"
res3[7, 1] <- "lambda[temporal]"

dcast(res, mean + upper + lower + model + param ~ .)

as.tibble(rbind(res1, res2, res3)) %>%
  mutate(mean = as.numeric(Mean),
         lower = as.numeric(`2.5%`),
         upper = as.numeric(`97.5%`),
         sd = as.numeric(`S.D.`)) %>%
  select(param, mean, sd, lower, upper, model) %>%
  mutate(param = gsub("\\[1\\]", "\\[sp9\\]", param),
         param = gsub("\\[2\\]", "\\[sp12\\]", param),
         param = gsub("\\[3\\]", "\\[sp15\\]", param),
         param = gsub("\\[4\\]", "\\[sp20\\]", param),
         param = gsub("\\[5\\]", "\\[other\\]", param)) -> res

ggplot(res, aes(col = param)) +
  facet_grid(. ~ model, scales = "free") +
  geom_hline(yintercept = 0, linetype = "dashed", col = "grey30") +
  geom_pointrange(aes(x = param, y = mean, ymin = lower, ymax = upper)) +
  coord_flip() +
  labs(x = NULL, y = "Estimate") +
  theme_bw() +
  theme(legend.position = "none")

## ---- table.btb ----
# Function to make results table for each model
# mod1, mod1a are spatial models (single + multiple)
# mod3, mod3a are spatio-temporal models
# mod3.1, mod3.1a are spatio-period models
rest_fn <- function(mod, remove.alpha = TRUE) {
  res <- summary(mod)$tab[, c(1, 2, 1)]
  res[, 3] <- abs(res[, 3] / summary(mod)$tab[, 2])
  tmp <- res
  tmp[, 1] <- decPlac(res[, 1])
  tmp[, 2] <- decPlac(res[, 2], 3)
  tmp[, 3] <- decPlac(res[, 3], 1)
  res <- cbind(tmp, symnum(2 * pnorm(-abs(res[, 3])), corr = FALSE, na = TRUE,
                           cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                           symbols = c("***", "**", "*", ".", " ")))
  colnames(res) <- c("Estimate", "S.D.", "|$Z$|-score", "")

  mod.name <- deparse(substitute(mod))
  if (mod.name == "mod1")
    res <- rbind(res, NA)
  if (mod.name == "mod1a")
    for (i in 1:5) res <- rbind(res, NA)
  # if (mod.name == "mod2")
  #   res <- rbind(res[-nrow(res), ], NA, res[nrow(res), ])
  # if (mod.name == "mod2a") {
  #   tmp <- res[1:5, ]
  #   for (i in 1:5) tmp <- rbind(tmp, NA)
  #   res <- rbind(tmp, res[6:10, ])
  # }
  if (isTRUE(remove.alpha)) res <- res[-(1:5), ]
  res
}

rest <- cbind(rest_fn(mod1), rest_fn(mod3), rest_fn(mod3.1))
resta <- cbind(rest_fn(mod1a), rest_fn(mod3a), rest_fn(mod3.1a))
tab.btb <- rbind(rest, resta)
rownames(tab.btb) <- c(
  "Spatial", "Temporal",
  c(paste0("Spatial (", levels(btb$sp), ")"),
    paste0("Temporal (", levels(btb$sp), ")"))
)

options(knitr.kable.NA = " ")
kable(tab.btb, booktabs = TRUE, format = "latex", linesep = "", escape = FALSE,
      align = rep(c("r", "r", "r", "l"), 3),
      caption = "Results of the fitted I-probit models.") %>%
  kable_styling() %>%
  add_header_above(c(" "               = 1,
                     "Spatial"         = 3,
                     " "               = 1,
                     "Spatio-temporal" = 3,
                     " "               = 1,
                     "Spatio-period"   = 3)) %>%
  add_header_above(c(" ", "Model" = 12), bold = TRUE) %>%
  group_rows("Shared scale model", 1, 2) %>%
  group_rows("Separate scale model", 3, 12) %>%
  add_footnote(c(
    paste0("Lower-bound values  (Brier scores) for the shared scale model are ",
           decPlac(logLik(mod1), 1), " (", decPlac(get_brier_score(mod1), 3), "), ",
           decPlac(logLik(mod3), 1), " (", decPlac(get_brier_score(mod3), 3), "), and ",
           decPlac(logLik(mod3.1), 1), " (", decPlac(get_brier_score(mod3.1), 3), ") respectively."),
    paste0("Lower-bound values  (Brier scores) for the separate scale model are ",
           decPlac(logLik(mod1a), 1), " (", decPlac(get_brier_score(mod1a), 3), "), ",
           decPlac(logLik(mod3a), 1), " (", decPlac(get_brier_score(mod3a), 3), "), and ",
           decPlac(logLik(mod3.1a), 1), " (", decPlac(get_brier_score(mod3.1a), 3), ") respectively.")
  ), notation = "symbol") %>%
  landscape()

## ---- plot.btb.prep ----
# Obtain points inside the polygon
rescalee <- function(x) {
  res <- x * rep(sd.x, each = nrow(x)) + rep(mu.x, each = nrow(x))
  colnames(res) <- c("x", "y")
  res
}
X.var <- 1:2
maxmin <- cbind(apply(X, 2, min), apply(X, 2, max))
xx <- list(NULL)
for (j in 1:2) {
  mm <- maxmin[X.var[j], ]
  xx[[j]] <- seq(from = mm[1] - 1, to = mm[2] + 1, length.out = 500)
}
mm <- maxmin[X.var, ]
x.df.full <- expand.grid(xx[[1]], xx[[2]])
tmp <- x.df.full * rep(sd.x, each = nrow(x.df.full)) +
  rep(mu.x, each = nrow(x.df.full))
isin <- inside.owin(tmp[, 1], tmp[, 2], W)
x.df <- x.df.full[isin, ]

# Calculate fitted probabilities
fill.col <- iprior::ggColPal(5)
N <- nrow(x.df)
a <- predict(mod1a, list(x.df))
a1 <- predict(mod3.1a, list(x.df, matrix(rep(levels(btb$period)[1], N))))
a2 <- predict(mod3.1a, list(x.df, matrix(rep(levels(btb$period)[2], N))))
a3 <- predict(mod3.1a, list(x.df, matrix(rep(levels(btb$period)[3], N))))
a4 <- predict(mod3.1a, list(x.df, matrix(rep(levels(btb$period)[4], N))))

## ---- plot.btb ----
# Function to plot for SPATIAL MODEL
plot_spatial_model <- function(m = 1, method = "top.pieces", points = FALSE,
                               contour.labels = FALSE) {
  current.label <- paste0("Spoligotype ", gsub("Sp", "", levels(btb$sp)[m]))
  contour.df <- cbind(rescalee(x.df), prob = a$prob[, m])
  ggplot(data = contour.df, aes(x, y)) +
    geom_raster(aes(fill = prob)) -> v

  if (isTRUE(points)) {
    v <- v + geom_point(data = subset(btb, sp == levels(btb$sp)[m]), aes(x, y),
                        col = "grey15", shape = 21, fill = fill.col[m], size = 1)
  }

  if (isTRUE(contour.labels)) {
    v <- v + geom_dl(aes(z = prob, label = ..level..), stat = "contour",
                     method = list(method, cex = 0.65),
                     colour = "grey10") +
      guides(fill = FALSE)
  } else {
    v <- v + guides(fill = guide_colorbar(barwidth = 0.45, barheight = 16))
  }

  v + geom_contour(aes(z = prob), col = "grey20", size = 0.3) +
    geom_polygon(data = plot.df, aes(x, y, group = id), fill = NA,
                 col = "grey25") +
    labs(x = "Eastings (1,000 km)", y = "Northings (1,000 km)") +
    scale_x_continuous(labels = function(x) x / 1000) +
    scale_y_continuous(labels = function(x) x / 1000) +
    scale_fill_continuous(name = NULL, low = "white", high = fill.col[m],
                          limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
    theme_bw() +
    annotate("text", x = min(plot.df$x), y = max(plot.df$y),
             label = current.label, vjust = 1, hjust = 0)
}
p1 <- plot_spatial_model(1, points = FALSE, contour.labels = TRUE)
p2 <- plot_spatial_model(2, points = FALSE, contour.labels = TRUE)
p3 <- plot_spatial_model(3, points = FALSE, contour.labels = TRUE)
p4 <- plot_spatial_model(4, points = FALSE, contour.labels = TRUE)
plot_grid(p1, p2, p3, p4, ncol = 2,
          labels = c("(a)", "(b)", "(c)", "(d)"), label_size = 10,
          label_fontface = "plain")

## ---- plot.temporal.btb ----
# Function to plot for SPATIO-TEMPORAL MODEL
plot_stemporal_model <- function(year = 1, points = TRUE) {
  current.period <- levels(btb$period)[year]
  current.label <- paste0("Year: ", current.period)
  alphaa <- 0.95
  if (year == 1)
    contour.df <- cbind(rescalee(x.df), prob = a1$prob)
  if (year == 2)
    contour.df <- cbind(rescalee(x.df), prob = a2$prob)
  if (year == 3)
    contour.df <- cbind(rescalee(x.df), prob = a3$prob)
  if (year == 4)
    contour.df <- cbind(rescalee(x.df), prob = a4$prob)

  # Add first layer ------------------------------------------------------------
  p <- ggplot(contour.df, aes(x, y)) +
    geom_raster(aes(alpha = alphaa * contour.df[, 2 + 1]), fill = fill.col[1]) +
    scale_alpha_continuous(range = c(0, alphaa))

  # Add subsequent layers ------------------------------------------------------
  for (j in 2:4) {
    p <- p +
      annotate(geom = "raster", x = contour.df$x, y = contour.df$y,
               alpha = alphaa * contour.df[, 2 + j], fill = fill.col[j])
  }

  # Add decision boundaries ----------------------------------------------------
  tmp.df <- melt(contour.df, id.vars = c("x", "y"))
  p <- p +
    geom_contour(data = tmp.df, aes(x, y, z = value, col = variable), size = 1,
                 linetype = "dashed", binwidth = 0.1,
                 breaks = seq(0.5, 0.5, by = 0.1)) +
    scale_colour_manual(values = iprior::ggColPal(5)[1:4])

  # # Add contour labels ---------------------------------------------------------
  # for (j in 1:4) {
  #   p <- p +
  #     geom_dl(data = contour.df, aes(x, y, z = contour.df[, 2 + j],
  #                                    label = ..level..), col = "grey20",
  #             stat = "contour", method = list("top.pieces", cex = 0.65),
  #             breaks = seq(0.5, 0.9, by = 0.1))
  # }

  # Add points -----------------------------------------------------------------
  if (isTRUE(points)) {
    as.tibble(btb) %>%
      subset(sp != "Other" & period == current.period) -> points.df
    p <- p +
      geom_point(data = points.df, aes(x, y, fill = sp), col = "grey15",
                 shape = 21, size = 1) +
      scale_fill_manual(values = iprior::ggColPal(5)[1:4])
  }

  p +
    geom_polygon(data = plot.df, aes(x, y, group = id), fill = NA,
                 col = "grey25") +
    labs(x = "Eastings (1,000 km)", y = "Northings (1,000 km)") +
    scale_x_continuous(labels = function(x) x / 1000) +
    scale_y_continuous(labels = function(x) x / 1000) +
    theme_bw() +
    guides(colour = FALSE, alpha = FALSE,
           fill = guide_legend(nrow = 2,
                               title = "Spoligotype",
                               direction = "horizontal",
                               title.position = "top",
                               override.aes = list(size = 2))) +
    theme(legend.position = c(0.99, 0.01), legend.justification = c(1, 0),
          legend.title = element_text(size = 8), legend.title.align = 0.5) +
    # legend.box.background = element_rect(fill = NA)) +
    annotate("text", x = min(plot.df$x), y = max(plot.df$y),
             label = current.label, vjust = 1, hjust = 0)
}

p1 <- plot_stemporal_model(1)
p2 <- plot_stemporal_model(2)
p3 <- plot_stemporal_model(3)
p4 <- plot_stemporal_model(4)
plot_grid(p1, p2, p3, p4, ncol = 2,
          labels = c("(a)", "(b)", "(c)", "(d)"), label_size = 10,
          label_fontface = "plain")
