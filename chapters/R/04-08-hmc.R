# Chapter 4
# Hamiltonian Monte Carlo example
chapter.no <- "04"

sho <- function(t, y, parameters) {
  x <- y[1]
  p <- y[2]
  m <- parameters
  dy <- numeric(2)
  dy[1] <- p / m
  dy[2] <- -x
  list(dy)
}
col_pal <- iprior::gg_col_hue(2)
is_even <- function(x) x %% 2 == 0
chapter.no <- "04"

# Phase portrait 1

pp <- function(no = 1) {
  tmp1 <- phaseR::flowField(sho, mgp = c(1.8, 0.7, 0), #main = "Initialise",
                            x.lim = c(-2.6, 2.1), xlab = "Position (x)",
                            y.lim = c(-3, 3), ylab = "Momentum (z)",
                            points = 19, add = FALSE, parameters = 2)
  grid()
  abline(h = 0, col = "gray", lwd = 2, lty = 3)
  points(2, 0, lwd = 3, col = col_pal[1])

  if (no > 1) {
    tmp2 <- phaseR::trajectory(sho, y0 = c(2,0), t.end = 1,
                               parameters = 1, colour = col_pal[1], lwd = 2)
    x <- tmp2$x[length(tmp2$x)]
    y <- tmp2$y[length(tmp2$x)]
    points(x, y, lwd = 2, col = col_pal[1])
    lines(c(x, x), c(y, 0), lwd = 2, lty = 3)
    points(x, 0, lwd = 2)
  }

  if (no > 2) {
    lines(c(x, x), c(0, 1), lwd = 2, col = col_pal[2])
    points(x, 1, lwd = 3, col = col_pal[1])
  }

  if (no > 3) {
    tmp2 <- phaseR::trajectory(sho,
                               y0 = c(x, 1),
                               t.end = 3,
                               parameters = 1,
                               colour = col_pal[1], lwd = 2)
    x <- tmp2$x[length(tmp2$x)]
    y <- tmp2$y[length(tmp2$x)]
    points(x, y, lwd = 2, col = col_pal[1])
    lines(c(x, x), c(y, 0), lwd = 2, lty = 3)
    points(x, 0, lwd = 2)
  }

  if (no > 4) {
    lines(c(x, x), c(0, 2.8), lwd = 2, col = col_pal[2])
    points(x, 2.8, lwd = 3, col = col_pal[1])
  }

  if (no > 5) {
    tmp2 <- phaseR::trajectory(sho,
                               y0 = c(x, 2.8),
                               t.end = 0.5,
                               parameters = 1,
                               colour = col_pal[1], lwd = 2)
    x <- tmp2$x[length(tmp2$x)]
    y <- tmp2$y[length(tmp2$x)]
    points(x, y, lwd = 2, col = col_pal[1])
    lines(c(x, x), c(y, 0), lwd = 2, lty = 3)
    points(x, 0, lwd = 2)
  }

  if (no > 6) {
    lines(c(x, x), c(0, -2.5), lwd = 2, col = col_pal[2])
    points(x, -2.5, lwd = 3, col = col_pal[1])
  }

  if (no > 7) {
    tmp2 <- phaseR::trajectory(sho,
                               y0 = c(x, -2.5),
                               t.end = 1.5,
                               parameters = 1,
                               colour = col_pal[1], lwd = 2)
    x <- tmp2$x[length(tmp2$x)]
    y <- tmp2$y[length(tmp2$x)]
    points(x, y, lwd = 2, col = col_pal[1])
    lines(c(x, x), c(y, 0), lwd = 2, lty = 3)
    points(x, 0, lwd = 2)
  }

  if (is_even(no)) {
    mtext(paste0(no - 1, ": Metropolis update"), side = 3, adj = 0, line = 1.2,
          cex = 1.25)
  } else {
    if (no == 1) {
      mtext("0: Initialise", side = 3, adj = 0, line = 1.2, cex = 1.25)
    } else {
      mtext(paste0(no - 1, ": Perturb momentum"), side = 3, adj = 0, line = 1.2, cex = 1.25)
    }
  }
}

ww <- 6; hh <- 5
pdf("figure/04-phase1.pdf", width = ww, height = hh)
pp(1)
dev.off()
pdf("figure/04-phase2.pdf", width = ww, height = hh)
pp(2)
dev.off()
pdf("figure/04-phase3.pdf", width = ww, height = hh)
pp(3)
dev.off()
pdf("figure/04-phase4.pdf", width = ww, height = hh)
pp(4)
dev.off()
pdf("figure/04-phase5.pdf", width = ww, height = hh)
pp(5)
dev.off()
pdf("figure/04-phase6.pdf", width = ww, height = hh)
pp(6)
dev.off()
pdf("figure/04-phase7.pdf", width = ww, height = hh)
pp(7)
dev.off()
pdf("figure/04-phase8.pdf", width = ww, height = hh)
pp(8)
dev.off()
move_fig_to_chapter()
