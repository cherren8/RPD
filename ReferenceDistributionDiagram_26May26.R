## Make a schematic of empirical reference distributions for new Fig. 1

library(ggplot2)

# Find the most abundant taxain GBH and make a plot of their ratio

br.ref <- ab.ref
gbh.abun <- apply(br.ref, 2, mean)
tax1 <- order(gbh.abun, decreasing = T)[3]
tax2 <- order(gbh.abun, decreasing = T)[6]

tax1.ab <- br.ref[, tax1]
tax2.ab <- br.ref[, tax2]

par(mfrow = c(1, 1))
hist(log(tax1.ab / tax2.ab), breaks = 30)

log.rat <- log(tax1.ab / tax2.ab)

tax.df <- cbind(tax1.ab, tax2.ab, log.rat)

tax.df <- tax.df[is.finite(log.rat), ]

# Make a density plot
ggplot(tax.df, aes(x = log.rat)) +
  geom_density( outline.type = "upper",
                lineend = "butt",
                linejoin = "round",
                linemitre = 10,
                fill = "blue",
                alpha = .2,
                adjust = .7) +
  xlim(-6, 6) +
  theme( axis.text.x = element_blank(),
         axis.ticks.x = element_blank(),
        # axis.title.x = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.85, 0.65),
        legend.title = element_text(size = 15), # Larger title
        legend.text = element_text(size = 12) ,
        panel.grid = element_blank(),
        plot.margin = margin(20, 20, 20, 20, "pt"),
        axis.title.y = element_text(size = 20),
       axis.title.x = element_text(size = 20, margin = margin(t = 15)),
        axis.text.y = element_text(size = 15)) +
  geom_segment(aes(x = median(log.rat[is.finite(log.rat)]),
                   xend = median(log.rat[is.finite(log.rat)]),
                   y = 0,
                   yend = .305),
               linetype = "dashed",
               color = "grey30", size = 1) +
  labs(y = "Density", x = "Reference Sample Values", colour = "Test Values") +
  # geom_segment(aes(x = quantile(log.rat[is.finite(log.rat)], .06),
  #                  xend = quantile(log.rat[is.finite(log.rat)], .06),
  #                  y = -0.003,
  #                  yend = 0.004),
  #              linetype = "dashed",
  #              color = "purple", size = 1) +
  # geom_segment(aes(x = quantile(log.rat[is.finite(log.rat)], .8),
  #                  xend = quantile(log.rat[is.finite(log.rat)], .8),
  #                  y = -0.003,
  #                  yend = 0.004),
  #              linetype = "dashed",
  #              color = "forestgreen", size = 1) +
  scale_y_continuous(expand = c(0.03, 0)) +
  # Optional: Customize the color manually
  scale_color_manual(values = c("Sample 1" = "purple", "Sample 2" = "maroon")) +
  geom_point(aes(x = quantile(log.rat[is.finite(log.rat)], .81),
                 y = 0, colour = "Sample 2"),
             shape = 8,
             size = 4,
             stroke = 1.2) +
  geom_point(aes(x = quantile(log.rat[is.finite(log.rat)], .06),
                 y = 0, colour = "Sample 1"),
             shape = 8,
             size = 4,
             stroke = 1.2) +
  geom_segment(aes(x = quantile(log.rat[is.finite(log.rat)], .09),
                   y = 0,
                   xend = quantile(log.rat[is.finite(log.rat)], .45),
                   yend = 0),
               arrow = arrow(length = unit(0.2, "cm")),
               colour = "purple") +
  geom_segment(aes(x = quantile(log.rat[is.finite(log.rat)], .76),
                   y = 0,
                   xend = quantile(log.rat[is.finite(log.rat)], .55),
                   yend = 0),
               arrow = arrow(length = unit(0.2, "cm")),
               colour = "maroon") +
  geom_ribbon(
    data = {
      d <- density(log.rat[is.finite(log.rat)], adjust = .7)
      q06 <- quantile(log.rat[is.finite(log.rat)], .06)
      q50 <- quantile(log.rat[is.finite(log.rat)], .50)
      data.frame(x = d$x, y = d$y) |> subset(x >= q06 & x <= q50)
    },
    aes(x = x, ymin = 0, ymax = y),
    fill = "purple",
    alpha = 0.32
  ) +
  # Shade 50th to 80th in blue
  geom_ribbon(
    data = {
      d <- density(log.rat[is.finite(log.rat)], adjust = .7)
      q50 <- quantile(log.rat[is.finite(log.rat)], .50)
      q81 <- quantile(log.rat[is.finite(log.rat)], .81)
      data.frame(x = d$x, y = d$y) |> subset(x >= q50 & x <= q81)
    },
    aes(x = x, ymin = 0, ymax = y),
    fill = "maroon",
    alpha = 0.32
  ) +
  annotate("text", x = quantile(log.rat[is.finite(log.rat)], .2), y = .05, label = "0.44", color = "purple", size = 6) +
  annotate("text", x = quantile(log.rat[is.finite(log.rat)], .7), y = .05, label = "0.31", color = "maroon", size = 6)



