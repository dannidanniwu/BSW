library(ggplot2)
dp <- dd[, .(avg = mean(y)), keyby = .(A, site, k)]

ggplot(data = dp, aes(x = k, y = avg)) +
  geom_line(aes(group = site, color = factor(A))) +
  scale_color_manual(values = c("#a53e3e", "#3ea5a5"),
                     guide = guide_legend(reverse = TRUE),
                     labels = c("no treatment", "treatment")) +
  ylab("average Y") +
  xlab("period (k)") +
  theme(panel.grid = element_blank(),
        legend.title = element_blank())

ggplot(data = dd, aes(x = k, y = y)) +
  geom_point(aes(color = factor(A, labels = c("Control", "Intervention"))),
             size = 0.1) +
  scale_color_manual(values = c("#d07b7c", "#7ba7d0")) +
  facet_wrap(~site, ncol = 8) +
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        axis.text.y = element_text(size = 7),
        strip.text = element_text(size = 8),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 2)))

dSum <- dd[, .(Y = mean(y)), keyby = .(site, k, A, startTrt)]

ggplot(data = dSum, 
       aes(x = k, y = Y, group = interaction(site, A))) +
  geom_line(aes(color = factor(A))) +
  facet_grid(startTrt~.) +
  scale_x_continuous(breaks = seq(0, 19, by = 1), name = "week") +
  scale_color_manual(values = c("#b8cce4", "#4e81ba")) +
  theme(panel.grid = element_blank(),
        legend.position = "none") 