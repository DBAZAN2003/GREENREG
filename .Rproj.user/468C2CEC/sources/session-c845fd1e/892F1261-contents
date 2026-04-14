
mariel <- read.csv("C:/Users/dbaza/Documents/datos_experimento.csv")
head(mariel)

mariel_limpio <- na.omit(mariel)
mariel_limpio

library(lme4)

# Asegurar formato correcto
mariel_limpio$sustrato <- as.factor(mariel_limpio$sustrato)
mariel_limpio$id_planta <- as.factor(mariel_limpio$id_planta)
mariel_limpio$medicion <- as.numeric(mariel_limpio$medicion)

modelo <- lmer(altura ~ sustrato * medicion + (1|id_planta),
               data = mariel_limpio)
modelo

summary(modelo)
anova(modelo)


library(ggplot2)

ggplot(mariel_limpio,
       aes(x = medicion, y = altura, color = sustrato)) +
  stat_summary(fun = mean, geom = "line", linewidth = 1.3) +
  stat_summary(fun = mean, geom = "point", size = 2) +
  theme_minimal() +
  labs(title = "Crecimiento promedio por sustrato",
       x = "Semana / Medición",
       y = "Altura")
