library(GREENREG)
data("datos_tala_ilegal_poisson")
modelo_poisson_chapi <- reg_poisson(numero_talas ~ densidad_poblacion + superficie_forestal_ha + dias_vigilancia, data = datos_tala_ilegal_poisson)
modelo_poisson_chapi
plot(modelo_poisson_chapi)


library(stats)
data("datos_tala_ilegal_poisson")
modelo_base <- glm(numero_talas ~ .,
                   data = datos_tala_ilegal_poisson,
                   family = poisson(link = "log"))
print(summary(modelo_base))
cat("\n--- Métricas del Modelo (glm) ---\n")
cat("Deviance:", deviance(modelo_base), "\n")
cat("AIC:", AIC(modelo_base), "\n")
