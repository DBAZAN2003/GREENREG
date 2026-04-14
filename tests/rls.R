library(GREENREG)
data("datos_rendimiento_maiz")
modelo <- rls(rendimiento_maiz_ton_ha ~ precipitacion_mm, data = datos_rendimiento_maiz)
modelo
plot(modelo)

?rls


library(lmtest)
library(car)
data("datos_rendimiento_maiz")
modelo_lm <- lm(rendimiento_maiz_ton_ha ~ precipitacion_mm, data = datos_rendimiento_maiz)
summary(modelo_lm)
shapiro_test <- shapiro.test(residuals(modelo_lm))
print(shapiro_test)
bp_test <- bptest(modelo_lm)
print(bp_test)
dw_test <- durbinWatsonTest(modelo_lm)
print(dw_test)
