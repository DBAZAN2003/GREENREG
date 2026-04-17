library(GREENREG)
data("datos_rendimiento_maiz")
modelo_rlm <- rlm(rendimiento_maiz_ton_ha ~., data = datos_rendimiento_maiz)
modelo_rlm
plot(modelo_rlm)

?rlm

library(lmtest)
library(car)
data("datos_rendimiento_maiz")
modelo_base <- lm(rendimiento_maiz_ton_ha ~ ., data = datos_rendimiento_maiz)
print(summary(modelo_base))
cat("\nShapiro-Wilk (Normalidad):\n")
print(shapiro.test(residuals(modelo_base)))
cat("\nBreusch-Pagan (Homocedasticidad):\n")
print(bptest(modelo_base))
cat("\nDurbin-Watson (Autocorrelación):\n")
print(durbinWatsonTest(modelo_base))
plot(modelo_base)

?lm
