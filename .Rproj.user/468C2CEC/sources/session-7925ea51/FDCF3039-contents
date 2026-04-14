library(CHAPIREG)
data("datos_anomalia_temperatura")
ts_data <- datos_anomalia_temperatura$anomalia_c

modelo_ARMAchapi <- modelo_arma(ts_data, p = 1, q = 1, include_mean = TRUE)
modelo_ARMAchapi
plot(modelo_ARMAchapi)

?modelo_arma

library(stats)
data("datos_anomalia_temperatura")
ts_data <- datos_anomalia_temperatura$anomalia_c
modelo_base <- stats::arima(ts_data,
                            order = c(1, 0, 1),    # p=1, d=0, q=1
                            method = "CSS",
                            include.mean = TRUE)
modelo_base
cat("\nLjung-Box (Residuos):\n")
print(Box.test(residuals(modelo_base), lag = 10, type = "Ljung-Box"))
cat("\nShapiro-Wilk (Residuos):\n")
print(shapiro.test(residuals(modelo_base)))
