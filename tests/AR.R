library(CHAPIREG)
data("datos_anomalia_temperatura")
modelo_ARchapi <- modelo_ar(datos_anomalia_temperatura$anomalia_c, p = 2)
modelo_ARchapi
plot(modelo_ARchapi)

?modelo_ar


library(stats)
data("datos_anomalia_temperatura")
ts_data <- datos_anomalia_temperatura$anomalia_c
modelo_base <- stats::ar(ts_data,
                         aic = FALSE,
                         order.max = 2,
                         method = "yule-walker",
                         intercept = FALSE)
print(modelo_base)
