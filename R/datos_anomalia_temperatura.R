#' Anomalías de temperatura media diaria
#'
#' Conjunto de datos que muestra las anomalías (diferencias respecto a la media histórica)
#' de la temperatura media diaria registradas a lo largo del tiempo. Este dataset puede
#' utilizarse para ilustrar análisis de series de tiempo y modelos autorregresivos o ARIMA.
#'
#' @format Un data frame con 730 observaciones y 2 variables:
#' \describe{
#'   \item{fecha}{Fecha del registro (AAAA-MM-DD)}
#'   \item{anomalia_c}{Anomalía de temperatura respecto al promedio histórico (°C)}
#' }
#'
#' @details Los datos fueron construidos de manera simulada a partir de patrones reales
#' observados en bases climáticas del INEGI y de la CONAGUA. La creación propia de esta
#' base permite garantizar que los valores sean adecuados para fines didácticos y
#' compatibles con las funciones del paquete \pkg{CHAPIREG}.
#'
#' @examples
#' # Probar un modelo Autorregresivo de orden 2
#' modelo_ar(anomalia_temperatura$anomalia_c, p = 2)
#'
#' # Probar un modelo ARIMA(1,1,1)
#' modelo_arima(anomalia_temperatura$anomalia_c, p = 1, d = 1, q = 1)
#'
#' @source Elaboración propia con fines didácticos a partir de información pública de INEGI y CONAGUA.
"datos_anomalia_temperatura"

