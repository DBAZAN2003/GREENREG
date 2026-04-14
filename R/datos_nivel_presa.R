#' Niveles de agua en presas de México
#'
#' Conjunto de datos que presenta los niveles de almacenamiento de agua
#' registrados en distintas presas y fechas, útil para análisis de series de tiempo
#' y modelos ARIMA.
#'
#' @format Un data frame con 1825 observaciones y 2 variables:
#' \describe{
#'   \item{fecha}{Fecha del registro (AAAA-MM-DD)}
#'   \item{nivel_m}{Nivel de almacenamiento de la presa (metros)}
#' }
#'
#' @details Los datos fueron simulados a partir de patrones observados en reportes
#' de CONAGUA, manteniendo coherencia temporal y realismo. Se generaron con fines
#' didácticos para demostrar la aplicación de modelos AR, MA y ARIMA en \pkg{CHAPIREG}.
#'
#' @examples
#' # Visualizar la serie
#' plot(nivel_presa$fecha, nivel_presa$nivel_m, type = "l")
#'
#' # Ajustar un modelo ARIMA(1,1,1)
#' modelo_complejo <- modelo_arima(nivel_presa$nivel_m, p = 1, d = 1, q = 1)
#'
#' @source Elaboración propia con fines didácticos basada en patrones de datos de CONAGUA.
"datos_nivel_presa"

