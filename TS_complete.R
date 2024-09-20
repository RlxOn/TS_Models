pacman::p_load(readr, readxl, tidyverse, lubridate, tseries, forecast, timeSeries,
               lmtest, nortest, ggplot2, gridExtra, broom, FinTS, fGarch, itsmr,
               astsa, rugarch, rmgarch, imputeTS, plotly, tsdl, babynames,
               mFilter, strucchange)

set.seed(1989) #Taylor's Version

#### Modulo para ARIMA =========================================================
source("C:/Users/raulb/Documents/Codigos_Taylors_Version/TS_ARIMA.R", 
       encoding = "UTF-8")


#### Series para trabajar ======================================================
data("AirPassengers")
data("babynames")
ts_12 <- tsdl[[12]]

df_spei <- read_xlsx("C:/Users/raulb/Documents/Datasets/SPEI.xlsx",
                     col_types = c("date", "numeric","numeric","numeric",
                                   "numeric","numeric","numeric"))

df_tarjetas <- read_xlsx("C:/Users/raulb/Documents/Datasets/Tarjetas.xlsx",
                         col_types = c("date","numeric","numeric"))

df_inpc <- read_xlsx("C:/Users/raulb/Documents/Datasets/Inflacion.xlsx",
                     col_types = c("date", "numeric"))


#Agrupamos las de SPEI y tarjetas por mes
df_spei_mes <- df_spei %>%
  group_by(Fecha = lubridate::floor_date(Fecha, "month")) %>%
  summarize(TaT = sum(TaT, na.rm = T),
            M8K = sum(M8K, na.rm = T),
            B300K = sum(B300K, na.rm = T),
            M300K = sum(M300K, na.rm = T),
            Bajo_Valor = sum(Bajo_Valor, na.rm = T),
            Alto_Valor = sum(Alto_Valor, na.rm = T)) %>%
  ungroup()

df_tarjeta_mes <- df_tarjetas %>%
  group_by(Fecha = lubridate::floor_date(Fecha, "month")) %>%
  summarize(Credito = sum(Credito, na.rm = T),
            Debito = sum(Debito, na.rm = T))

#### Imputacion de datos =======================================================
imputer <- function(ts, method = NULL, option = NULL, k = 1){
  #Inicializamos
  ts <- ts
  method <- method
  option <- option
  k <- k
  
  #Interpolacion
  if(method == "interpolation"){
    ts_complete <- imputeTS::na_interpolation(ts, option = option)
  }
  
  #Kalman
  if(method == "kalman"){
    ts_complete <- imputeTS::na_kalman(ts)
  }
  
  #Medias moviles
  if(method == "ma"){
    ts_complete <- imputeTS::na_ma(ts, k = k)
  }
  
  print(imputeTS::ggplot_na_imputations(ts, ts_complete))

  return(ts_complete)
}

spei_diario_tat <- ts(df_spei$TaT, start = 2009, frequency = 365)


spei_diario_tat <- imputer(spei_diario_tat, method = "ma", k = 7)





#### Filtro de Holdrick-Prescott ===============================================
hp_filter <- function(tseries, freq = 1000){
  #Inicializamos
  tseries <- tseries
  freq <- freq
  
  hp <- hpfilter(tseries, type = "lambda", freq = freq)
  
  time <- getTime(tseries)
  trend <- hp$trend[,1]
  cycle <- hp$cycle
  
  df <- data.frame(
    time = getTime(tseries),
    serie = tseries,
    trend = trend,
    cycle = cycle)
  
  #Transformacion para hacer el ggplot
  df_alter <- df %>%
    select(-cycle) %>%
    gather("Type", "Values", -time)
  
  print(
    ggplot(df_alter, aes(x = time, y = Values, color = Type)) +
      geom_line(linewidth = 1.3) +
      scale_color_manual(values = c(serie = "#eea38c", trend = "#8ceea3")) +
      labs(title = "Serie y filtro HP", x = NULL, y = NULL, color = NULL) +
      theme_minimal() + theme(legend.position = "bottom")
  )
  
  print(
    ggplot(df, aes(x = time, y = cycle)) + geom_line(color = "#eea38c", linewidth = 1.2) +
      labs(title = "Componente cíclica", x = NULL, y = NULL) + theme_minimal() +
      geom_hline(yintercept = 0, color = "#8cd7ee", linewidth = 1.2, linetype = "dashed")
  )
  
  return(df)

}


#### Quiebres de Bai-Perron ====================================================
bp_breaks <- function(tseries){
  #Inicializamos
  tseries <- tseries
  
  #Quiebres de bp
  bp_model <- breakpoints(tseries ~ 1)
  b_points <- c(1, bp_model$breakpoints, length(tseries))
  
  #Realizamos las rectas de regresion
  df_aux <- data.frame()

  for (i in 2:length(b_points)-1) {
    s = b_points[i]
    f = b_points[i+1]
    cat("from", s, "to", f,"\n")
    
    ltime <- getTime(tseries)
    
    serie_aux <- tseries[s:f]
    
    lmodel <- lm(serie_aux ~ getTime(serie_aux))
    
    df_aux2 <- data.frame(
      time = ltime[s:f],
      serie = serie_aux,
      line = lmodel$fitted.values
    )
    
    df_aux2 <- df_aux2 %>%
      mutate(lower = line - sd(serie_aux),
             upper = line + sd(serie_aux))
    
    df_aux <- rbind(df_aux, df_aux2)
  }
  
  #Graficamos
  print(
    ggplot() + geom_line(data = df_aux, aes(x = time, y = serie), color = "#004080",
                         linewidth = 1.3) + 
      geom_line(data = df_aux, aes(x = time, y = line), color = "#800040",
                linewidth = 1.2) + 
      geom_line(data = df_aux, aes(x = time, y = lower), color = "gray40",
                linetype = "dashed", linewidth = 1) + 
      geom_line(data = df_aux, aes(x = time, y = upper), color = "gray40",
                linetype = "dashed", linewidth = 1) + 
      theme_minimal() + xlab(NULL) + ylab(NULL) + ggtitle("Serie con quiebres de BP")
  )
  
  #Devolvemos los objetos
  return_list <- list(bp_model, b_points, df_aux)
  names(return_list) <- c("model", "breakpoints", "results")
  
  return(return_list)
  
}

#### Suavizamiento con MA ======================================================
roll_ma_sd <- function(tseries, k = 12, fill = NA, align = "center"){
  #Inicializamos los valores
  k <- k
  fill <- fill
  align <- align
  
  #Aplicamos medias moviles y sd
  ts_ma <- rollmean(tseries, k = k, fill = fill, align = align)
  ts_sd <- rollapply(tseries, width = k, FUN = sd, align = align, fill = fill)

  #Guardamos todo en un df
  df_series <- data.frame(
    Time = getTime(tseries),
    Serie = tseries,
    MA = ts_ma,
    Msd = ts_sd) %>%
    mutate(lower = MA - 1.96*Msd,
           upper = MA + 1.96*Msd)
  
  
  #Graficamos la serie con los moviles
  print(
    ggplot(data = df_series) +
      geom_line(aes(x = Time, y = Serie), linewidth = 1.2, color = "#813c60") +
      geom_line(aes(x = Time, y = MA), linewidth = 1.4, color = "#e0cc6b") +
      geom_ribbon(aes(x = Time, ymin = lower, ymax = upper), fill = "#778899", alpha = 0.3) + 
      theme_minimal() + xlab(NULL) + ylab(NULL) + ggtitle("Serie suavizada")
  )
  
  #Regresamos el df con los valores obtenidos
  return(df_series)
}

#### Suavizamiento exponencial y Holt-Winters===================================
smooth_ht <- function(tseries, h = 12, method = NULL, seasonal = "additive"){
  #Inicializamos
  tseries <- tseries
  h <- h
  method <- method
  seasonal <- seasonal
  
  #Suavizamiento exponencial
  if(method == "smooth"){
    model <- forecast::ses(tseries, h = h)
  }
  
  #Holt
  if(method == "holt"){
    model <- forecast::holt(tseries, h = h)
  }
  
  #Holt-Winters
  if(method == "hw"){
    model <- forecast::hw(tseries, h = h, seasonal = seasonal)
  }
  
  #Imprimimos el plot
  print(
    autoplot(model) + xlab(NULL) + ylab(NULL) + theme_minimal()
    )
  
  return(model)
}

