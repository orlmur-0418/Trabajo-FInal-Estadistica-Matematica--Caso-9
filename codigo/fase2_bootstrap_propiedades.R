# Numerales 4 y 5: distribución muestral bootstrap, intervalos al 98 % y
# propiedades de los estimadores de tau_j = theta_j / theta_2

if (!file.exists("fase1_modelo_simulacion_estimacion.R"))
  stop("Fije como directorio de trabajo la carpeta codigo (setwd).")
source("fase1_modelo_simulacion_estimacion.R", encoding = "UTF-8")

set.seed(2048)
B      <- 10000
alpha  <- 0.02
probs  <- c(alpha / 2, 1 - alpha / 2)
estimadores <- names(est_obs)
colores <- c(MV = "#2a78d6", Bayes_media = "#eb6834",
             Bayes_mediana = "#1baf7a", Bayes_moda = "#eda100")

# Numeral 4. Distribución muestral bootstrap e intervalos al 98 %

# Remuestras no paramétricas (NP) y paramétricas (P)
Ystar_np <- replicate(B, tabulate(sample(muestra$codigo, n, replace = TRUE), nbins = k))
Ystar_p  <- rmultinom(B, size = n, prob = theta_hat_MV)
rownames(Ystar_np) <- rownames(Ystar_p) <- delitos

remuestras_Y2_cero <- c(NP = sum(Ystar_np[2, ] == 0), P = sum(Ystar_p[2, ] == 0))
print(remuestras_Y2_cero)
stopifnot(all(remuestras_Y2_cero == 0))

verif_esquemas <- rbind(media_teorica = n * theta_hat_MV,
                        media_NP      = rowMeans(Ystar_np),
                        media_P       = rowMeans(Ystar_p),
                        var_teorica   = n * theta_hat_MV * (1 - theta_hat_MV),
                        var_NP        = apply(Ystar_np, 1, var),
                        var_P         = apply(Ystar_p, 1, var))
colnames(verif_esquemas) <- delitos
print(round(verif_esquemas, 2))

boot <- list(NP = estimar_tau(Ystar_np, a0), P = estimar_tau(Ystar_p, a0))
theta_star <- list(NP = t(Ystar_np) / n, P = t(Ystar_p) / n)

resumir_boot <- function(tipo) {
  do.call(rbind, lapply(estimadores, function(e) {
    do.call(rbind, lapply(seq_along(nombres_tau), function(l) {
      t_star <- boot[[tipo]][[e]][, l]
      t_obs  <- est_obs[[e]][1, l]
      q      <- quantile(t_star, probs, names = FALSE, type = 7)
      data.frame(bootstrap = tipo, estimador = e, parametro = nombres_tau[l],
                 estimacion = t_obs, media_boot = mean(t_star),
                 sesgo = mean(t_star) - tau_hat[l], EE = sd(t_star),
                 ECM = mean((t_star - tau_hat[l])^2),
                 perc_inf = q[1], perc_sup = q[2],
                 basico_inf = 2 * t_obs - q[2], basico_sup = 2 * t_obs - q[1],
                 cubre_tau_verdadero = q[1] <= tau_verdadero[l] & tau_verdadero[l] <= q[2])
    }))
  }))
}
resumen_boot <- rbind(resumir_boot("NP"), resumir_boot("P"))
rownames(resumen_boot) <- NULL
print(resumen_boot, digits = 4)

ic98_MV <- subset(resumen_boot, estimador == "MV",
                  c(bootstrap, parametro, estimacion, perc_inf, perc_sup, basico_inf, basico_sup))
print(ic98_MV, digits = 4)

cuantiles_MV <- do.call(rbind, lapply(seq_along(nombres_tau), function(l)
  rbind(quantile(boot$NP$MV[, l], c(.01, .05, .25, .5, .75, .95, .99)),
        quantile(boot$P$MV[, l],  c(.01, .05, .25, .5, .75, .95, .99)))))
rownames(cuantiles_MV) <- paste(rep(nombres_tau, each = 2), c("NP", "P"))
print(round(cuantiles_MV, 4))

write.csv(resumen_boot, "../resultados/n4_bootstrap_resumen_ic98.csv", row.names = FALSE)

png("../resultados/n4_bootstrap_MV_NP_vs_P.png", width = 1500, height = 520, res = 130)
par(mfrow = c(1, 3), mar = c(4.5, 4, 3, 1))
for (l in seq_along(nombres_tau)) {
  d_np <- density(boot$NP$MV[, l]); d_p <- density(boot$P$MV[, l])
  hist(boot$NP$MV[, l], breaks = 60, freq = FALSE, col = "grey90", border = "white",
       main = sprintf("%s = %s / hurto", nombres_tau[l], delitos[idx_tau[l]]),
       xlab = expression(hat(tau)^"*"), ylab = "Densidad", ylim = c(0, 1.1 * max(d_np$y, d_p$y)))
  lines(d_np, col = colores["MV"], lwd = 2)
  lines(d_p, col = colores["Bayes_media"], lwd = 2, lty = 2)
  q <- quantile(boot$NP$MV[, l], probs)
  abline(v = q, col = colores["MV"], lty = 3)
  abline(v = tau_hat[l], col = "black", lwd = 2)
  abline(v = tau_verdadero[l], col = "grey40", lwd = 2, lty = 4)
  legend("topright", bty = "n", cex = 0.8,
         legend = c("NP", "P", "IC 98% percentil (NP)", "tau_hat", "tau generador"),
         col = c(colores[c("MV", "Bayes_media", "MV")], "black", "grey40"),
         lty = c(1, 2, 3, 1, 4), lwd = c(2, 2, 1, 2, 2))
}
invisible(dev.off())

png("../resultados/n4_bootstrap_estimadores_NP.png", width = 1500, height = 520, res = 130)
par(mfrow = c(1, 3), mar = c(4.5, 4, 3, 1))
for (l in seq_along(nombres_tau)) {
  ds <- lapply(estimadores, function(e) density(boot$NP[[e]][, l]))
  plot(NA, xlim = range(sapply(ds, function(d) range(d$x))),
       ylim = c(0, max(sapply(ds, function(d) max(d$y)))),
       main = sprintf("%s: estimadores (bootstrap NP)", nombres_tau[l]),
       xlab = expression(hat(tau)^"*"), ylab = "Densidad")
  for (i in seq_along(estimadores)) lines(ds[[i]], col = colores[i], lwd = 2)
  abline(v = tau_verdadero[l], col = "grey40", lwd = 2, lty = 4)
  legend("topright", bty = "n", cex = 0.8, legend = c(estimadores, "tau generador"),
         col = c(colores, "grey40"), lty = c(1, 1, 1, 1, 4), lwd = 2)
}
invisible(dev.off())

# Numeral 5. Propiedades de los estimadores

# log10 P(Y_2 = 0): probabilidad de que tau_hat no esté definido
log10_P_Y2_cero <- c(theta_generador = n * log10(1 - unname(theta_verdadero[2])),
                     theta_hat       = n * log10(1 - unname(theta_hat_MV[2])))
print(log10_P_Y2_cero)

# Propiedades exactas por suma finita: Y_2 ~ Bin(n, theta_2), Y_j | Y_2 = m ~ Bin(n - m, theta_j / (1 - theta_2))
rejilla <- expand.grid(y = 0:n, m = 1:n)
rejilla <- rejilla[rejilla$y + rejilla$m <= n, ]

distribucion_exacta <- function(theta, l) {
  j  <- idx_tau[l]
  pr <- dbinom(rejilla$m, n, theta[2]) * dbinom(rejilla$y, n - rejilla$m, theta[j] / (1 - theta[2]))
  pr <- pr / (1 - (1 - theta[2])^n)                 # condicionado a Y_2 >= 1
  Yg <- matrix(0, k, nrow(rejilla))
  Yg[2, ] <- rejilla$m; Yg[j, ] <- rejilla$y
  Yg[setdiff(idx_tau, j)[1], ] <- n - rejilla$y - rejilla$m
  est <- estimar_tau(Yg, a0)
  list(prob = pr, valores = sapply(est, function(M) M[, l]))
}

propiedades_exactas <- function(theta, etiqueta) {
  tau <- theta[idx_tau] / theta[2]
  do.call(rbind, lapply(seq_along(nombres_tau), function(l) {
    de <- distribucion_exacta(theta, l)
    do.call(rbind, lapply(estimadores, function(e) {
      v <- de$valores[, e]; E <- sum(de$prob * v); V <- sum(de$prob * (v - E)^2)
      data.frame(fuente = etiqueta, estimador = e, parametro = nombres_tau[l],
                 esperanza = E, sesgo = E - tau[l], EE = sqrt(V), ECM = V + (E - tau[l])^2)
    }))
  }))
}
exacto_theta  <- propiedades_exactas(theta_verdadero, "exacto_theta_generador")
exacto_plugin <- propiedades_exactas(theta_hat_MV,    "exacto_plugin_theta_hat")

# E[tau_hat_j | Y_2 >= 1] en forma cerrada para el MV
E_MV_cerrada <- sapply(idx_tau, function(j) {
  m <- 1:n
  theta_verdadero[j] / (1 - theta_verdadero[2]) *
    sum((n - m) / m * dbinom(m, n, theta_verdadero[2])) / (1 - (1 - theta_verdadero[2])^n)
})
stopifnot(max(abs(E_MV_cerrada - subset(exacto_theta, estimador == "MV")$esperanza)) < 1e-10)

# Sesgo aproximado del MV (Taylor de segundo orden): tau_j / (n theta_2)
sesgo_aprox_MV <- rbind(theta_generador = tau_verdadero / (n * theta_verdadero[2]),
                        theta_hat       = tau_hat / (n * theta_hat_MV[2]))

# Monte Carlo bajo el modelo generador
set.seed(20262)
M_mc <- 20000
mc <- estimar_tau(rmultinom(M_mc, n, theta_verdadero), a0)
mc_tabla <- do.call(rbind, lapply(estimadores, function(e) do.call(rbind, lapply(seq_along(nombres_tau), function(l) {
  v <- mc[[e]][, l]
  data.frame(fuente = "monte_carlo_theta_generador", estimador = e, parametro = nombres_tau[l],
             esperanza = mean(v), sesgo = mean(v) - tau_verdadero[l], EE = sd(v),
             ECM = mean((v - tau_verdadero[l])^2))
}))))

# Cota de Cramér-Rao para tau
info_fisher_tau <- function(tau, n) {
  u <- c(tau[1], 1, tau[2], tau[3]); s <- sum(u); th <- u / s
  J <- sapply(1:3, function(l) { e <- numeric(k); e[idx_tau[l]] <- 1; (e - th) / s })
  n * t(J) %*% diag(1 / th) %*% J
}
CR_tau <- function(theta) {
  tau <- theta[idx_tau] / theta[2]
  rbind(inversa_info = diag(solve(info_fisher_tau(tau, n))),
        forma_cerrada = tau * (1 + tau) / (n * theta[2]))
}
CR_generador <- CR_tau(theta_verdadero); CR_datos <- CR_tau(theta_hat_MV)
colnames(CR_generador) <- colnames(CR_datos) <- nombres_tau
stopifnot(max(abs(CR_generador[1, ] - CR_generador[2, ])) < 1e-12,
          max(abs(CR_datos[1, ] - CR_datos[2, ])) < 1e-12)

# Cota para theta_j: coincide con Var(theta_hat_j) = theta_j (1 - theta_j) / n
CR_theta <- function(theta) {
  I_n <- n * (diag(1 / theta[1:3]) + matrix(1 / theta[4], 3, 3))
  cr  <- diag(solve(I_n))
  stopifnot(max(abs(cr - theta[1:3] * (1 - theta[1:3]) / n)) < 1e-12)
  theta * (1 - theta) / n
}
tabla_theta <- data.frame(
  delito = delitos, theta_generador = as.vector(theta_verdadero), theta_hat = as.vector(theta_hat_MV),
  EE_exacto_generador = sqrt(CR_theta(theta_verdadero)),
  EE_CR_con_datos     = sqrt(CR_theta(theta_hat_MV)),
  EE_boot_NP = apply(theta_star$NP, 2, sd), EE_boot_P = apply(theta_star$P, 2, sd),
  sesgo_boot_NP = colMeans(theta_star$NP) - theta_hat_MV,
  sesgo_boot_P  = colMeans(theta_star$P) - theta_hat_MV)
print(tabla_theta, digits = 4)

tabla_CR <- data.frame(
  parametro = nombres_tau,
  CR_theta_generador = CR_generador[2, ], sqrt_CR_generador = sqrt(CR_generador[2, ]),
  CR_con_datos = CR_datos[2, ], sqrt_CR_datos = sqrt(CR_datos[2, ]),
  Var_exacta_MV = subset(exacto_theta, estimador == "MV")$EE^2,
  eficiencia_CR_sobre_Var = CR_generador[2, ] / subset(exacto_theta, estimador == "MV")$EE^2)
print(tabla_CR, digits = 4)

# Comparación de sesgo, EE y ECM según la fuente
boot_tabla <- with(resumen_boot, data.frame(fuente = paste0("bootstrap_", bootstrap),
                   estimador, parametro, esperanza = media_boot, sesgo, EE, ECM))
comparacion <- rbind(exacto_theta, mc_tabla, exacto_plugin, boot_tabla)
comparacion$fuente <- factor(comparacion$fuente,
  levels = c("exacto_theta_generador", "monte_carlo_theta_generador",
             "exacto_plugin_theta_hat", "bootstrap_NP", "bootstrap_P"))
comparacion <- comparacion[order(comparacion$parametro, comparacion$estimador, comparacion$fuente), ]
rownames(comparacion) <- NULL
for (p in nombres_tau) {
  cat("\n====", p, "====\n")
  print(subset(comparacion, parametro == p, -parametro), digits = 4, row.names = FALSE)
}
print(round(sesgo_aprox_MV, 6))
write.csv(comparacion, "../resultados/n5_comparacion_sesgo_EE_ECM.csv", row.names = FALSE)
write.csv(tabla_CR, "../resultados/n5_cramer_rao_tau.csv", row.names = FALSE)
write.csv(tabla_theta, "../resultados/n5_propiedades_theta.csv", row.names = FALSE)

png("../resultados/n5_distribucion_exacta_vs_bootstrap_MV.png", width = 1500, height = 520, res = 130)
par(mfrow = c(1, 3), mar = c(4.5, 4, 3, 1))
for (l in seq_along(nombres_tau)) {
  de <- distribucion_exacta(theta_verdadero, l)
  o  <- order(de$valores[, "MV"]); x <- de$valores[o, "MV"]; Fx <- cumsum(de$prob[o])
  rango <- quantile(boot$NP$MV[, l], c(.0005, .9995))
  plot(x, Fx, type = "s", col = "grey40", lwd = 2, xlim = rango,
       main = sprintf("%s: FDA de tau_hat_MV", nombres_tau[l]), xlab = expression(hat(tau)), ylab = "F")
  lines(ecdf(boot$NP$MV[, l]), col = colores["MV"], do.points = FALSE, lwd = 1.5)
  lines(ecdf(boot$P$MV[, l]), col = colores["Bayes_media"], do.points = FALSE, lwd = 1.5, lty = 2)
  abline(v = c(tau_verdadero[l], tau_hat[l]), lty = c(4, 1), col = c("grey40", "black"))
  legend("bottomright", bty = "n", cex = 0.8,
         legend = c("exacta (theta generador)", "bootstrap NP", "bootstrap P"),
         col = c("grey40", colores[c("MV", "Bayes_media")]), lty = c(1, 1, 2), lwd = 2)
}
invisible(dev.off())

resultados_fase2 <- list(B = B, alpha = alpha, resumen_boot = resumen_boot,
                         comparacion = comparacion, cramer_rao = tabla_CR,
                         propiedades_theta = tabla_theta, sesgo_aprox_MV = sesgo_aprox_MV)
saveRDS(resultados_fase2, "../resultados/resultados_fase2.rds")
