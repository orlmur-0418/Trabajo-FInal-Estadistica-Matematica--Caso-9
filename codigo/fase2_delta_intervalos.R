# Numerales 6 y 7: método delta frente a bootstrap e intervalos de confianza al 98 %

if (!file.exists("fase2_bootstrap_propiedades.R"))
  stop("Fije como directorio de trabajo la carpeta codigo (setwd).")
source("fase2_bootstrap_propiedades.R", encoding = "UTF-8")

z <- qnorm(1 - alpha / 2)

# Numeral 6. Método delta

# theta_hat es un promedio de indicadoras: Var = diag(theta) - theta theta'
# g_j(theta) = theta_j / theta_2, con gradiente (1/theta_2 en j, -theta_j/theta_2^2 en 2)
var_delta <- function(theta) {
  Sigma <- diag(theta) - theta %o% theta
  G <- sapply(idx_tau, function(j) {
    g <- numeric(k); g[j] <- 1 / theta[2]; g[2] <- -theta[j] / theta[2]^2; g })
  V <- t(G) %*% Sigma %*% G / n
  dimnames(V) <- list(nombres_tau, nombres_tau)
  V
}
# Forma cerrada: (diag(tau) + tau tau') / (n theta_2), igual a la inversa de la información de Fisher
forma_cerrada_delta <- function(theta) {
  tau <- theta[idx_tau] / theta[2]
  (diag(tau) + tau %o% tau) / (n * theta[2])
}
V_delta_generador <- var_delta(theta_verdadero)
V_delta_datos     <- var_delta(theta_hat_MV)
stopifnot(max(abs(V_delta_generador - forma_cerrada_delta(theta_verdadero))) < 1e-12,
          max(abs(V_delta_datos - forma_cerrada_delta(theta_hat_MV))) < 1e-12,
          max(abs(V_delta_datos - solve(info_fisher_tau(tau_hat, n)))) < 1e-12)

ee_delta <- sqrt(diag(V_delta_datos))
# Con los datos: Var_hat = Y_j (Y_j + Y_2) / Y_2^3
s_obs  <- Y[idx_tau] + Y[2]
pi_obs <- Y[idx_tau] / s_obs
stopifnot(max(abs(ee_delta^2 - Y[idx_tau] * s_obs / Y[2]^3)) < 1e-12,
          max(abs(ee_delta^2 - pi_obs / (s_obs * (1 - pi_obs)^3))) < 1e-12)

# Escala logarítmica: Var(log tau_hat_j) estimada con 1/Y_j + 1/Y_2
ee_log <- sqrt(1 / Y[idx_tau] + 1 / Y[2])
names(ee_delta) <- names(ee_log) <- nombres_tau

asimetria <- function(x, p = rep(1 / length(x), length(x))) {
  m <- sum(p * x); sum(p * (x - m)^3) / sum(p * (x - m)^2)^1.5
}
# Distancia de Kolmogorov entre una distribución discreta y la normal delta
dist_KS <- function(v, p, mu, sigma) {
  o <- order(v); v <- v[o]; Fv <- cumsum(p[o])
  Fn_izq <- c(0, head(Fv, -1)); Phi <- pnorm(v, mu, sigma)
  max(abs(Fv - Phi), abs(Fn_izq - Phi))
}
comparar_delta <- function(l) {
  tb <- c(NP = list(boot$NP$MV[, l]), P = list(boot$P$MV[, l]))
  de_gen <- distribucion_exacta(theta_verdadero, l); de_hat <- distribucion_exacta(theta_hat_MV, l)
  v_gen <- de_gen$valores[, "MV"]; v_hat <- de_hat$valores[, "MV"]
  sd_gen <- sqrt(V_delta_generador[l, l]); vlog <- ee_log[l]^2
  q_norm <- tau_hat[l] + c(-1, 1) * z * ee_delta[l]
  data.frame(
    parametro = nombres_tau[l],
    EE_delta_datos = ee_delta[l], EE_delta_generador = sd_gen,
    EE_exacto_plugin = subset(exacto_plugin, estimador == "MV")$EE[l],
    EE_exacto_generador = subset(exacto_theta, estimador == "MV")$EE[l],
    EE_boot_NP = sd(tb$NP), EE_boot_P = sd(tb$P),
    asim_normal = 0, asim_lognormal = (exp(vlog) + 2) * sqrt(exp(vlog) - 1),
    asim_boot_NP = asimetria(tb$NP), asim_boot_P = asimetria(tb$P),
    asim_exacta_plugin = asimetria(v_hat, de_hat$prob),
    cola_izq_NP = mean(tb$NP < q_norm[1]), cola_der_NP = mean(tb$NP > q_norm[2]),
    cola_izq_exacta_gen = sum(de_gen$prob[v_gen < tau_verdadero[l] - z * sd_gen]),
    cola_der_exacta_gen = sum(de_gen$prob[v_gen > tau_verdadero[l] + z * sd_gen]),
    KS_plugin_vs_normal = dist_KS(v_hat, de_hat$prob, tau_hat[l], ee_delta[l]),
    KS_generador_vs_normal = dist_KS(v_gen, de_gen$prob, tau_verdadero[l], sd_gen),
    row.names = NULL)
}
tabla_delta <- do.call(rbind, lapply(seq_along(nombres_tau), comparar_delta))
print(tabla_delta, digits = 4)

cuantiles_delta <- do.call(rbind, lapply(seq_along(nombres_tau), function(l) {
  p <- c(.01, .05, .5, .95, .99)
  rbind(normal_delta = qnorm(p, tau_hat[l], ee_delta[l]),
        lognormal_delta = tau_hat[l] * exp(qnorm(p) * ee_log[l]),
        bootstrap_NP = quantile(boot$NP$MV[, l], p, names = FALSE),
        bootstrap_P  = quantile(boot$P$MV[, l], p, names = FALSE))
}))
colnames(cuantiles_delta) <- paste0("q", c(1, 5, 50, 95, 99))
rownames(cuantiles_delta) <- paste(rep(nombres_tau, each = 4), rownames(cuantiles_delta))
print(round(cuantiles_delta, 4))

write.csv(tabla_delta, "../resultados/n6_delta_vs_bootstrap.csv", row.names = FALSE)
write.csv(data.frame(fila = rownames(cuantiles_delta), cuantiles_delta),
          "../resultados/n6_cuantiles_delta_bootstrap.csv", row.names = FALSE)

png("../resultados/n6_delta_vs_bootstrap.png", width = 1500, height = 520, res = 130)
par(mfrow = c(1, 3), mar = c(4.5, 4, 3, 1))
for (l in seq_along(nombres_tau)) {
  x  <- seq(min(boot$NP$MV[, l]), max(boot$NP$MV[, l]), length.out = 400)
  fn <- dnorm(x, tau_hat[l], ee_delta[l])
  fl <- dlnorm(x, log(tau_hat[l]), ee_log[l])
  d_p <- density(boot$P$MV[, l])
  hist(boot$NP$MV[, l], breaks = 60, freq = FALSE, col = "grey88", border = "white",
       main = sprintf("%s = %s / hurto", nombres_tau[l], delitos[idx_tau[l]]),
       xlab = expression(hat(tau)), ylab = "Densidad", ylim = c(0, 1.1 * max(fn, fl, d_p$y)))
  lines(d_p, col = "grey35", lwd = 1.5, lty = 2)
  lines(x, fn, col = colores["MV"], lwd = 2.5)
  lines(x, fl, col = colores["Bayes_media"], lwd = 2, lty = 4)
  abline(v = tau_hat[l], lwd = 1.5)
  legend("topright", bty = "n", cex = 0.78,
         legend = c("bootstrap NP (histograma)", "bootstrap P", "normal (delta)",
                    "lognormal (delta en log)", "tau_hat"),
         col = c("grey70", "grey35", colores[c("MV", "Bayes_media")], "black"),
         lty = c(1, 2, 1, 4, 1), lwd = c(8, 1.5, 2.5, 2, 1.5))
}
invisible(dev.off())

png("../resultados/n6_qq_bootstrap_vs_normal_delta.png", width = 1500, height = 520, res = 130)
par(mfrow = c(1, 3), mar = c(4.5, 4.5, 3, 1))
pp <- ppoints(200)
for (l in seq_along(nombres_tau)) {
  qn <- qnorm(pp, tau_hat[l], ee_delta[l])
  plot(qn, quantile(boot$NP$MV[, l], pp), pch = 16, cex = 0.5, col = colores["MV"],
       main = sprintf("%s: cuantiles bootstrap vs normal delta", nombres_tau[l]),
       xlab = "Cuantil normal delta", ylab = "Cuantil bootstrap")
  points(qn, quantile(boot$P$MV[, l], pp), pch = 1, cex = 0.5, col = "grey40")
  lines(qn, tau_hat[l] * exp(qnorm(pp) * ee_log[l]), col = colores["Bayes_media"], lwd = 2, lty = 4)
  abline(0, 1, lty = 2)
  abline(v = qnorm(probs, tau_hat[l], ee_delta[l]), lty = 3, col = "grey50")
  legend("topleft", bty = "n", cex = 0.78,
         legend = c("bootstrap NP", "bootstrap P", "lognormal (delta en log)", "identidad", "cuantiles 1 % y 99 %"),
         col = c(colores["MV"], "grey40", colores["Bayes_media"], "black", "grey50"),
         pch = c(16, 1, NA, NA, NA), lty = c(NA, NA, 4, 2, 3), lwd = c(NA, NA, 2, 1, 1))
}
invisible(dev.off())

# Numeral 7. Intervalos de confianza al 98 %

# Condicional a S_j = Y_j + Y_2 = s:  Y_j | s ~ Bin(s, pi_j), pi_j = tau_j / (1 + tau_j)
ic_tau <- function(y, m) {
  s <- y + m; t <- y / m; p <- y / s
  ee <- sqrt(y * s / m^3); el <- sqrt(1 / y + 1 / m)
  # Score (Wilson) para pi_j, llevado a tau = pi / (1 - pi)
  cen <- (p + z^2 / (2 * s)) / (1 + z^2 / s)
  rad <- z * sqrt(p * (1 - p) / s + z^2 / (4 * s^2)) / (1 + z^2 / s)
  odds <- function(q) q / (1 - q)
  # Exacto condicional (Clopper-Pearson) en la escala de tau
  cp_inf <- ifelse(y == 0, 0, y / ((m + 1) * qf(1 - alpha / 2, 2 * (m + 1), 2 * pmax(y, 1))))
  cp_sup <- (y + 1) * qf(1 - alpha / 2, 2 * (y + 1), 2 * m) / m
  list(Wald_delta = cbind(t - z * ee, t + z * ee),
       Wald_log   = cbind(ifelse(y > 0, t * exp(-z * el), NA), ifelse(y > 0, t * exp(z * el), NA)),
       Score_condicional  = cbind(odds(cen - rad), odds(cen + rad)),
       Exacto_condicional = cbind(cp_inf, cp_sup))
}

for (l in seq_along(nombres_tau)) {
  y <- Y[idx_tau[l]]; s <- s_obs[l]
  ic <- ic_tau(y, Y[2])
  pi_cp <- ic$Exacto_condicional / (1 + ic$Exacto_condicional)
  pi_sc <- ic$Score_condicional / (1 + ic$Score_condicional)
  stopifnot(abs(pbinom(y, s, pi_cp[2]) - alpha / 2) < 1e-10,
            abs(pbinom(y - 1, s, pi_cp[1], lower.tail = FALSE) - alpha / 2) < 1e-10,
            max(abs(abs(y / s - pi_sc) / sqrt(pi_sc * (1 - pi_sc) / s) - z)) < 1e-10)
}

# Percentil-t bootstrap: T* = (tau_hat* - tau_hat) / ee_delta*
pivote_t <- function(Ystar, l) {
  y <- Ystar[idx_tau[l], ]; m <- Ystar[2, ]
  (y / m - tau_hat[l]) / sqrt(y * (y + m) / m^3)
}
ic_percentil_t <- function(tipo, l) {
  q <- quantile(pivote_t(if (tipo == "NP") Ystar_np else Ystar_p, l), probs, names = FALSE)
  tau_hat[l] - rev(q) * ee_delta[l]
}

intervalos <- do.call(rbind, lapply(seq_along(nombres_tau), function(l) {
  ic <- ic_tau(Y[idx_tau[l]], Y[2])
  lims <- rbind(t(sapply(ic, as.vector)),
                Percentil_t_NP = ic_percentil_t("NP", l), Percentil_t_P = ic_percentil_t("P", l))
  data.frame(parametro = nombres_tau[l], metodo = rownames(lims), estimacion = tau_hat[l],
             inferior = lims[, 1], superior = lims[, 2], longitud = lims[, 2] - lims[, 1],
             cubre_tau_generador = lims[, 1] <= tau_verdadero[l] & tau_verdadero[l] <= lims[, 2],
             row.names = NULL)
}))
print(intervalos, digits = 4)
write.csv(intervalos, "../resultados/n7_intervalos_98.csv", row.names = FALSE)

# Cobertura exacta y longitud esperada bajo theta generador
cobertura_exacta <- function(theta, l) {
  tau <- theta[idx_tau[l]] / theta[2]
  pr  <- distribucion_exacta(theta, l)$prob
  ic  <- ic_tau(rejilla$y, rejilla$m)
  do.call(rbind, lapply(names(ic), function(met) {
    L <- ic[[met]][, 1]; U <- ic[[met]][, 2]; ok <- !is.na(L)
    data.frame(parametro = nombres_tau[l], metodo = met,
               cobertura = sum(pr[ok & L <= tau & tau <= U]),
               falla_izq = sum(pr[ok & U < tau]), falla_der = sum(pr[ok & L > tau]),
               longitud_esperada = sum(pr[ok] * (U - L)[ok]), P_no_definido = sum(pr[!ok]))
  }))
}
cobertura <- do.call(rbind, lapply(seq_along(nombres_tau), function(l) cobertura_exacta(theta_verdadero, l)))

# Percentil-t: cobertura por Monte Carlo anidado (esquema P)
set.seed(20263)
M_cob <- 2000; B_cob <- 999
cobertura_t <- do.call(rbind, lapply(seq_along(nombres_tau), function(l) {
  j <- idx_tau[l]; tau <- tau_verdadero[l]
  pr3 <- c(theta_verdadero[j], theta_verdadero[2], 1 - theta_verdadero[j] - theta_verdadero[2])
  res <- t(replicate(M_cob, {
    y3 <- rmultinom(1, n, pr3)[, 1]
    th <- y3 / n; th_t <- y3[1] / y3[2]; ee_t <- sqrt(y3[1] * (y3[1] + y3[2]) / y3[2]^3)
    Ys <- rmultinom(B_cob, n, th)
    Tst <- (Ys[1, ] / Ys[2, ] - th_t) / sqrt(Ys[1, ] * (Ys[1, ] + Ys[2, ]) / Ys[2, ]^3)
    q <- quantile(Tst, probs, names = FALSE)
    c(th_t - q[2] * ee_t, th_t - q[1] * ee_t)
  }))
  data.frame(parametro = nombres_tau[l], metodo = "Percentil_t (Monte Carlo)",
             cobertura = mean(res[, 1] <= tau & tau <= res[, 2]),
             falla_izq = mean(res[, 2] < tau), falla_der = mean(res[, 1] > tau),
             longitud_esperada = mean(res[, 2] - res[, 1]), P_no_definido = 0)
}))
cobertura <- rbind(cobertura, cobertura_t)
cobertura <- cobertura[order(cobertura$parametro), ]; rownames(cobertura) <- NULL
error_mc_cobertura <- sqrt((1 - alpha) * alpha / M_cob)
print(cobertura, digits = 4)
write.csv(cobertura, "../resultados/n7_cobertura_exacta_98.csv", row.names = FALSE)

png("../resultados/n7_intervalos_98.png", width = 1500, height = 560, res = 130)
par(mfrow = c(1, 3), mar = c(4.5, 9.5, 3, 1))
for (l in seq_along(nombres_tau)) {
  d <- subset(intervalos, parametro == nombres_tau[l]); h <- nrow(d):1
  plot(NA, xlim = range(d$inferior, d$superior, tau_verdadero[l]), ylim = c(0.5, nrow(d) + 0.5),
       yaxt = "n", ylab = "", xlab = expression(tau), main = sprintf("%s: IC al 98 %%", nombres_tau[l]))
  axis(2, at = h, labels = gsub("_", " ", d$metodo), las = 1, cex.axis = 0.8)
  segments(d$inferior, h, d$superior, h, lwd = 3, col = colores["MV"])
  points(d$estimacion, h, pch = 19, cex = 0.8)
  abline(v = tau_verdadero[l], lty = 4, col = "grey40", lwd = 2)
}
invisible(dev.off())

resultados_fase2b <- list(z = z, V_delta_datos = V_delta_datos, V_delta_generador = V_delta_generador,
                          ee_delta = ee_delta, ee_log = ee_log, tabla_delta = tabla_delta,
                          cuantiles_delta = cuantiles_delta, intervalos = intervalos,
                          cobertura = cobertura, M_cob = M_cob, B_cob = B_cob)
saveRDS(resultados_fase2b, "../resultados/resultados_fase2_delta_intervalos.rds")
