# Figuras 1 y 2 del artículo (escala de grises, PDF)

source("fase3_pruebas_hipotesis.R", encoding = "UTF-8")
out <- "../resultados"

leyenda_inferior <- function(legend, cex, seg.len, ncol = length(legend), ...) {
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  ancho <- sapply(split(strwidth(legend, cex = cex), ceiling(seq_along(legend) * ncol / length(legend))), max)
  legend("bottom", legend = legend, ncol = ncol, bty = "n", xpd = NA, cex = cex,
         seg.len = seg.len, text.width = ancho, x.intersp = 0.6, y.intersp = 0.95, ...)
}

# Figura 1: distribución bootstrap del MV y aproximaciones normal y lognormal
tit <- c("Asesinato / hurto", "Secuestro / hurto", "Extorsión / hurto")
ejex <- c(expression(hat(tau)[1]), expression(hat(tau)[3]), expression(hat(tau)[4]))
pdf(file.path(out, "fig1_distribucion_muestral.pdf"), width = 4.92, height = 1.7, pointsize = 8)
par(mfrow = c(1, 3), mar = c(3.6, 3.1, 1.6, 0.4), mgp = c(2.5, 0.7, 0), oma = c(2.2, 0, 0, 0),
    cex.axis = 1.05, cex.lab = 1.3, cex.main = 1.25)
for (l in 1:3) {
  b <- boot$NP$MV[, l]
  x <- seq(min(b), max(b), length.out = 400)
  fn <- dnorm(x, tau_hat[l], ee_delta[l]); fl <- dlnorm(x, log(tau_hat[l]), ee_log[l])
  hist(b, breaks = 30, freq = FALSE, col = "grey82", border = "white", main = tit[l],
       xlab = ejex[l], ylab = "", ylim = c(0, 1.1 * max(fn, fl)))
  if (l == 1) title(ylab = "Densidad", line = 2.0)
  lines(x, fn, lwd = 1.6, lty = 1); lines(x, fl, lwd = 1.6, lty = "42")
  abline(v = tau_hat[l], lty = 3, lwd = 1.2)
}
leyenda_inferior(c("Bootstrap", "Normal", "Lognormal", "Estimación MV"), cex = 1.3, seg.len = 2.4,
                 col = c("grey75", "black", "black", "black"), lty = c("solid", "solid", "42", "13"),
                 lwd = c(7, 1.6, 1.6, 1.2))
dev.off()

# Figura 2: potencia para tau_4 y error tipo I en la frontera
lty_p <- c("13", "63", "solid", "33"); lwd_p <- c(2.2, 1.5, 1.3, 1.3)
col_g <- c("black", "grey45", "black", "black")
etq <- c("Wald", "Wald log", "UMP aleatorizada", "UMP no aleatorizada")
pdf(file.path(out, "fig2_potencia_tamano.pdf"), width = 4.92, height = 2.15, pointsize = 8)
par(mfrow = c(1, 2), mar = c(3.5, 3.4, 3.1, 0.5), mgp = c(2.2, 0.7, 0), oma = c(3.3, 0, 0, 0),
    cex.axis = 1.05, cex.lab = 1.2, cex.main = 1.2)
matplot(tau_rejilla[[3]], curvas[[3]], type = "l", lty = lty_p, col = col_g, lwd = lwd_p,
        xlim = c(0.5, 1.15), ylim = c(0, 1), xlab = expression(tau[4]), ylab = "Potencia", main = "")
title(expression(bold("(a) Potencia para ") * tau[4]), line = 2.0)
abline(h = alpha, col = "grey55"); abline(v = c(tau_verdadero[3], tau_hat[3], 1), lty = 3, col = "grey35")
mtext(c(expression(tau[4]), expression(hat(tau)[4]), "1"),
      side = 3, at = c(tau_verdadero[3] - 0.012, tau_hat[3] + 0.012, 1), line = 0.15, cex = 1.05)
text(0.5, alpha + 0.06, expression(alpha == 0.02), adj = 0, cex = 1.1)
matplot(psi_rejilla, tam_psi, type = "l", lty = lty_p, col = col_g, lwd = lwd_p,
        ylim = c(0, max(tam_psi) * 1.05), xlab = expression(psi[j] == theta[j] + theta[2]),
        ylab = "Error tipo I", main = "")
title(expression(bold("(b) Error tipo I en ") * tau[j] == 1), line = 2.0)
abline(h = alpha, col = "grey55"); abline(v = psi_hat, lty = 3, col = "grey35")
mtext(c(expression(hat(psi)[3]), expression(hat(psi)[1]), expression(hat(psi)[4])), side = 3,
      at = c(psi_hat[2], psi_hat[1] - 0.025, psi_hat[3] + 0.03), line = 0.1, cex = 1.05)
text(1, alpha - 0.013, expression(alpha == 0.02), adj = 1, cex = 1.1)
leyenda_inferior(etq, cex = 1.25, seg.len = 3, ncol = 2, lty = lty_p, col = col_g, lwd = lwd_p)
dev.off()
