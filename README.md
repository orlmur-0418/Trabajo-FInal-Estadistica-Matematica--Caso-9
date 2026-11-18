# Inferencia sobre razones de proporciones de delitos

<b>Trabajo final · Estadística Matemática · Caso 9</b><br>
Universidad del Valle · Escuela de Estadística<br>
Orlaín Murillo · Leandro Vélez

---

## Anexo: video de presentación

<p align="center">
  <a href="https://youtu.be/tekDNd4a32k">
    <img src="https://img.youtube.com/vi/tekDNd4a32k/hqdefault.jpg" alt="Ver el video del trabajo final en YouTube" width="640">
  </a>
</p>

<p align="center">
  <a href="https://youtu.be/tekDNd4a32k"><b>▶ Ver el video en YouTube: youtu.be/tekDNd4a32k</b></a>
</p>

Enlace directo: https://youtu.be/tekDNd4a32k

| Minuto | Contenido |
|---|---|
| [0:00](https://youtu.be/tekDNd4a32k?t=0) | Apertura y parámetro de interés |
| [0:31](https://youtu.be/tekDNd4a32k?t=31) | Numeral 1. Modelo estadístico |
| [1:46](https://youtu.be/tekDNd4a32k?t=106) | Numeral 2. Simulación |
| [3:10](https://youtu.be/tekDNd4a32k?t=190) | Numeral 3. Estimadores |
| [3:53](https://youtu.be/tekDNd4a32k?t=233) | Numeral 4. Bootstrap e intervalos |
| [4:34](https://youtu.be/tekDNd4a32k?t=274) | Numeral 5. Propiedades de los estimadores |
| [4:53](https://youtu.be/tekDNd4a32k?t=293) | Numeral 6. Método delta |
| [5:40](https://youtu.be/tekDNd4a32k?t=340) | Numeral 7. Intervalos de confianza del 98 % |
| [6:27](https://youtu.be/tekDNd4a32k?t=387) | Numeral 8. Pruebas de hipótesis |
| [7:09](https://youtu.be/tekDNd4a32k?t=429) | Numeral 9. Prueba uniformemente más potente |
| [8:04](https://youtu.be/tekDNd4a32k?t=484) | Comparación de las cuatro pruebas |
| [8:49](https://youtu.be/tekDNd4a32k?t=529) | Conclusiones |

---

## Estructura del repositorio

```
codigo/        scripts en R (la cadena se ejecuta en orden: fase 1 -> fase 2 -> fase 3)
resultados/    tablas (.csv), gráficas (.png) y objetos de R (.rds)
```

| Script | Numerales |
|---|---|
| `fase1_modelo_simulacion_estimacion.R` | 1 a 3: modelo, simulación y estimadores (MV, momentos y Bayes) |
| `fase2_bootstrap_propiedades.R` | 4 y 5: bootstrap, intervalos al 98 %, sesgo, error estándar, ECM y Cramér-Rao |
| `fase2_delta_intervalos.R` | 6 y 7: método delta e intervalos de confianza con su cobertura exacta |
| `fase3_pruebas_hipotesis.R` | 8 y 9: pruebas de Wald, potencia, errores tipo I y II y prueba UMP |
| `figuras_articulo.R` | Figuras 1 y 2 del artículo |
| `razon_extorsion_asesinato.R` | Cálculos complementarios citados en el artículo |
