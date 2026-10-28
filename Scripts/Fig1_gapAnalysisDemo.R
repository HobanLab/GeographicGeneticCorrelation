## ---------------------------------------------------------------------------
## Gap analysis figure for Quercus georgiana
##   Panel A  occurrence map: 50-km buffers, hatched where the occurrence is
##            represented in ex situ collections
##   Panel B  geographic and ecological coverage, with the unfilled part of each
##            column showing the conservation gap
##
## Requires: sf.  Occurrence coordinates, ecoregion outlines and coverage values
## below are illustrative placeholders - swap in real data before publishing.
## ---------------------------------------------------------------------------

library(sf)
sf_use_s2(FALSE)

out_tiff <- "Fig1_gapAnalysisDemo.tiff"
dpi      <- 600

## ---- inputs ---------------------------------------------------------------
buffer_km <- 50
cov_geo   <- 27      # geographic coverage (%)
cov_eco   <- 44      # ecological coverage (%)
crs_map   <- 5070    # NAD83 / Conus Albers (equal area, metres)

ext_ll <- c(xmin = -88.6, ymin = 29.9, xmax = -78.6, ymax = 36.6)

## occurrences: sampled = wild source of an ex situ living specimen
occ <- data.frame(
  lon = c(-84.17, -84.06, -83.98, -84.35, -83.72, -84.52, -83.55, -82.55, -82.20, -83.10, -84.62, -84.90, -85.80, -86.15, -86.55, -85.35, -81.55, -81.05, -80.35, -84.75),
  lat = c(33.81, 33.63, 33.90, 33.98, 33.72, 34.12, 34.05, 33.55, 33.48, 33.35, 33.55, 33.90, 33.48, 33.65, 33.35, 32.95, 35.20, 34.55, 34.90, 34.75),
  sampled = c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE)
)

## schematic ecoregion outlines (lon, lat pairs)
eco_defs <- list(
  list(name = "Interior Plateau", col = "#E7A6AE", xy = c(
    -88.239,33.896, -91.5,36.936, -91.5,39.5, -80.825,39.5, -88.239,33.896)),
  list(name = "Ridge and Valley", col = "#C6A8DA", xy = c(
    -87.21,33.396, -87.967,34.101, -80.825,39.5, -79.134,39.5, -87.21,33.396)),
  list(name = "Blue Ridge", col = "#A98BC8", xy = c(
    -84.425,34.284, -85.146,34.956, -79.134,39.5, -77.524,39.5, -84.425,34.284)),
  list(name = "Piedmont", col = "#E4897F", xy = c(
    -88.5,32.05, -89.0,32.043, -85.146,34.956, -84.425,34.284, -77.524,39.5, -76.0,39.5,
    -76.0,36.6, -79.6,34.95, -81.0,34.0, -82.2,33.4, -83.6,32.9, -84.95,32.55, -85.9,32.45,
    -86.093,32.424, -86.448,32.755, -87.063,32.289, -88.5,32.05)),
  list(name = "Southeastern Plains", col = "#F2C489", xy = c(
    -91.5,36.936, -88.239,33.896, -80.825,39.5, -87.967,34.101, -87.21,33.396, -89.0,32.043,
    -88.5,32.05, -87.063,32.289, -86.448,32.755, -86.093,32.424, -85.9,32.45, -84.95,32.55,
    -83.6,32.9, -82.2,33.4, -81.0,34.0, -79.6,34.95, -76.0,36.6, -76.0,36.0, -79.2,34.55,
    -80.3,33.55, -81.4,32.65, -82.3,32.0, -83.5,31.2, -85.0,31.0, -86.5,31.1, -88.6,31.0,
    -91.5,31.0, -91.5,36.936)),
  list(name = "Southern Coastal Plain", col = "#DCA9D4", xy = c(
    -75.5,37.666, -76.0,37.22, -76.0,39.5, -75.5,39.5, -75.5,37.666)),
  list(name = "Southern Coastal Plain", col = "#DCA9D4", xy = c(
    -91.5,27.0, -91.5,31.0, -88.6,31.0, -86.5,31.1, -85.0,31.0, -83.5,31.2, -82.3,32.0,
    -81.4,32.65, -80.3,33.55, -79.2,34.55, -78.727,34.764, -80.1,33.3, -81.05,32.1,
    -81.265,27.0, -91.5,27.0)),
  list(name = "Middle Atlantic Coastal Plain", col = "#C79BDB", xy = c(
    -75.5,27.0, -81.265,27.0, -81.05,32.1, -80.1,33.3, -78.727,34.764, -76.0,36.0,
    -76.0,37.22, -75.5,37.666, -75.5,27.0))
)

eco_defs_name <- unique(vapply(eco_defs, `[[`, "", "name"))
eco_defs_col  <- unique(vapply(eco_defs, `[[`, "", "col"))

## colours
col_geo <- "#1F4E79"   # geographic coverage bar
col_eco <- "#29ABE2"   # ecological coverage bar

## ---- basemap --------------------------------------------------------------
## state outlines (any polygon layer of US states works; this one is a small
## public GeoJSON). Replace with a local file if you prefer.
states_src <- "https://raw.githubusercontent.com/PublicaMundi/MappingAPI/master/data/geojson/us-states.json"
states <- st_read(states_src, quiet = TRUE)
states <- states[states$name %in% c("Alabama", "Georgia", "Tennessee",
                                    "North Carolina", "South Carolina",
                                    "Florida", "Mississippi", "Kentucky",
                                    "Virginia", "Arkansas", "Missouri",
                                    "Louisiana", "West Virginia", "Indiana",
                                    "Illinois", "Ohio"), ]
states <- st_make_valid(states)

ext     <- st_as_sfc(st_bbox(ext_ll, crs = 4326))
ext_big <- st_as_sfc(st_bbox(ext_ll + c(-2.8, -2.2, 2.8, 2.2), crs = 4326))
land    <- st_intersection(st_union(states), ext_big)

## ---- ecoregions -----------------------------------------------------------
## Schematic Level III ecoregion outlines (see eco_defs above). To use the real
## EPA layer instead:  eco <- st_read("us_eco_l3.shp"); eco <- st_transform(...)
ring_to_poly <- function(v) {
  m <- matrix(v, ncol = 2, byrow = TRUE)
  if (!identical(m[1, ], m[nrow(m), ])) m <- rbind(m, m[1, ])
  st_polygon(list(m))
}
eco <- st_sf(
  name     = vapply(eco_defs, `[[`, "", "name"),
  col      = vapply(eco_defs, `[[`, "", "col"),
  geometry = st_sfc(lapply(eco_defs, function(d) ring_to_poly(d$xy)), crs = 4326)
)
eco <- st_intersection(st_make_valid(eco), land)

## ---- project everything ---------------------------------------------------
p <- function(x) st_transform(x, crs_map)
frame_p  <- st_as_sfc(st_bbox(p(ext)))   # target map frame
land_p   <- p(land)
## buffer closes hairline seams between adjacent bands, then clip back to land
eco_p    <- st_intersection(st_buffer(p(eco), 9000), land_p)
states_p <- p(st_intersection(st_geometry(states), ext_big))

occ_p <- p(st_as_sf(occ, coords = c("lon", "lat"), crs = 4326))
buf_p <- st_buffer(occ_p, buffer_km * 1000)

bb <- st_bbox(frame_p)

## window that fills the panel exactly at true scale (no distortion)
fit_window <- function(bb) {
  pin <- par("pin"); w <- diff(bb[c(1, 3)]); h <- diff(bb[c(2, 4)])
  if (pin[1] / pin[2] > w / h) {
    w2 <- h * pin[1] / pin[2]
    list(x = mean(bb[c(1, 3)]) + c(-1, 1) * w2 / 2, y = bb[c(2, 4)])
  } else {
    h2 <- w * pin[2] / pin[1]
    list(x = bb[c(1, 3)], y = mean(bb[c(2, 4)]) + c(-1, 1) * h2 / 2)
  }
}

## ---- device ---------------------------------------------------------------
tiff(out_tiff, width = 12, height = 6.8, units = "in", res = dpi,
     pointsize = 13, compression = "lzw", type = "cairo")
par(family = "sans", xpd = NA)

## =============================== panel A ===================================
par(fig = c(0.030, 0.560, 0.290, 0.955), mar = c(0, 0, 1.6, 0), new = FALSE)
plot.new()
win <- fit_window(bb)
plot.window(xlim = win$x, ylim = win$y, xaxs = "i", yaxs = "i")
par(xpd = FALSE)          # clip map layers to the panel

ord <- c(which(eco_p$name == "Southeastern Plains"),
         which(eco_p$name != "Southeastern Plains"))   # background band first
for (i in ord)
  plot(st_geometry(eco_p)[i], col = eco_p$col[i], border = eco_p$col[i],
       lwd = 0.5, add = TRUE)
plot(states_p, col = NA, border = "grey25", lwd = 0.7, add = TRUE)

## buffers: hatched = sampled, dashed outline = unsampled
draw_buffer <- function(g, sampled) {
  m <- st_coordinates(g)[, 1:2]
  if (sampled) {
    polygon(m, col = "white", border = NA)
    polygon(m, density = 14, angle = 45, col = "grey10", border = "grey10",
            lwd = 0.9)
  } else {
    polygon(m, col = adjustcolor("white", 0.55), border = "grey25", lty = 2,
            lwd = 0.5)
  }
}
for (i in which(!occ$sampled)) draw_buffer(st_geometry(buf_p)[i], FALSE)
for (i in which(occ$sampled))  draw_buffer(st_geometry(buf_p)[i], TRUE)

## occurrences
xy <- st_coordinates(occ_p)
points(xy[!occ$sampled, , drop = FALSE], pch = 21, bg = "white", col = "black",
       cex = 0.95, lwd = 1)
points(xy[occ$sampled, , drop = FALSE], pch = 17, col = "black", cex = 1.35)

## state labels
st_lab <- data.frame(
  lab = c("TN", "AL", "GA", "SC", "NC", "FL", "MS"),
  lon = c(-86.9, -87.3, -83.4, -80.6, -79.5, -84.2, -89.6),
  lat = c(35.9, 33.1, 32.3, 33.9, 35.7, 30.3, 32.5))
sl <- st_coordinates(p(st_as_sf(st_lab, coords = c("lon", "lat"), crs = 4326)))
text(sl, labels = st_lab$lab, font = 2, cex = 1.0, col = "grey15")

## scale bar (km)
x0 <- bb["xmin"] + 0.05 * diff(bb[c(1, 3)]); y0 <- bb["ymin"] + 0.04 * diff(bb[c(2, 4)])
rect(x0, y0, x0 + 1e5, y0 + 9e3, col = "black", border = "black")
rect(x0 + 1e5, y0, x0 + 2e5, y0 + 9e3, col = "white", border = "black")
text(c(x0, x0 + 1e5, x0 + 2e5), y0 + 2.2e4, c("0", "100", "200"), cex = 0.7)
text(x0 + 2.2e5, y0 + 4e3, "km", cex = 0.7, adj = 0)

## north arrow
nx <- bb["xmax"] - 0.06 * diff(bb[c(1, 3)]); ny <- bb["ymin"] + 0.06 * diff(bb[c(2, 4)])
arrows(nx, ny, nx, ny + 4.5e4, length = 0.06, lwd = 1.2, col = "grey15")
text(nx, ny + 6.5e4, "N", font = 2, cex = 0.95)

par(xpd = NA)
box(col = "grey20")
title(main = bquote(bold("A") ~ ~ italic("Quercus georgiana") ~ "occurrences and" ~
                      italic("ex situ") ~ "samples"),
      adj = 0, cex.main = 1.15, font.main = 1, line = 0.4)

## =============================== legends ===================================
par(fig = c(0.030, 0.620, 0.005, 0.275), mar = c(0, 0, 0, 0), new = TRUE)
plot.new(); plot.window(c(0, 1), c(0, 1))

legend(0, 1.04, bty = "n", cex = 0.95, pt.cex = c(1.4, 1.05), y.intersp = 1.3,
       pch = c(17, 21), pt.bg = c(NA, "white"), col = "black",
       legend = c(expression(paste("Occurrence sampled for ", italic("ex situ"),
                                   " collection")),
                  "Unsampled occurrence"))
sw <- function(y, hatched) {
  xs <- c(0.004, 0.030, 0.030, 0.004); ys <- c(y - 0.045, y - 0.045, y + 0.045, y + 0.045)
  if (hatched) {
    polygon(xs, ys, col = "white", border = "grey10")
    polygon(xs, ys, density = 22, angle = 45, col = "grey10", border = "grey10")
  } else {
    polygon(xs, ys, col = "white", border = "grey25", lty = 2)
  }
}
sw(0.42, TRUE);  text(0.045, 0.42, "Buffer of a sampled occurrence",
                      adj = c(0, 0.5), cex = 0.95)
sw(0.16, FALSE); text(0.045, 0.16, "Buffer of an unsampled occurrence",
                      adj = c(0, 0.5), cex = 0.95)
legend(0.60, 1.04, bty = "n", cex = 0.9, y.intersp = 1.1, title = "Ecoregions",
       title.adj = 0, fill = eco_defs_col, border = NA, legend = eco_defs_name)

## =============================== panel B ===================================
par(fig = c(0.600, 0.880, 0.330, 0.930), mar = c(4.0, 4.2, 2.2, 0.2), new = TRUE)
plot.new()
plot.window(xlim = c(0.4, 2.6), ylim = c(0, 105), xaxs = "i", yaxs = "i")

par(xpd = FALSE)
abline(h = seq(0, 100, 25), col = "grey88", lwd = 0.6)
w <- 0.30
for (k in 1:2) {
  v  <- c(cov_geo, cov_eco)[k]
  cc <- c(col_geo, col_eco)[k]
  rect(k - w, 0, k + w, 100, col = "grey95", border = "grey60", lty = 2, lwd = 0.7)
  rect(k - w, 0, k + w, v, col = cc, border = NA)
  text(k, v + 3, paste0(round(v), "%"), font = 2, cex = 1.05, col = cc, adj = c(0.5, 0))
  text(k, v + (100 - v) / 2, "gap", font = 3, cex = 0.95, col = "grey45")
}
par(xpd = NA)
axis(2, at = seq(0, 100, 25), labels = c("0", "25%", "50%", "75%", "100%"),
     las = 1, cex.axis = 0.95, tcl = -0.25, mgp = c(2, 0.5, 0))
axis(1, at = 1:2, labels = FALSE, tcl = 0)
text(1:2, -8, c("Geographic\ncoverage", "Ecological\ncoverage"), cex = 1.0, adj = c(0.5, 1))
mtext(bquote("Estimated diversity represented" ~ italic("ex situ")), side = 2,
      line = 2.8, cex = 1.0)
title(main = bquote(bold("B") ~ ~ "Coverage and gaps"), adj = 0, cex.main = 1.15,
      font.main = 1, line = 0.8)

## ---- footnotes ------------------------------------------------------------
par(fig = c(0, 1, 0, 1), mar = c(0, 0, 0, 0), new = TRUE)
plot.new(); plot.window(c(0, 1), c(0, 1))
text(0.600, 0.235, "Geographic coverage = hatched buffer area / total buffer area",
     adj = c(0, 1), cex = 0.85, col = col_geo)
text(0.600, 0.170, paste("Ecological coverage = ecoregions in hatched buffers /",
                         "\n           ecoregions in all buffers"),
     adj = c(0, 1), cex = 0.85, col = "#0E7FAE")

dev.off()
