# ==============================================================================
# APPENDIX S1: collection points, elevations and map for the 2021 experiment.
# Source: data/raw_data/WC Garden Collections_Coords.xlsx, sheet "Sheet1" (Macrosystems 2021 Hybrid
# Garden Collection Record): six populations, one GPS point per cross type collected. These are the
# six population codes in the data (BR, IC, JC, LC, OC = Oak Creek, SF). Weber River (sheet
# "Sources for GH Expt") was collected in 2023 for a later experiment and is not used.
# Coordinates are degree-decimal-minute strings (N38° 04.184') converted to decimal degrees here;
# the workbook's decimal Lat/Lon columns (on other sheets) are misentered and not used.
# Elevations: elevatr::get_elev_point (USGS Elevation Point Query Service; network access needed).
# Map: elevatr terrain tiles, Natural Earth rivers and USGS NHD flowlines (network access needed).
# ==============================================================================
# Run from the project root (R 4.6.1):  Rscript r_code/12_map_appendix_S1.R
# Writes: results/appendices/appendix_S1_collection_points.csv, appendix_S1_elevation_ranges.csv,
#         figures/AppendixS1_map.pdf/.png
# ==============================================================================
suppressPackageStartupMessages({
  library(readxl); library(tidyverse); library(sf); library(terra); library(elevatr)
  library(ggspatial); library(rnaturalearth); library(nhdplusTools); library(ggrepel); library(patchwork) })
source("r_code/helpers/utilities_data_prep.R")
out <- "results/appendices"; dir.create(out, showWarnings = FALSE, recursive = TRUE)

# ---- Sheet1 -> one row per collection point ----
raw <- read_excel("data/raw_data/WC Garden Collections_Coords.xlsx", sheet = "Sheet1", col_names = FALSE, .name_repair = "minimal")
names(raw)[1:5] <- c("population", "dates", "cross_type", "n_genotypes", "gps")
pts <- raw[-(1:2), 1:5] %>% fill(population, dates) %>% filter(!is.na(gps)) %>%
  mutate(n_genotypes = as.integer(n_genotypes))
dm <- function(s, hemi) {          # "N38° 04.184'" -> 38.06973; west longitudes negative
  m <- str_match(s, paste0(hemi, "\\s*(\\d+)°\\s*([0-9.]+)'"))
  sgn <- ifelse(hemi %in% c("S", "W"), -1, 1); sgn * (as.numeric(m[, 2]) + as.numeric(m[, 3]) / 60) }
pts <- pts %>% mutate(lat = dm(gps, "N"), lon = dm(gps, "W"),
                      site = str_trim(sub(",?\\s*\\(.*$", "", population)),
                      state = str_match(population, ",\\s*([A-Z]{2})\\)")[, 2],
                      code = recode(site, "Blue River" = "BR", "Indian Creek" = "IC", "Jack's Canyon" = "JC",
                                    "Little Colorado River" = "LC", "Oak Creek" = "OC", "San Francisco River" = "SF"))
stopifnot(nrow(pts) == 18, !anyNA(c(pts$lat, pts$lon, pts$state)), abs(dm("N38° 04.184'", "N") - 38.06973) < 1e-5)
ds <- load_and_prepare_data()
stopifnot(setequal(unique(pts$code), unique(ds$pop)))
cat("Sheet1 populations match the six population codes in the data:", paste(sort(unique(pts$code)), collapse = ", "), "\n")
# every mother-tree ID in the data is within the number of genotypes collected for that population x cross type
idchk <- ds %>% group_by(pop, type) %>% summarise(max_geno = max(geno), .groups = "drop") %>%
  left_join(pts %>% transmute(pop = code, type = substr(cross_type, 1, 1), n_genotypes), by = c("pop", "type"))
cat("Mother-tree numbers above the Sheet1 genotype count:", sum(idchk$max_geno > idchk$n_genotypes), "\n")

# ---- elevations ----
el <- get_elev_point(st_as_sf(pts, coords = c("lon", "lat"), crs = 4326, remove = FALSE), src = "epqs")
for (k in 1:5) {   # the service sometimes times out (HTTP 504) for single points: query those again
  miss <- is.na(el$elevation); if (!any(miss)) break
  Sys.sleep(15)
  el$elevation[miss] <- get_elev_point(st_as_sf(pts[miss, ], coords = c("lon", "lat"), crs = 4326, remove = FALSE), src = "epqs")$elevation }
pts$elevation_m <- round(el$elevation, 0)
stopifnot(!anyNA(pts$elevation_m))
cp <- pts %>% transmute(Population = site, State = state, Code = code, `Cross type collected` = cross_type, `Genotypes collected` = n_genotypes,
                        `Collection dates` = dates, `GPS (as recorded)` = gps, Latitude = round(lat, 5), Longitude = round(lon, 5), `Elevation (m)` = elevation_m)
readr::write_excel_csv(cp, file.path(out, "appendix_S1_collection_points.csv"))
rng <- pts %>% group_by(Population = site, State = state, Code = code) %>%
  summarise(`Elevation min (m)` = min(elevation_m), `Elevation max (m)` = max(elevation_m), .groups = "drop")
readr::write_excel_csv(rng, file.path(out, "appendix_S1_elevation_ranges.csv"))
cat("\nCollection points:\n"); print(as.data.frame(cp %>% select(Code, `Cross type collected`, Latitude, Longitude, `Elevation (m)`)), row.names = FALSE)
cat("\nElevation range by population (m):\n"); print(as.data.frame(rng), row.names = FALSE)
cat(sprintf("Overall: %d-%d m (collection points)\n", min(pts$elevation_m), max(pts$elevation_m)))

# ---- map (no title) ----
CM <- 1 / 2.54; FONT <- "Arial"
pad <- 0.35
xlim <- range(pts$lon) + c(-pad, pad); ylim <- range(pts$lat) + c(-pad, pad)
w_km <- as.numeric(st_distance(st_sfc(st_point(c(xlim[1], mean(ylim))), crs = 4326), st_sfc(st_point(c(xlim[2], mean(ylim))), crs = 4326))) / 1000
h_km <- as.numeric(st_distance(st_sfc(st_point(c(mean(xlim), ylim[1])), crs = 4326), st_sfc(st_point(c(mean(xlim), ylim[2])), crs = 4326))) / 1000
bbox_sf <- st_as_sf(data.frame(x = xlim, y = ylim), coords = c("x", "y"), crs = 4326)
dem <- get_elev_raster(locations = bbox_sf, z = 8, clip = "bbox") %>% rast() %>% project("EPSG:4326")
hill <- shade(terrain(dem, v = "slope", unit = "radians"), terrain(dem, v = "aspect", unit = "radians"), angle = 40, direction = 315)
ev <- values(dem)[, 1]; hv <- values(hill)[, 1]; xy <- crds(dem, na.rm = FALSE); keep <- !is.na(ev) & !is.na(hv)
erng <- quantile(ev[keep], c(0.005, 0.995), names = FALSE)
pal <- colorRampPalette(c("#EDE3C8", "#D9C398", "#B98F58", "#8B5E34", "#6B4226"))
cols <- pal(256)[pmin(256, pmax(1, as.integer((ev[keep] - erng[1]) / diff(erng) * 255) + 1))]
rgbm <- sweep(col2rgb(cols) / 255, 2, 0.45 + 0.55 * hv[keep], `*`); rgbm[rgbm > 1] <- 1
base <- data.frame(x = xy[keep, 1], y = xy[keep, 2], fill = rgb(rgbm[1, ], rgbm[2, ], rgbm[3, ]))
sf_use_s2(FALSE)
states <- st_as_sf(maps::map("state", plot = FALSE, fill = TRUE)) %>% st_make_valid()
sf_use_s2(TRUE)
ne_riv <- bind_rows(ne_download(scale = 10, type = "rivers_lake_centerlines", category = "physical", returnclass = "sf"),
                    ne_download(scale = 10, type = "rivers_north_america", category = "physical", returnclass = "sf")) %>%
  transmute(lw = pmax(0.15, (7 - coalesce(scalerank, 6)) * 0.08))
nhd <- map(split(pts, pts$code), function(p) {
  aoi <- st_as_sfc(st_bbox(c(xmin = min(p$lon), xmax = max(p$lon), ymin = min(p$lat), ymax = max(p$lat)), crs = st_crs(4326))) %>%
    st_transform(5070) %>% st_buffer(8000) %>% st_transform(4326)
  # stop rather than draw the map without rivers: a silent gap would make the figure differ from the committed one
  nhd_fail <- function(why) stop("USGS NHD flowline download failed for population ", p$code[1], " (", why, "). ",
                                 "Rerun when the service is available, or run 'Rscript run_all.R --skip-map' to keep the committed map.", call. = FALSE)
  fl <- tryCatch(get_nhdplus(AOI = aoi, realization = "flowline"), error = function(e) nhd_fail(conditionMessage(e)))
  if (is.null(fl) || !nrow(fl)) nhd_fail("no flowlines returned")
  fl %>% st_transform(4326) %>% st_geometry() %>% st_sf() %>% mutate(lw = 0.2) })
rivers <- bind_rows(ne_riv, bind_rows(nhd))
site_cols <- setNames(RColorBrewer::brewer.pal(6, "Dark2"), sort(unique(pts$site)))
# labels placed by hand (lx, ly) clear of the points, scale bar and insets, with a leader line to the site centroid
lab <- pts %>% group_by(site) %>% summarise(lon = mean(lon), lat = mean(lat), .groups = "drop") %>%
  left_join(tribble(~site, ~lx, ~ly,
                    "Indian Creek", -109.57, 37.78, "Oak Creek", -111.85, 35.45, "Jack's Canyon", -110.55, 35.50,
                    "Little Colorado River", -109.30, 34.62, "San Francisco River", -109.20, 34.08, "Blue River", -109.80, 33.55), by = "site")
stopifnot(!anyNA(lab$lx))
ct_shape <- c(Fremont = 21, Hybrid = 22, Narrowleaf = 24)
mp <- ggplot() +
  geom_raster(data = base, aes(x, y), fill = base$fill) +
  geom_sf(data = rivers, aes(linewidth = lw), colour = "#2E6F95", alpha = 0.8, inherit.aes = FALSE, show.legend = FALSE) + scale_linewidth_identity() +
  geom_sf(data = states, fill = NA, colour = "grey20", linewidth = 0.35, inherit.aes = FALSE) +
  geom_segment(data = lab, aes(x = lx, y = ly, xend = lon, yend = lat), colour = "grey30", linewidth = 0.3) +
  geom_point(data = pts, aes(lon, lat, fill = site, shape = cross_type), colour = "black", size = 1.8, stroke = 0.4) +
  geom_label(data = lab, aes(lx, ly, label = site, colour = site), fill = "white", alpha = 0.9, fontface = "bold", size = 2.4,
             family = FONT, label.size = 0, label.padding = unit(0.12, "lines"), show.legend = FALSE) +
  scale_fill_manual(values = site_cols, guide = "none") + scale_colour_manual(values = site_cols, guide = "none") +
  scale_shape_manual(values = ct_shape, name = "Cross type collected") +
  annotate("text", x = c(-111.6, -108.75, -111.2), y = c(37.4, 35.6, 36.2), label = c("UT", "NM", "AZ"), colour = "grey25",
           fontface = "italic", size = 3, family = FONT) +
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE, crs = 4326) +
  annotation_scale(location = "br", width_hint = 0.3, text_cex = 0.6, height = unit(0.15, "cm"), text_family = FONT) +
  annotation_north_arrow(location = "tr", which_north = "true", style = north_arrow_orienteering(text_size = 6, text_family = FONT),
                         height = unit(0.8, "cm"), width = unit(0.8, "cm")) +
  labs(x = NULL, y = NULL) + theme_bw(base_size = 7, base_family = FONT) +
  theme(legend.position = "bottom", legend.key.size = unit(0.3, "cm"), panel.grid = element_blank(), axis.text = element_text(size = 6, colour = "black"),
        plot.background = element_rect(fill = "white", colour = NA), plot.margin = margin(3, 3, 3, 3))
leg <- ggplot(data.frame(x = 1, y = seq(erng[1], erng[2], length.out = 256)), aes(x, y, fill = y)) + geom_tile() +
  scale_fill_gradientn(colours = pal(256), guide = "none") + scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(position = "right", expand = c(0, 0), breaks = seq(1000, 3500, by = 500), labels = function(b) paste0(format(b, big.mark = ","), " m")) +
  ggtitle("Elevation") + theme_void(base_family = FONT) +
  theme(axis.text.y = element_text(size = 5.5, margin = margin(l = 2)), plot.title = element_text(size = 6, face = "bold", margin = margin(b = 5)),
        plot.background = element_rect(fill = "white", colour = "grey60", linewidth = 0.2), plot.margin = margin(3, 5, 3, 5))
us <- map_data("state")
loc <- ggplot() + geom_polygon(data = us, aes(long, lat, group = group), fill = "grey88", colour = "white", linewidth = 0.1) +
  annotate("rect", xmin = xlim[1], xmax = xlim[2], ymin = ylim[1], ymax = ylim[2], fill = NA, colour = "red", linewidth = 0.4) +
  coord_fixed(1.3) + theme_void() + theme(plot.background = element_rect(fill = "white", colour = "grey60", linewidth = 0.2))
final <- mp + inset_element(loc, left = 0.02, bottom = 0.86, right = 0.42, top = 0.995) +
  inset_element(leg, left = 0.02, bottom = 0.03, right = 0.2, top = 0.24)
W <- 8.9 * CM; H <- min(23 * CM, (W - 1.3 * CM) * h_km / w_km + 1.3 * CM)   # panel width excludes the latitude labels
ggsave("figures/AppendixS1_map.pdf", final, width = W, height = H, units = "in", device = cairo_pdf, bg = "white")
ragg::agg_png("figures/AppendixS1_map.png", width = W, height = H, units = "in", res = 600, background = "white"); plot(final); invisible(dev.off())
cat(sprintf("wrote figures/AppendixS1_map (.pdf, .png): %.1f x %.1f cm\n", W * 2.54, H * 2.54))
