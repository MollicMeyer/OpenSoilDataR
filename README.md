# OpenSoilDataR Tutorial 🚜🌱

#### by Meyer Bohn, Geospatial Laboratory for Soil Informatics --- Department of Agronomy, Iowa State University

## 📚 Citation

If you use **OpenSoilDataR** in your work, please cite it as:

> Bohn, M.P. 2025. *OpenSoilDataR: Tools for querying, harmonizing, and summarizing public soil datasets*. GitHub repository, https://github.com/MollicMeyer/OpenSoilDataR

This repository provides a tutorial on using **OpenSoilDataR** to fetch soil property rasters from:
- **POLARIS (PSP)** ```Chaney, N.W., Minasny, B., Herman, J.D., Nauman, T.W., et al. 2019. POLARIS Soil Properties: 30-m Probabilistic Maps of Soil Properties Over the Contiguous United States. Water Resour Res 55, 2916–2938. https://doi.org/10.1029/2018WR022797```
- **SoilGrids v2 (SG2)** ```Hengl, T., De Jesus, J.M., Heuvelink, G.B.M., Gonzalez, et al. 2017. SoilGrids250m: Global gridded soil information based on machine learning. PLoS One 12, e0169748. https://doi.org/10.1371/JOURNAL.PONE.0169748```
- **SOLUS100 (SOL)**  ```Nauman, T.W., Kienast-Brown, S., Roecker, S.M., Brungard, et al. 2024. Soil landscapes of the United States (SOLUS): Developing predictive soil property maps of the conterminous United States using hybrid training sets. Soil Science Society of America Journal 88, 2046–2065. https://doi.org/10.1002/SAJ2.20769```
- **Cal Soil Res Lab Soil Props 800m  (CSRL)**  **⚠️ Note:** This dataset is intended for **statewide and regional soil analysis**.
For **local-scale, fine-resolution soil data**, refer to **SSURGO (Soil Survey Geographic Database)**.```Walkinshaw, Mike, A.T. O'Geen, D.E. Beaudette. "Soil Properties." California Soil Resource Lab, 1 Oct. 2020,
casoilresource.lawr.ucdavis.edu/soil-properties/.```
- **SSURGO/gNATSGO/RSS (SGO)**  ```USDA Natural Resources Conservation Service, 2024. Soil Survey Geographic (SSURGO) Database. United States Department of Agriculture. Available at: https://sdmdataaccess.nrcs.usda.gov```

- compute **zonal statistics** for soil properties
- create **SoilProfileColletion** objects from **points**,**zones**, and **rasters**
- plot comparison soil depth functions for different public soil dataset map products

---

## 📥 Installation

Ensure you have the necessary **R packages** installed:

```r
install.packages(c("terra", "sf", "dplyr", "tidyr", "rgeedim", "soilDB", "httr"))
devtools::install_github("MollicMeyer/OpenSoilDataR")  

````
##### Load the libraries:
```r
library(rgeedim)
library(OpenSoilDataR)
library(terra)
library(soilDB)
library(sf)

````

🔑 Authenticate & Initialize Google Earth Engine
Before fetching data, authenticate Google Earth Engine:

```r
gd_authenticate(auth_mode = "notebook")  # Authenticate
gd_initialize()  # Initialize

````

📍 Define Area of Interest (AOI)
Set your area of interest (AOI) using a raster file:

```r
aoi_path <- "T:/IA-Kitchen/R/sabR/SABRplots.shp"
aoi <- terra::vect(aoi_path)
```

---

## 🌍 Fetch Soil Property Data

### PSP (POLARIS)

```r
psp <- fetch_PSP(
  aoi = aoi,
  properties = c("sand_mean", "clay_mean", "om_mean", "ph_mean"),
  depths = c("0_5", "5_15", "15_30", "30_60", "60_100", "100_200"),
  crs = "EPSG:4326",
  scale = 30,
  tosoc = TRUE,
  convertOM = TRUE,
  export = FALSE
)
```

### SG2 (SoilGrids v2)

```r
sg2 <- fetch_SG2(
  aoi = aoi,
  properties = c("sand", "clay", "soc", "phh2o"),
  depths = c("0-5cm", "5-15cm", "15-30cm", "30-60cm", "60-100cm", "100-200cm"),
  crs = "EPSG:4326",
  scale = 250,
  export = FALSE
)
```

### SOLUS100

```r
sol <- fetch_SOL(
  aoi = aoi,
  properties = c("sandtotal", "claytotal", "soc", "ph1to1h2o"),
  depths = c("0_cm", "5_cm", "15_cm", "30_cm", "60_cm", "100_cm", "150_cm"),
  measures = c("p"),
  export = FALSE
)
```

### SGO (SSURGO)

```r
sgo <- fetch_SGO(
  aoi = aoi,
  properties = c("sandtotal_r", "claytotal_r", "om_r", "ph1to1h2o_r"),
  depths = list(c(0, 5), c(5, 15), c(15, 30), c(30, 60), c(60, 100), c(100, 200)),
  crs = "EPSG:4326",
  res = 30,
  export = FALSE,
  tosoc = TRUE
)
```

---

## 📊 Zonal Statistics

```r
plots <- c("NE_plot", "NW_plot", "SE_plot", "SW_plot")
std_props <- c("sand", "clay", "soc", "ph")
tdepth <- 0
bdepth <- 30

z_psp <- s.zonalstats(psp$stack, tdepth, bdepth, std_props, aoi, plots, c("mean", "sd"), wtd.mean = TRUE)
z_sg2 <- s.zonalstats(sg2$stack, tdepth, bdepth, std_props, aoi, plots, c("mean", "sd"), wtd.mean = TRUE)
z_sol <- s.zonalstats(sol$stack, tdepth, bdepth, std_props, aoi, plots, c("mean", "sd"), wtd.mean = TRUE)
z_sgo <- s.zonalstats(sgo$stack, tdepth, bdepth, std_props, aoi, plots, c("mean", "sd"), wtd.mean = TRUE)
```

Combine and visualize:

```r
z_all <- bind_rows(z_psp$Weighted, z_sg2$Weighted, z_sol$Weighted, z_sgo$Weighted)
z_wide <- z_all %>% pivot_wider(names_from = Statistic, values_from = ZonalStats)

ggplot(z_wide, aes(x = Plot, y = mean, color = Product)) +
  geom_point(position = position_dodge(0.6), size = 4) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.5, position = position_dodge(0.6)) +
  facet_wrap(~SoilProperty, scales = "free_y") +
  theme_minimal(base_size = 13) +
  labs(x = "Plot", y = "Mean ± SD", title = "Zonal Statistics by Plot and Product", color = "Product") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
```

---

## 💼 SoilProfileCollection (SPC) Conversions

### From Raster Stack (Global)

```r
sgo_spc <- ras_to_SPC(sgo$stack, source = "SGO")
psp_spc <- ras_to_SPC(psp$stack, source = "PSP")
sol_spc <- ras_to_SPC(sol$stack, source = "SOL")
sg2_spc <- ras_to_SPC(sg2$stack, source = "SG2")

spc_list <- list(sgo_spc, psp_spc, sg2_spc, sol_spc)
plot_depth_functions(spc_list, c("SGO", "PSP", "SG2", "SOL"), c("sand", "clay", "soc", "ph"), c(0,5,15,30,60,100,150))
```

### From Polygons (Zone-Averaged)

```r
subset_ids <- c("NW_plot", "SW_plot", "NE_plot", "SE_plot")

sgo_spc <- zones_to_SPC(sgo$stack, aoi, stat = "mean", id_column = "Name", subset_ids = subset_ids)
psp_spc <- zones_to_SPC(psp$stack, aoi, stat = "mean", id_column = "Name", subset_ids = subset_ids)
sol_spc <- zones_to_SPC(sol$stack, aoi, stat = "mean", id_column = "Name", subset_ids = subset_ids)
sg2_spc <- zones_to_SPC(sg2$stack, aoi, stat = "mean", id_column = "Name", subset_ids = subset_ids)

spc_list <- list(sgo_spc, psp_spc, sg2_spc, sol_spc)
plot_depth_functions(spc_list, c("SGO", "PSP", "SG2", "SOL"), c("sand", "clay", "soc", "ph"), c(0,5,15,30,60,100,150))
```

### From Random Points

```r
set.seed(42)
aoi_sf <- if (!inherits(aoi, "sf")) sf::st_as_sf(aoi) else aoi
rand_pts <- sf::st_sample(aoi_sf, size = 25, type = "random") %>% sf::st_as_sf()
rand_pts$Name <- paste0("pt_", seq_len(nrow(rand_pts)))

sgo_spc <- pts_to_SPC(sgo$stack, rand_pts, id_column = "Name")
psp_spc <- pts_to_SPC(psp$stack, rand_pts, id_column = "Name")
sol_spc <- pts_to_SPC(sol$stack, rand_pts, id_column = "Name")
sg2_spc <- pts_to_SPC(sg2$stack, rand_pts, id_column = "Name")

spc_list <- list(sgo_spc, psp_spc, sg2_spc, sol_spc)
plot_depth_functions(spc_list, c("SGO", "PSP", "SG2", "SOL"), c("sand", "clay", "soc", "ph"), c(0,5,15,30,60,100,150))
```

---
