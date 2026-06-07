library(rstac)
library(terra)
library(sf)

print("=== [DAPUR DATA CORE - V2] EKSTRAKSI AUTONOMOUS DENGAN FALLBACK LOOP ===")

if(!file.exists("batas_desa_pasanggrahan.geojson")) {
  stop("Eror Mutlak: Berkas 'batas_desa_pasanggrahan.geojson' wajib ada di folder proyek!")
}
batas_desa_sf <- st_read("batas_desa_pasanggrahan.geojson", quiet = TRUE)
bbox_query    <- unname(as.numeric(st_bbox(batas_desa_sf)))

setGDALconfig("GDAL_HTTP_MAX_RETRY", "5")     
setGDALconfig("GDAL_HTTP_RETRY_DELAY", "3")   
setGDALconfig("GDAL_HTTP_VERSION", "1.1")     

s_api <- stac("https://planetarycomputer.microsoft.com/api/stac/v1")
list_fusion_historis <- list()

tahun_sekarang      <- as.numeric(format(Sys.Date(), "%Y"))
tahun_akhir_histori <- tahun_sekarang - 1
tahun_loop           <- 2021:tahun_akhir_histori

for (thn in tahun_loop) {
  print(paste("------------------------------------------------------------"))
  print(paste("--> Memproses Siklus Fenologi Juni Tahun:", thn))
  
  # 1. Kueri Metadata Satelit
  pencarian_s2 <- s_api %>%
    stac_search(collections = "sentinel-2-l2a", bbox = bbox_query,
                datetime = paste0(thn, "-06-01/", thn, "-06-30"), limit = 10) %>%
    get_request()
  
  pencarian_s1 <- s_api %>%
    stac_search(collections = "sentinel-1-rtc", bbox = bbox_query,
                datetime = paste0(thn, "-06-01/", thn, "-06-30"), limit = 5) %>%
    get_request()
  
  if(length(pencarian_s2$features) == 0 || length(pencarian_s1$features) == 0) {
    print(paste("      [Lewat] Data tahun", thn, "tidak lengkap di server."))
    next
  }
  
  # 2. URUTKAN INDEKS AWAN DARI TERKECIL KE TERBESAR
  semua_awan   <- sapply(pencarian_s2$features, function(x) x$properties$`eo:cloud_cover`)
  indeks_urut  <- order(semua_awan)
  
  s2_terpilih_berhasil <- FALSE
  
  # INTERNAL FALLBACK LOOP: Menguji ubin satu per satu hingga menemukan yang overlap
  for (idx in indeks_urut) {
    tryCatch({
      s2_signed   <- items_sign(pencarian_s2$features[[idx]], sign_planetary_computer())
      crs_satelit <- crs(rast(s2_signed$assets$B04$href))
      poly_desa   <- project(vect(batas_desa_sf), crs_satelit)
      
      # Tes Pemotongan Optik
      r_red <- crop(rast(paste0("/vsicurl/", s2_signed$assets$B04$href)), poly_desa)
      r_nir <- crop(rast(paste0("/vsicurl/", s2_signed$assets$B08$href)), poly_desa)
      ndvi  <- (r_nir - r_red) / (r_nir + r_red)
      
      # Jika sampai baris ini tidak eror, berarti ubin ini overlap sempurna dengan desa!
      s2_terpilih_berhasil <- TRUE
      break
    }, error = function(e) {
      # Jika meleset, abaikan erornya dan biarkan loop berlanjut mencari ubin cadangan berikutnya
    })
  }
  
  if (!s2_terpilih_berhasil) {
    print(paste("      [GAGAL] Semua ubin S2 tahun", thn, "meleset dari poligon desa."))
    next
  }
  
  # 3. PROSES DATA RADAR SENTINEL-1 DENGAN FALLBACK LOOP SEJENIS
  s1_terpilih_berhasil <- FALSE
  for (idx_s1 in 1:length(pencarian_s1$features)) {
    tryCatch({
      s1_item   <- pencarian_s1$features[[idx_s1]]
      s1_signed <- items_sign(s1_item, sign_planetary_computer())
      
      r_vv <- crop(rast(paste0("/vsicurl/", s1_signed$assets$vv$href)), poly_desa)
      r_vh <- crop(rast(paste0("/vsicurl/", s1_signed$assets$vh$href)), poly_desa)
      
      r_vv_res <- resample(r_vv, ndvi, method = "bilinear")
      r_vh_res <- resample(r_vh, ndvi, method = "bilinear")
      radar_ratio <- r_vh_res / r_vv_res
      
      sawah_mask <- ndvi >= 0.2 & ndvi <= 0.85
      ndvi_clean <- ndvi
      ndvi_clean[sawah_mask == 0] <- (radar_ratio[sawah_mask == 0] * 0.5) + 0.2 
      
      list_fusion_historis[[as.character(thn)]] <- mask(ndvi_clean, poly_desa)
      s1_terpilih_berhasil <- TRUE
      break
    }, error = function(e) {
      # Coba orbit lintasan radar berikutnya jika lintasan pertama meleset
    })
  }
  
  if (s1_terpilih_berhasil) {
    print(paste("      [SUKSES] Data Fusion Optik+Radar tahun", thn, "berhasil dirakit."))
  } else {
    print(paste("      [Lewat] Data Radar tahun", thn, "tidak ada yang overlap."))
  }
}

print("------------------------------------------------------------")
print("--> Sinkronisasi Resolusi Matriks Spasial Global...")
raster_acuan <- list_fusion_historis[[length(list_fusion_historis)]]
for (nama_thn in names(list_fusion_historis)) {
  list_fusion_historis[[nama_thn]] <- resample(list_fusion_historis[[nama_thn]], raster_acuan, method = "near")
}

tumpukan_raster    <- rast(list_fusion_historis)
mean_historis_june <- app(tumpukan_raster, fun = "mean", na.rm = TRUE)
sd_historis_june   <- app(tumpukan_raster, fun = "sd", na.rm = TRUE)

nama_file_save <- paste0("baseline_2021_to_", tahun_akhir_histori, ".rds")
saveRDS(list(mean_raster = terra::wrap(mean_historis_june), sd_raster = terra::wrap(sd_historis_june)), nama_file_save)

print(paste("=== KESUKSESAN MUTLAK:", nama_file_save, "BERHASIL DIKUNCI MATANG! ==="))