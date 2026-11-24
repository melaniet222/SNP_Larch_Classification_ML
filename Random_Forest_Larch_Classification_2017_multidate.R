# ============================================================
# LARCH CLASSIFICATION - SENTINEL-2 - SNP - 2017
# MULTI-DATE + ASPECT-IMPROVED MODEL
# Dates: 04 Oct, 14 Oct, 19 Oct 2017
# ============================================================


# ------------------------------------------------------------
# 1. Load packages
# ------------------------------------------------------------

library(terra)
library(sf)
library(ranger)


# ------------------------------------------------------------
# 2. Define paths
# ------------------------------------------------------------

raster_dir <- "S:/course/geo715/student/phenology/gis_work/RastersGIS"

gdb_file <- "S:/course/geo715/student/phenology/gis_work/Phenology.gdb"

output_dir <- "S:/course/geo715/student/phenology/gis_work/RastersGIS"


# ------------------------------------------------------------
# 3. Load reference raster
# ------------------------------------------------------------

template <- rast(
  file.path(
    raster_dir,
    "S2_20171014_10m_B08_SNP_Forest.tif"
  )
)


# ------------------------------------------------------------
# 4. Function to load and align one raster
# ------------------------------------------------------------

load_align <- function(filename) {
  
  r <- rast(
    file.path(
      raster_dir,
      filename
    )
  )
  
  r_aligned <- resample(
    r,
    template,
    method = "bilinear"
  )
  
  return(r_aligned)
}


# ------------------------------------------------------------
# 5. Load 04 October
# ------------------------------------------------------------

b02_04 <- load_align("S2_20171004_10m_B02_SNP_Forest.tif")
b03_04 <- load_align("S2_20171004_10m_B03_SNP_Forest.tif")
b04_04 <- load_align("S2_20171004_10m_B04_SNP_Forest.tif")
b05_04 <- load_align("S2_20171004_10m_B05_SNP_Forest.tif")
b06_04 <- load_align("S2_20171004_10m_B06_SNP_Forest.tif")
b07_04 <- load_align("S2_20171004_10m_B07_SNP_Forest.tif")
b08_04 <- load_align("S2_20171004_10m_B08_SNP_Forest.tif")
b8a_04 <- load_align("S2_20171004_10m_B8A_SNP_Forest.tif")
b11_04 <- load_align("S2_20171004_10m_B11_SNP_Forest.tif")
b12_04 <- load_align("S2_20171004_10m_B12_SNP_Forest.tif")


# ------------------------------------------------------------
# 6. Load 14 October
# ------------------------------------------------------------

b02_14 <- load_align("S2_20171014_10m_B02_SNP_Forest.tif")
b03_14 <- load_align("S2_20171014_10m_B03_SNP_Forest.tif")
b04_14 <- load_align("S2_20171014_10m_B04_SNP_Forest.tif")
b05_14 <- load_align("S2_20171014_10m_B05_SNP_Forest.tif")
b06_14 <- load_align("S2_20171014_10m_B06_SNP_Forest.tif")
b07_14 <- load_align("S2_20171014_10m_B07_SNP_Forest.tif")

b08_14 <- template

b8a_14 <- load_align("S2_20171014_10m_B8A_SNP_Forest.tif")
b11_14 <- load_align("S2_20171014_10m_B11_SNP_Forest.tif")
b12_14 <- load_align("S2_20171014_10m_B12_SNP_Forest.tif")


# ------------------------------------------------------------
# 7. Load 19 October
# ------------------------------------------------------------

b02_19 <- load_align("S2_20171019_10m_B02_SNP_Forest.tif")
b03_19 <- load_align("S2_20171019_10m_B03_SNP_Forest.tif")
b04_19 <- load_align("S2_20171019_10m_B04_SNP_Forest.tif")
b05_19 <- load_align("S2_20171019_10m_B05_SNP_Forest.tif")
b06_19 <- load_align("S2_20171019_10m_B06_SNP_Forest.tif")
b07_19 <- load_align("S2_20171019_10m_B07_SNP_Forest.tif")
b08_19 <- load_align("S2_20171019_10m_B08_SNP_Forest.tif")
b8a_19 <- load_align("S2_20171019_10m_B8A_SNP_Forest.tif")
b11_19 <- load_align("S2_20171019_10m_B11_SNP_Forest.tif")
b12_19 <- load_align("S2_20171019_10m_B12_SNP_Forest.tif")


# ------------------------------------------------------------
# 8. Combine all 30 predictors
# ------------------------------------------------------------

sentinel_multidate <- c(
  b02_04, b03_04, b04_04, b05_04, b06_04,
  b07_04, b08_04, b8a_04, b11_04, b12_04,
  
  b02_14, b03_14, b04_14, b05_14, b06_14,
  b07_14, b08_14, b8a_14, b11_14, b12_14,
  
  b02_19, b03_19, b04_19, b05_19, b06_19,
  b07_19, b08_19, b8a_19, b11_19, b12_19
)


# ------------------------------------------------------------
# 9. Name predictors
# ------------------------------------------------------------

names(sentinel_multidate) <- c(
  "B02_04", "B03_04", "B04_04", "B05_04", "B06_04",
  "B07_04", "B08_04", "B8A_04", "B11_04", "B12_04",
  
  "B02_14", "B03_14", "B04_14", "B05_14", "B06_14",
  "B07_14", "B08_14", "B8A_14", "B11_14", "B12_14",
  
  "B02_19", "B03_19", "B04_19", "B05_19", "B06_19",
  "B07_19", "B08_19", "B8A_19", "B11_19", "B12_19"
)


# ------------------------------------------------------------
# 10. Check stack
# ------------------------------------------------------------

cat("\nNumber of predictors:\n")
print(nlyr(sentinel_multidate))

cat("\nGeometry check:\n")
print(
  compareGeom(
    sentinel_multidate[[1]],
    sentinel_multidate[[30]]
  )
)


# ------------------------------------------------------------
# 11. Load aspect-improved training polygons
# ------------------------------------------------------------

training <- st_read(
  gdb_file,
  layer = "Training_Larches_2017_aspect_improved2"
)


# ------------------------------------------------------------
# 12. Check reference polygons
# ------------------------------------------------------------

cat("\nReference polygons:\n")

print(
  table(
    training$Class,
    training$TrainValid_First
  )
)


# ------------------------------------------------------------
# 13. Separate training and validation polygons
# ------------------------------------------------------------

train_poly <- training[
  training$TrainValid_First == "t",
]

valid_poly <- training[
  training$TrainValid_First == "v",
]

cat("\nTraining polygons:\n")
print(table(train_poly$Class))

cat("\nValidation polygons:\n")
print(table(valid_poly$Class))


# ------------------------------------------------------------
# 14. Convert polygons to terra vectors
# ------------------------------------------------------------

train_vect <- vect(train_poly)

valid_vect <- vect(valid_poly)


# ------------------------------------------------------------
# 15. Extract training pixels
# ------------------------------------------------------------

train_values <- extract(
  sentinel_multidate,
  train_vect
)

train_df <- as.data.frame(
  train_values
)

train_df$Class <- train_poly$Class[
  train_df$ID
]

train_df <- na.omit(
  train_df
)

train_df$Class <- factor(
  train_df$Class,
  levels = c(0, 1)
)

cat("\nTraining pixels before balancing:\n")
print(table(train_df$Class))


# ------------------------------------------------------------
# 16. Extract validation pixels
# ------------------------------------------------------------

valid_values <- extract(
  sentinel_multidate,
  valid_vect
)

valid_df <- as.data.frame(
  valid_values
)

valid_df$Class <- valid_poly$Class[
  valid_df$ID
]

valid_df <- na.omit(
  valid_df
)

valid_df$Class <- factor(
  valid_df$Class,
  levels = c(0, 1)
)

cat("\nValidation pixels:\n")
print(table(valid_df$Class))


# ------------------------------------------------------------
# 17. Balance training pixels
# ------------------------------------------------------------

set.seed(123)

class_0 <- train_df[
  train_df$Class == "0",
]

class_1 <- train_df[
  train_df$Class == "1",
]

n_min <- min(
  nrow(class_0),
  nrow(class_1)
)

class_0_balanced <- class_0[
  sample(
    nrow(class_0),
    n_min
  ),
]

class_1_balanced <- class_1[
  sample(
    nrow(class_1),
    n_min
  ),
]

train_balanced <- rbind(
  class_0_balanced,
  class_1_balanced
)

cat("\nBalanced training pixels:\n")
print(table(train_balanced$Class))


# ------------------------------------------------------------
# 18. Train multi-date Random Forest
# ------------------------------------------------------------

predictor_names <- names(
  sentinel_multidate
)

rf_formula <- as.formula(
  paste(
    "Class ~",
    paste(
      predictor_names,
      collapse = " + "
    )
  )
)

rf_model_multidate <- ranger(
  formula = rf_formula,
  data = train_balanced,
  num.trees = 500,
  importance = "permutation",
  probability = TRUE,
  seed = 123
)

print(
  rf_model_multidate
)


# ------------------------------------------------------------
# 19. Save Random Forest model
# ------------------------------------------------------------

saveRDS(
  rf_model_multidate,
  file.path(
    output_dir,
    "RF_model_2017_multidate_aspect_improved.rds"
  )
)


# ------------------------------------------------------------
# 20. Variable importance
# ------------------------------------------------------------

variable_importance <- sort(
  rf_model_multidate$variable.importance,
  decreasing = TRUE
)

cat("\n--- VARIABLE IMPORTANCE ---\n")

print(
  variable_importance
)


# ------------------------------------------------------------
# 21. Predict validation pixels
# ------------------------------------------------------------

validation_prediction <- predict(
  rf_model_multidate,
  data = valid_df
)$predictions


# ------------------------------------------------------------
# 22. Test thresholds 0.50 - 0.99
# ------------------------------------------------------------

thresholds <- seq(
  0.50,
  0.99,
  by = 0.01
)

threshold_results <- data.frame(
  Threshold = thresholds,
  Accuracy = NA,
  Recall = NA,
  Precision = NA,
  F1 = NA
)

for (i in seq_along(thresholds)) {
  
  threshold <- thresholds[i]
  
  predicted <- ifelse(
    validation_prediction[, "1"] >= threshold,
    "1",
    "0"
  )
  
  predicted <- factor(
    predicted,
    levels = c("0", "1")
  )
  
  observed <- factor(
    valid_df$Class,
    levels = c("0", "1")
  )
  
  cm <- table(
    Observed = observed,
    Predicted = predicted
  )
  
  accuracy <- sum(
    diag(cm)
  ) / sum(cm)
  
  recall <- cm["1", "1"] /
    sum(cm["1", ])
  
  precision <- cm["1", "1"] /
    sum(cm[, "1"])
  
  if (
    is.na(precision) ||
    is.na(recall) ||
    (precision + recall) == 0
  ) {
    
    f1 <- NA
    
  } else {
    
    f1 <- 2 *
      precision *
      recall /
      (
        precision +
          recall
      )
  }
  
  threshold_results$Accuracy[i] <- accuracy * 100
  threshold_results$Recall[i] <- recall * 100
  threshold_results$Precision[i] <- precision * 100
  threshold_results$F1[i] <- f1 * 100
}


# ------------------------------------------------------------
# 23. Save threshold results
# ------------------------------------------------------------

write.csv(
  threshold_results,
  file.path(
    output_dir,
    "Threshold_results_2017_multidate_aspect_improved.csv"
  ),
  row.names = FALSE
)

print(
  threshold_results
)


# ------------------------------------------------------------
# 24. Select best threshold using highest F1
# ------------------------------------------------------------

best_row <- which.max(
  threshold_results$F1
)

threshold_final <- threshold_results$Threshold[
  best_row
]

cat("\n--- BEST THRESHOLD ---\n")

cat(
  "Threshold:",
  threshold_final,
  "\n"
)

cat(
  "Accuracy:",
  round(
    threshold_results$Accuracy[best_row],
    2
  ),
  "%\n"
)

cat(
  "Recall:",
  round(
    threshold_results$Recall[best_row],
    2
  ),
  "%\n"
)

cat(
  "Precision:",
  round(
    threshold_results$Precision[best_row],
    2
  ),
  "%\n"
)

cat(
  "F1:",
  round(
    threshold_results$F1[best_row],
    2
  ),
  "%\n"
)


# ------------------------------------------------------------
# 25. Final validation classification
# ------------------------------------------------------------

valid_df$Predicted <- ifelse(
  validation_prediction[, "1"] >= threshold_final,
  "1",
  "0"
)

valid_df$Predicted <- factor(
  valid_df$Predicted,
  levels = c("0", "1")
)

confusion <- table(
  Observed = valid_df$Class,
  Predicted = valid_df$Predicted
)

cat("\n--- CONFUSION MATRIX ---\n")

print(
  confusion
)


# ------------------------------------------------------------
# 26. Calculate validation metrics
# ------------------------------------------------------------

overall_accuracy <- sum(
  diag(confusion)
) / sum(confusion)

larch_recall <- confusion["1", "1"] /
  sum(confusion["1", ])

larch_precision <- confusion["1", "1"] /
  sum(confusion[, "1"])

larch_f1 <- 2 *
  larch_precision *
  larch_recall /
  (
    larch_precision +
      larch_recall
  )


# ------------------------------------------------------------
# 27. Print validation summary
# ------------------------------------------------------------

cat("\n--- MULTI-DATE RANDOM FOREST VALIDATION ---\n")

cat(
  "Threshold:",
  threshold_final,
  "\n"
)

cat(
  "Overall accuracy:",
  round(
    overall_accuracy * 100,
    2
  ),
  "%\n"
)

cat(
  "Larch recall:",
  round(
    larch_recall * 100,
    2
  ),
  "%\n"
)

cat(
  "Larch precision:",
  round(
    larch_precision * 100,
    2
  ),
  "%\n"
)

cat(
  "Larch F1:",
  round(
    larch_f1 * 100,
    2
  ),
  "%\n"
)


# ------------------------------------------------------------
# 28. Function for full-raster probability prediction
# ------------------------------------------------------------

predict_larch <- function(
    model,
    data
) {
  
  predictions <- predict(
    model,
    data = data
  )$predictions
  
  return(
    predictions[, "1"]
  )
}


# ------------------------------------------------------------
# 29. Predict larch probability for entire SNP forest
# ------------------------------------------------------------

larch_probability_multidate <- terra::predict(
  sentinel_multidate,
  rf_model_multidate,
  fun = predict_larch,
  na.rm = TRUE
)

names(
  larch_probability_multidate
) <- "Larch_Probability"


# ------------------------------------------------------------
# 30. Save probability raster
# ------------------------------------------------------------

writeRaster(
  larch_probability_multidate,
  file.path(
    output_dir,
    "Larch_Probability_2017_multidate_aspect_improved.tif"
  ),
  overwrite = TRUE
)


# ------------------------------------------------------------
# 31. Create binary classification
# ------------------------------------------------------------

larch_binary_multidate <- larch_probability_multidate >= threshold_final

names(
  larch_binary_multidate
) <- "Larch"


# ------------------------------------------------------------
# 32. Save binary classification
# ------------------------------------------------------------

writeRaster(
  larch_binary_multidate,
  file.path(
    output_dir,
    "Larch_Classification_2017_multidate_aspect_improved.tif"
  ),
  overwrite = TRUE
)


# ------------------------------------------------------------
# 33. Calculate preliminary larch area
# ------------------------------------------------------------

n_larch_pixels <- global(
  larch_binary_multidate,
  "sum",
  na.rm = TRUE
)[1, 1]

# 10 m x 10 m = 100 m2

larch_area_m2 <- n_larch_pixels * 100

larch_area_ha <- larch_area_m2 / 10000


# ------------------------------------------------------------
# 34. Final summary
# ------------------------------------------------------------

cat("\n========================================\n")

cat("MULTI-DATE ASPECT-IMPROVED CLASSIFICATION\n")

cat("========================================\n")

cat(
  "\nReference polygons:\n"
)

print(
  table(
    training$Class,
    training$TrainValid_First
  )
)

cat(
  "\nSelected threshold:",
  threshold_final,
  "\n"
)

cat(
  "Overall accuracy:",
  round(
    overall_accuracy * 100,
    2
  ),
  "%\n"
)

cat(
  "Larch recall:",
  round(
    larch_recall * 100,
    2
  ),
  "%\n"
)

cat(
  "Larch precision:",
  round(
    larch_precision * 100,
    2
  ),
  "%\n"
)

cat(
  "Larch F1:",
  round(
    larch_f1 * 100,
    2
  ),
  "%\n"
)

cat(
  "Pixels classified as larch:",
  n_larch_pixels,
  "\n"
)

cat(
  "Preliminary classified larch area:",
  round(
    larch_area_ha,
    2
  ),
  "ha\n"
)

cat(
  "\nTop variable importance:\n"
)

print(
  head(
    variable_importance,
    15
  )
)