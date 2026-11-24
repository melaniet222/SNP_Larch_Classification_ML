# ============================================================
# LARCH CLASSIFICATION - SENTINEL-2 - SNP - 2017
# ASPECT-IMPROVED MODEL 2
# ============================================================


# ------------------------------------------------------------
# 1. Load packages
# ------------------------------------------------------------

library(terra)
library(sf)
library(ranger)


# ------------------------------------------------------------
# 2. Define file paths
# ------------------------------------------------------------

raster_dir <- "S:/course/geo715/student/phenology/gis_work/RastersGIS"

gdb_file <- "S:/course/geo715/student/phenology/gis_work/Phenology.gdb"

output_dir <- "S:/course/geo715/student/phenology/gis_work/RastersGIS"


# ------------------------------------------------------------
# 3. Load Sentinel-2 bands
# ------------------------------------------------------------

b02 <- rast(file.path(
  raster_dir,
  "S2_20171014_10m_B02_SNP_Forest.tif"
))

b03 <- rast(file.path(
  raster_dir,
  "S2_20171014_10m_B03_SNP_Forest.tif"
))

b04 <- rast(file.path(
  raster_dir,
  "S2_20171014_10m_B04_SNP_Forest.tif"
))

b05 <- rast(file.path(
  raster_dir,
  "S2_20171014_10m_B05_SNP_Forest.tif"
))

b06 <- rast(file.path(
  raster_dir,
  "S2_20171014_10m_B06_SNP_Forest.tif"
))

b07 <- rast(file.path(
  raster_dir,
  "S2_20171014_10m_B07_SNP_Forest.tif"
))

b08 <- rast(file.path(
  raster_dir,
  "S2_20171014_10m_B08_SNP_Forest.tif"
))

b8a <- rast(file.path(
  raster_dir,
  "S2_20171014_10m_B8A_SNP_Forest.tif"
))

b11 <- rast(file.path(
  raster_dir,
  "S2_20171014_10m_B11_SNP_Forest.tif"
))

b12 <- rast(file.path(
  raster_dir,
  "S2_20171014_10m_B12_SNP_Forest.tif"
))


# ------------------------------------------------------------
# 4. Align all bands to B08
# ------------------------------------------------------------

template <- b08

b02 <- resample(b02, template, method = "bilinear")
b03 <- resample(b03, template, method = "bilinear")
b04 <- resample(b04, template, method = "bilinear")
b05 <- resample(b05, template, method = "bilinear")
b06 <- resample(b06, template, method = "bilinear")
b07 <- resample(b07, template, method = "bilinear")

# B08 is already the reference raster

b8a <- resample(b8a, template, method = "bilinear")
b11 <- resample(b11, template, method = "bilinear")
b12 <- resample(b12, template, method = "bilinear")


# ------------------------------------------------------------
# 5. Combine bands into one raster stack
# ------------------------------------------------------------

sentinel <- c(
  b02,
  b03,
  b04,
  b05,
  b06,
  b07,
  b08,
  b8a,
  b11,
  b12
)

names(sentinel) <- c(
  "B02",
  "B03",
  "B04",
  "B05",
  "B06",
  "B07",
  "B08",
  "B8A",
  "B11",
  "B12"
)


# ------------------------------------------------------------
# 6. Check raster geometry
# ------------------------------------------------------------

print(sentinel)

compareGeom(
  b02, b03, b04, b05, b06,
  b07, b08, b8a, b11, b12
)


# ------------------------------------------------------------
# 7. Load UPDATED training polygons
# ------------------------------------------------------------

training <- st_read(
  gdb_file,
  layer = "Training_Larches_2017_aspect_improved2"
)

print(training)


# ------------------------------------------------------------
# 8. Check coordinate system and field names
# ------------------------------------------------------------

st_crs(training)

names(training)


# ------------------------------------------------------------
# 9. Check reference polygon counts
# ------------------------------------------------------------

cat("\n--- UPDATED REFERENCE POLYGONS ---\n")

cat("\nTotal polygons per class:\n")
print(
  table(training$Class)
)

cat("\nTraining / Validation split:\n")
print(
  table(training$TrainValid_First)
)

cat("\nClass x Training/Validation:\n")
print(
  table(
    training$Class,
    training$TrainValid_First
  )
)


# ------------------------------------------------------------
# 10. Separate Training and Validation polygons
# ------------------------------------------------------------

train_poly <- training[
  training$TrainValid_First == "t",
]

valid_poly <- training[
  training$TrainValid_First == "v",
]

cat("\nTraining polygons per class:\n")
print(
  table(train_poly$Class)
)

cat("\nValidation polygons per class:\n")
print(
  table(valid_poly$Class)
)


# ------------------------------------------------------------
# 11. Convert polygons to terra vectors
# ------------------------------------------------------------

train_vect <- vect(train_poly)

valid_vect <- vect(valid_poly)


# ------------------------------------------------------------
# 12. Extract Sentinel values from training polygons
# ------------------------------------------------------------

train_values <- extract(
  sentinel,
  train_vect
)

train_df <- as.data.frame(
  train_values
)


# ------------------------------------------------------------
# 13. Add class information using polygon ID
# ------------------------------------------------------------

train_df$Class <- train_poly$Class[
  train_df$ID
]


# ------------------------------------------------------------
# 14. Remove NA values
# ------------------------------------------------------------

train_df <- na.omit(
  train_df
)


# ------------------------------------------------------------
# 15. Convert Class to factor
# ------------------------------------------------------------

train_df$Class <- factor(
  train_df$Class,
  levels = c(0, 1)
)


# ------------------------------------------------------------
# 16. Check training pixel counts
# ------------------------------------------------------------

cat("\nTraining pixels per class before balancing:\n")

print(
  table(train_df$Class)
)


# ------------------------------------------------------------
# 17. Extract Sentinel values from validation polygons
# ------------------------------------------------------------

valid_values <- extract(
  sentinel,
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

cat("\nValidation pixels per class:\n")

print(
  table(valid_df$Class)
)


# ------------------------------------------------------------
# 18. Balance training pixels
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

print(
  table(train_balanced$Class)
)


# ------------------------------------------------------------
# 19. Train Random Forest
# ------------------------------------------------------------

rf_model <- ranger(
  Class ~ B02 + B03 + B04 + B05 + B06 +
    B07 + B08 + B8A + B11 + B12,
  data = train_balanced,
  num.trees = 500,
  importance = "permutation",
  probability = TRUE,
  seed = 123
)

print(
  rf_model
)


# ------------------------------------------------------------
# 20. Save Random Forest model
# ------------------------------------------------------------

saveRDS(
  rf_model,
  file.path(
    output_dir,
    "RF_model_2017_aspect_improved2.rds"
  )
)


# ------------------------------------------------------------
# 21. Variable importance
# ------------------------------------------------------------

variable_importance <- sort(
  rf_model$variable.importance,
  decreasing = TRUE
)

cat("\n--- VARIABLE IMPORTANCE ---\n")

print(
  variable_importance
)


# ------------------------------------------------------------
# 22. Predict validation data
# ------------------------------------------------------------

validation_prediction <- predict(
  rf_model,
  data = valid_df
)$predictions

head(
  validation_prediction
)


# ------------------------------------------------------------
# 23. Test thresholds from 0.50 to 0.99
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
    sum(
      cm["1", ]
    )
  
  precision <- cm["1", "1"] /
    sum(
      cm[, "1"]
    )
  
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
      (precision + recall)
  }
  
  threshold_results$Accuracy[i] <- accuracy * 100
  
  threshold_results$Recall[i] <- recall * 100
  
  threshold_results$Precision[i] <- precision * 100
  
  threshold_results$F1[i] <- f1 * 100
}


# ------------------------------------------------------------
# 24. Show threshold results
# ------------------------------------------------------------

print(
  threshold_results
)


# ------------------------------------------------------------
# 25. Save threshold results
# ------------------------------------------------------------

write.csv(
  threshold_results,
  file.path(
    output_dir,
    "Threshold_results_2017_aspect_improved2.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 26. Select threshold with highest F1 score
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
  "Larch recall:",
  round(
    threshold_results$Recall[best_row],
    2
  ),
  "%\n"
)

cat(
  "Larch precision:",
  round(
    threshold_results$Precision[best_row],
    2
  ),
  "%\n"
)

cat(
  "F1 score:",
  round(
    threshold_results$F1[best_row],
    2
  ),
  "%\n"
)


# ------------------------------------------------------------
# 27. Final validation classification
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


# ------------------------------------------------------------
# 28. Confusion matrix
# ------------------------------------------------------------

confusion <- table(
  Observed = valid_df$Class,
  Predicted = valid_df$Predicted
)

cat("\n--- FINAL CONFUSION MATRIX ---\n")

print(
  confusion
)


# ------------------------------------------------------------
# 29. Final accuracy metrics
# ------------------------------------------------------------

overall_accuracy <- sum(
  diag(confusion)
) / sum(confusion)

larch_recall <- confusion["1", "1"] /
  sum(
    confusion["1", ]
  )

larch_precision <- confusion["1", "1"] /
  sum(
    confusion[, "1"]
  )

larch_f1 <- 2 *
  larch_precision *
  larch_recall /
  (
    larch_precision +
      larch_recall
  )


# ------------------------------------------------------------
# 30. Print validation results
# ------------------------------------------------------------

cat("\n--- RANDOM FOREST VALIDATION ---\n")

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
# 31. Function to predict larch probability
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
# 32. Predict larch probability for entire SNP forest
# ------------------------------------------------------------

larch_probability <- terra::predict(
  sentinel,
  rf_model,
  fun = predict_larch,
  na.rm = TRUE
)

names(
  larch_probability
) <- "Larch_Probability"


# ------------------------------------------------------------
# 33. Plot probability raster
# ------------------------------------------------------------

plot(
  larch_probability,
  main = "Larch probability - 2017 - aspect improved 2"
)


# ------------------------------------------------------------
# 34. Save probability raster
# ------------------------------------------------------------

writeRaster(
  larch_probability,
  file.path(
    output_dir,
    "Larch_Probability_2017_aspect_improved2.tif"
  ),
  overwrite = TRUE
)


# ------------------------------------------------------------
# 35. Create binary classification
# ------------------------------------------------------------

larch_binary <- larch_probability >= threshold_final

names(
  larch_binary
) <- "Larch"


# ------------------------------------------------------------
# 36. Plot binary classification
# ------------------------------------------------------------

plot(
  larch_binary,
  main = paste(
    "Larch classification - threshold",
    threshold_final
  )
)


# ------------------------------------------------------------
# 37. Save binary classification
# ------------------------------------------------------------

writeRaster(
  larch_binary,
  file.path(
    output_dir,
    "Larch_Classification_2017_aspect_improved2.tif"
  ),
  overwrite = TRUE
)


# ------------------------------------------------------------
# 38. Calculate preliminary classified larch area
# ------------------------------------------------------------

n_larch_pixels <- global(
  larch_binary,
  "sum",
  na.rm = TRUE
)[1, 1]

# 10 m x 10 m = 100 m2

larch_area_m2 <- n_larch_pixels * 100

larch_area_ha <- larch_area_m2 / 10000


# ------------------------------------------------------------
# 39. Final summary
# ------------------------------------------------------------

cat("\n========================================\n")

cat("ASPECT-IMPROVED LARCH CLASSIFICATION 2\n")

cat("========================================\n")


cat("\nReference polygons:\n")

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

cat("\nVariable importance:\n")

print(
  variable_importance
)