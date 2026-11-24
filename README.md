# SNP_Larch_Classification_ML
Random Forest classification of larches in the Swiss National Park using Sentinel-2 imagery

# Larch classification from Sentinel-2

Mapping European larch (Larix decidua) in the Swiss National 
Park (SNP) using Sentinel-2 imagery and Random Forest machine 
learning. The resulting map serves as the area of interest 
(AOI) for the subsequent autumn phenology analysis.

## Workflow

- Preparation of Sentinel-2 imagery from 4, 14 and 19 October 2017
- Manual digitisation of training and validation polygons
- Random Forest classification (500 trees, 30 spectral predictors)
- Threshold evaluation and selection (0.99)
- Generation of a binary larch map at 10 m resolution
- Validation and spatial comparison with the SNP forest-type map

## Setup

Required R packages:

```r
install.packages(c("terra", "sf", "ranger"))
```

Input datasets must be obtained separately. File paths 
in the scripts need to be adapted before running the code.

## Data sources

- Sentinel-2 L2A (ESA / Copernicus), via Microsoft Planetary Computer
- Copernicus DEM GLO-30
- swisstopo SwissTLM3D forest layer
- Manually digitised training and validation polygons
- SNP forest-type map (HABITALP project, 2002–2006)
- BAFU: Swiss National Park boundary

## Output

Binary larch classification GeoTIFF at 10 m resolution, 
covering approximately 529 ha of classified larch.

Validation accuracy: 98.58%, precision: 88.84%, 
recall: 96.98%, F1 score: 92.73%.

Validation is based on manually labelled polygons rather 
than independent ground observations.

## Related repository

[Larch autumn phenology from Sentinel-2](https://github.com/hannahbeu/phenology)

