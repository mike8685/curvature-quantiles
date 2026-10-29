# Raw data for the real-data benchmark

Raw files are not committed. Put them under `data/raw/` (or point `REALBENCH_DATA`
elsewhere) with this layout. `python -m realbench datasets` shows what it finds.
At freeze, every file used is hashed (SHA-256) into `realbench/registry.json`, so
anyone can check they hold exactly the same bytes.

The six confirmatory datasets are fixed in `realbench/config.py`. Download each corpus
**whole**, exactly as published. Never hand-pick files or curves.

| path | dataset id | where to get it | notes |
| --- | --- | --- | --- |
| `natural_earth/ne_10m_coastline.geojson` | `ne_coastline` (confirmatory) | github.com/nvkelso/natural-earth-vector, `geojson/` | public domain |
| `natural_earth/ne_10m_rivers_lake_centerlines.geojson` | `ne_rivers` (confirmatory) | same | |
| `natural_earth/ne_10m_roads.geojson` | `ne_roads` (confirmatory) | same | |
| `geolife/Data/<user>/Trajectory/*.plt` | `geolife` (confirmatory) | Microsoft Research "GeoLife GPS Trajectories 1.3" (unzip as-is) | 182 users = clusters |
| `quickdraw/<category>.ndjson` | `quickdraw_raw` (confirmatory) | `gs://quickdraw_dataset/full/raw/<category>.ndjson` for exactly: moustache, sun, binoculars, stop sign, paint can, strawberry, crocodile, bandage, ambulance, nose | **raw** files only; `simplified` files are refused; other categories are ignored |
| `fusion360/**/*.json` | `fusion360` (confirmatory) | Fusion 360 Gallery "reconstruction" dataset (github.com/AutodeskAILab/Fusion360GalleryDataset, download link in its README) | all design JSON files; loops with splines are skipped and counted |
| `gpx/**/*.gpx` | `gpx` (exploratory) | gpxpy `test_files/` | parser fixtures, not a sample |
| `omniglot/strokes_{background,evaluation}/<alphabet>/<character>/*.txt` | `omniglot` (exploratory) | github.com/brendenlake/omniglot, `python/strokes_background.zip` and `strokes_evaluation.zip` (MIT) | 1623 characters x 20 drawers; cluster = character; strokes are pixel-quantised |
| `inkml/**/*.inkml` | `inkml` (exploratory) | CROHME (Zenodo record 8428035) or Google MathWriting | |
| `uji/*.txt` | `uji` (exploratory) | UCI "UJI Pen Characters (Version 2)", `ujipenchars2.txt` | |
| `dxf/**/*.dxf` | `dxf` (exploratory) | any DXF corpus | |

Quick Draw file names: the bucket uses the category name with spaces (for example
`stop sign.ndjson`). Spaces or underscores are both accepted.
