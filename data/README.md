# Raw data for the real-data benchmark

Raw files are not committed. Put them under `data/raw/` (or point `REALBENCH_DATA`
elsewhere) with this layout. `python -m realbench datasets` shows what it finds.
At freeze, every file used is hashed (SHA-256) into `realbench/registry.json`, so
anyone can check they hold exactly the same bytes.

| path | dataset id | where to get it | notes |
| --- | --- | --- | --- |
| `natural_earth/ne_10m_coastline.geojson` | `ne_coastline` | github.com/nvkelso/natural-earth-vector, `geojson/` | public domain |
| `natural_earth/ne_10m_rivers_lake_centerlines.geojson` | `ne_rivers` | same | |
| `natural_earth/ne_10m_roads.geojson` | `ne_roads` | same | |
| `gpx/**/*.gpx` | `gpx` | any real recorded GPX: your own watch exports, OpenStreetMap public GPS traces, gpxpy `test_files/` | one cluster per file; need 30+ files in the test split for a confirmatory result |
| `quickdraw/*.ndjson` | `quickdraw_raw` | Quick Draw **raw** files: `gs://quickdraw_dataset/full/raw/<category>.ndjson` | never the `simplified` files (RDP-processed; the loader refuses them) |
| `inkml/**/*.inkml` | `inkml` | CROHME (Zenodo record 8428035) or Google MathWriting | one cluster per file |
| `uji/*.txt` | `uji` | UCI "UJI Pen Characters (Version 2)", `ujipenchars2.txt` | parser written from the format description; check its counts on first run |
| `fusion360/**/*.json` | `fusion360` | Fusion 360 Gallery "reconstruction" dataset (github.com/AutodeskAILab/Fusion360GalleryDataset) | arcs/circles parsed by Fusion API field names; loops with splines are skipped and counted |
| `dxf/**/*.dxf` | `dxf` | any DXF corpus | fallback CAD source |

Choose which categories to download before looking at any results, and download
whole files or whole categories. Never hand-pick curves.
