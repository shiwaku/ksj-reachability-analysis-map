#!/usr/bin/env python3
"""
parquet → GeoJSON 変換スクリプト（tippecanoe への入力用）

output/arrival_map_*.parquet を読み込み、到達可能メッシュのみ
output/geojson/arrival_map_{name}.geojson として出力する。
origins_*.geojson が存在すれば output/geojson/origins.geojson にコピーする。
"""

import shutil
from pathlib import Path

import geopandas as gpd

REPO_ROOT = Path(__file__).parent.parent
output_dir = REPO_ROOT / "output"
geojson_dir = output_dir / "geojson"
geojson_dir.mkdir(parents=True, exist_ok=True)

for p in sorted(output_dir.glob("arrival_map_*.parquet")):
    name = p.stem.replace("arrival_map_", "")
    gdf = gpd.read_parquet(p)
    gdf = gdf[gdf["dist_min"].notna()].copy()
    gdf["origin"] = name
    gdf = gdf[["mesh_code", "dist_min", "dist_rank", "origin", "geometry"]]
    out = geojson_dir / f"arrival_map_{name}.geojson"
    gdf.to_file(out, driver="GeoJSON")
    print(f"  {out.name}: {len(gdf):,} features")

origins_files = sorted(output_dir.glob("origins_*.geojson"))
if origins_files:
    dest = geojson_dir / "origins.geojson"
    shutil.copy(origins_files[0], dest)
    print(f"  origins.geojson: {origins_files[0].name} からコピー")
else:
    print("  警告: origins_*.geojson が見つかりません")
