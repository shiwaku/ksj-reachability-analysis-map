#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
到達圏分析スクリプト（複数始点対応）

複数の始点から各 L6 メッシュへの最短到達時間を Multi-source Dijkstra で算出し、
GeoParquet 形式で出力する。

【アルゴリズム】
  スーパーソースノードを追加し全始点を同時キューに投入。
  scipy.sparse.csgraph.dijkstra を 1 回実行して最近傍始点からの最短到達時間を算出。

【使い方】
  python3 src/reachability_search.py
  python3 src/reachability_search.py --orig 35.8578,139.6490,埼玉県庁
  python3 src/reachability_search.py \\
      --orig 35.8578,139.6490,埼玉県庁 \\
      --orig 36.0420,139.4006,東松山市役所
  python3 src/reachability_search.py --orig-csv origins.csv
"""

import argparse
import csv
import json
import time
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import dijkstra as sp_dijkstra
from scipy.spatial import KDTree
from shapely.geometry import box

SCRIPT_DIR     = Path(__file__).parent
REPO_ROOT      = SCRIPT_DIR.parent
SAMPLE_DIR     = REPO_ROOT / "network" / "saitama_pref"
DEFAULT_LINKS  = str(SAMPLE_DIR / "KSJ_N13-24_saitama_pref_道路リンク.parquet")
DEFAULT_NODES  = str(SAMPLE_DIR / "KSJ_N13-24_saitama_pref_道路ノード.parquet")
DEFAULT_ACCESS = str(SAMPLE_DIR / "KSJ_N13-24_saitama_pref_アクセスリンク_L6.parquet")
DEFAULT_OUT    = str(REPO_ROOT / "output")

VEHICLE_BINS   = list(range(0, 100, 10)) + [float("inf")]
VEHICLE_LABELS = list(range(len(VEHICLE_BINS) - 1))
VEHICLE_COLORS = [
    (255,   0,   0),
    (255,  64,   0),
    (255, 128,   0),
    (255, 192,   0),
    (255, 255,   0),
    (192, 255,   0),
    (  0, 204,   0),
    (  0, 204, 128),
    (  0, 204, 204),
    ( 68,   0,  85),
]


def compute_l6_polygons(codes_series: pd.Series) -> gpd.GeoDataFrame:
    codes = codes_series.astype(str).str.zfill(11)
    pp = codes.str[0:2].astype(int).values
    qq = codes.str[2:4].astype(int).values
    r  = codes.str[4].astype(int).values
    s  = codes.str[5].astype(int).values
    t  = codes.str[6].astype(int).values
    u  = codes.str[7].astype(int).values
    v  = codes.str[8].astype(int).values
    w  = codes.str[9].astype(int).values
    x  = codes.str[10].astype(int).values

    lat = pp / 1.5
    lon = (qq + 100).astype(float)
    dlat2 = (2.0/3.0)/8;  dlon2 = 1.0/8
    lat += r * dlat2;     lon += s * dlon2
    dlat3 = dlat2/10;     dlon3 = dlon2/10
    lat += t * dlat3;     lon += u * dlon3
    dlat4 = dlat3/2;      dlon4 = dlon3/2
    lat += ((v-1)//2) * dlat4;  lon += ((v-1)%2) * dlon4
    dlat5 = dlat4/2;      dlon5 = dlon4/2
    lat += ((w-1)//2) * dlat5;  lon += ((w-1)%2) * dlon5
    dlat6 = dlat5/2;      dlon6 = dlon5/2
    lat += ((x-1)//2) * dlat6;  lon += ((x-1)%2) * dlon6

    geoms = [box(lo, la, lo+dlon6, la+dlat6) for lo, la in zip(lon, lat)]
    return gpd.GeoDataFrame({"mesh_code": codes.values}, geometry=geoms, crs="EPSG:4326")


def build_sparse_graph(links: gpd.GeoDataFrame):
    n1 = links["node1"].astype(int).to_numpy()
    n2 = links["node2"].astype(int).to_numpy()
    ws = (links["time_001min"].astype(float) * 0.01).to_numpy(dtype=np.float64)

    src = np.concatenate([n1, n2])
    dst = np.concatenate([n2, n1])
    w   = np.concatenate([ws, ws])

    unique = np.unique(np.concatenate([src, dst]))
    n2i    = {int(n): i for i, n in enumerate(unique.tolist())}
    n_v    = len(unique)

    rows = np.array([n2i[int(n)] for n in src], dtype=np.int32)
    cols = np.array([n2i[int(n)] for n in dst], dtype=np.int32)
    G    = csr_matrix((w, (rows, cols)), shape=(n_v, n_v))
    return unique, n2i, G


def nearest_road_node(nodes: gpd.GeoDataFrame, lat: float, lon: float):
    coords = np.array([[p.y, p.x] for p in nodes.geometry])
    dist, i = KDTree(coords).query([lat, lon])
    row = nodes.iloc[i]
    return int(row["node_id"]), float(row.geometry.y), float(row.geometry.x), float(dist * 111000)


def write_arrival_qml(qml_path: Path) -> None:
    n = len(VEHICLE_LABELS)

    def label_text(i):
        lo = int(VEHICLE_BINS[i])
        hi = VEHICLE_BINS[i + 1]
        return f"{lo}分超" if hi == float("inf") else f"{lo}〜{int(hi)}分"

    cats = "\n".join(
        f'      <category symbol="{i}" value="{VEHICLE_LABELS[i]}" '
        f'label="{label_text(i)}" render="true"/>'
        for i in range(n)
    )
    cats += f'\n      <category symbol="{n}" value="" label="到達不能" render="true"/>'

    syms = []
    for i, (r, g, b) in enumerate(VEHICLE_COLORS):
        syms.append(
            f'      <symbol name="{i}" type="fill" alpha="0.75" '
            f'clip_to_extent="1" is_animated="0" frame_rate="10">\n'
            f'        <data_defined_properties><Option type="Map">'
            f'<Option name="name" type="QString" value=""/>'
            f'<Option name="properties"/>'
            f'<Option name="type" type="QString" value="collection"/>'
            f'</Option></data_defined_properties>\n'
            f'        <layer class="SimpleFill" enabled="1" pass="0" locked="0">\n'
            f'          <Option type="Map">\n'
            f'            <Option name="color" type="QString" value="{r},{g},{b},255"/>\n'
            f'            <Option name="outline_style" type="QString" value="no"/>\n'
            f'            <Option name="style" type="QString" value="solid"/>\n'
            f'          </Option>\n'
            f'        </layer>\n'
            f'      </symbol>'
        )
    syms.append(
        f'      <symbol name="{n}" type="fill" alpha="0.75" '
        f'clip_to_extent="1" is_animated="0" frame_rate="10">\n'
        f'        <data_defined_properties><Option type="Map">'
        f'<Option name="name" type="QString" value=""/>'
        f'<Option name="properties"/>'
        f'<Option name="type" type="QString" value="collection"/>'
        f'</Option></data_defined_properties>\n'
        f'        <layer class="SimpleFill" enabled="1" pass="0" locked="0">\n'
        f'          <Option type="Map">\n'
        f'            <Option name="color" type="QString" value="170,170,170,180"/>\n'
        f'            <Option name="outline_style" type="QString" value="no"/>\n'
        f'            <Option name="style" type="QString" value="solid"/>\n'
        f'          </Option>\n'
        f'        </layer>\n'
        f'      </symbol>'
    )

    content = (
        "<!DOCTYPE qgis PUBLIC 'http://mrcc.com/qgis.dtd' 'SYSTEM'>\n"
        '<qgis version="3.34.0" styleCategories="Symbology">\n'
        '  <renderer-v2 type="categorizedSymbol" attr="dist_rank" '
        'forceraster="0" symbollevels="0" usingSymbolLevels="0" enableorderby="0">\n'
        '    <categories>\n'
        f'{cats}\n'
        '    </categories>\n'
        '    <symbols>\n'
        + "\n".join(syms) + "\n"
        '    </symbols>\n'
        '    <rotation/>\n'
        '    <sizescale/>\n'
        '  </renderer-v2>\n'
        '  <blendMode>0</blendMode>\n'
        '  <featureBlendMode>0</featureBlendMode>\n'
        '  <layerOpacity>1</layerOpacity>\n'
        '</qgis>\n'
    )
    qml_path.write_text(content, encoding="utf-8")


def save_arrival_map(mesh_codes, dist_mesh, out_parquet: Path, out_qml: Path) -> None:
    reachable = np.isfinite(dist_mesh)
    dist_safe = np.where(reachable, dist_mesh, 0.0)
    rank_idx  = np.minimum(np.digitize(dist_safe, VEHICLE_BINS[1:]), len(VEHICLE_LABELS) - 1)
    rank_str  = np.where(reachable, rank_idx.astype(str), "")

    df = pd.DataFrame({
        "mesh_code": mesh_codes,
        "dist_min":  np.where(reachable, np.round(dist_mesh, 2), np.nan),
        "dist_rank": rank_str,
    })

    mesh_gdf = compute_l6_polygons(pd.Series(mesh_codes))
    gdf = mesh_gdf.merge(df, on="mesh_code").set_crs("EPSG:4326")
    gdf.to_parquet(out_parquet)
    write_arrival_qml(out_qml)

    reachable_cnt = int(reachable.sum())
    sz = out_parquet.stat().st_size // 1024
    print(f"  {out_parquet.name}  ({reachable_cnt:,}メッシュ到達可能, {sz}KB)")
    print(f"  {out_qml.name}")


def write_od_qml(qml_path: Path) -> None:
    content = (
        '<!DOCTYPE qgis PUBLIC \'http://mrcc.com/qgis.dtd\' \'SYSTEM\'>\n'
        '<qgis version="3.34.0" styleCategories="Symbology|Labeling">\n'
        '  <renderer-v2 type="singleSymbol" forceraster="0" symbollevels="0"'
        ' usingSymbolLevels="0" enableorderby="0">\n'
        '    <symbols>\n'
        '      <symbol name="0" type="marker" alpha="1" clip_to_extent="1"'
        ' is_animated="0" frame_rate="10">\n'
        '        <data_defined_properties><Option type="Map">'
        '<Option name="name" type="QString" value=""/>'
        '<Option name="properties"/>'
        '<Option name="type" type="QString" value="collection"/>'
        '</Option></data_defined_properties>\n'
        '        <layer class="SimpleMarker" enabled="1" pass="0" locked="0">\n'
        '          <Option type="Map">\n'
        '            <Option name="color" type="QString" value="255,255,255,255"/>\n'
        '            <Option name="outline_color" type="QString" value="0,0,0,255"/>\n'
        '            <Option name="outline_width" type="QString" value="0.4"/>\n'
        '            <Option name="outline_width_unit" type="QString" value="MM"/>\n'
        '            <Option name="size" type="QString" value="3"/>\n'
        '            <Option name="size_unit" type="QString" value="MM"/>\n'
        '            <Option name="name" type="QString" value="circle"/>\n'
        '          </Option>\n'
        '        </layer>\n'
        '      </symbol>\n'
        '    </symbols>\n'
        '    <rotation/>\n'
        '    <sizescale/>\n'
        '  </renderer-v2>\n'
        '  <labeling type="simple">\n'
        '    <settings calloutType="simple">\n'
        '      <text-style fieldName="name" fontFamily="sans-serif" fontSize="10"'
        ' fontWeight="75" textColor="0,0,0,255" namedStyle="Bold"'
        ' textOpacity="1" blendMode="0" isExpression="0"'
        ' fontLetterSpacing="0" fontWordSpacing="0"'
        ' fontUnderline="0" fontStrikeout="0" fontItalic="0"'
        ' fontSizeUnit="Point" fontSizeMapUnitScale="3x:0,0,0,0,0,0">\n'
        '        <text-buffer bufferDraw="1" bufferSize="1" bufferSizeUnits="MM"'
        ' bufferColor="255,255,255,255" bufferOpacity="1" bufferBlendMode="0"'
        ' bufferNoFill="0" bufferJoinStyle="128"/>\n'
        '        <background shapeDraw="0"/>\n'
        '        <shadow shadowDraw="0"/>\n'
        '      </text-style>\n'
        '      <text-format autoWrapLength="0" useMaxLineLengthForAutoWrap="1"'
        ' addDirectionSymbol="0" leftDirectionSymbol="&lt;"'
        ' rightDirectionSymbol="&gt;" reverseDirectionSymbol="0"'
        ' placeDirectionSymbol="0" formatNumbers="0" decimals="3"'
        ' plusSign="0" multilineAlign="3"/>\n'
        '      <placement placement="2" offsetType="0" xOffset="0" yOffset="0"'
        ' offsetUnits="MM" dist="2" distUnits="MM" distMapUnitScale="3x:0,0,0,0,0,0"'
        ' repeatDistance="0" repeatDistanceUnits="MM"'
        ' repeatDistanceMapUnitScale="3x:0,0,0,0,0,0"'
        ' maxCurvedCharAngleIn="25" maxCurvedCharAngleOut="-25"'
        ' priority="5" predefinedPositionOrder="TR,TL,BR,BL,R,L,TSR,BSR"'
        ' fitInPolygonOnly="0" overrunDistance="0" overrunDistanceUnit="MM"'
        ' overrunDistanceMapUnitScale="3x:0,0,0,0,0,0"'
        ' labelOffsetMapUnitScale="3x:0,0,0,0,0,0"'
        ' polygonPlacementFlags="2" allowDegraded="0"'
        ' geometryGenerator="" geometryGeneratorEnabled="0"'
        ' geometryGeneratorType="PointGeometry" layerType="PointGeometry"'
        ' centroidWhole="0" centroidInside="0"'
        ' overlapHandling="PreventOverlap" zIndex="0"/>\n'
        '      <rendering obstacle="1" obstacleFactor="1" obstacleType="1"'
        ' scaleVisibility="0" minScale="1" maxScale="0"'
        ' limitNumLabels="0" maxNumLabels="2000"'
        ' displayAll="0" upsidedownLabels="0"'
        ' fontMinPixelSize="3" fontMaxPixelSize="10000"'
        ' mergeLines="0" drawLabels="1" labelPerPart="0"'
        ' scaleMin="0" scaleMax="0"/>\n'
        '      <dd_properties>\n'
        '        <Option type="Map"><Option name="name" type="QString" value=""/>'
        '<Option name="properties"/>'
        '<Option name="type" type="QString" value="collection"/></Option>\n'
        '      </dd_properties>\n'
        '      <callout type="simple">\n'
        '        <Option type="Map"><Option name="anchorPoint" type="QString" value="pole_of_inaccessibility"/>'
        '<Option name="blendMode" type="int" value="0"/>'
        '<Option name="ddProperties" type="Map">'
        '<Option name="name" type="QString" value=""/>'
        '<Option name="properties"/>'
        '<Option name="type" type="QString" value="collection"/>'
        '</Option>'
        '<Option name="drawToAllParts" type="bool" value="false"/>'
        '<Option name="enabled" type="QString" value="0"/>'
        '<Option name="labelAnchorPoint" type="QString" value="point_on_exterior"/>'
        '<Option name="lineSymbol" type="QString" value="&lt;symbol name=&quot;_&quot; type=&quot;line&quot; alpha=&quot;1&quot; clip_to_extent=&quot;1&quot; is_animated=&quot;0&quot; frame_rate=&quot;10&quot;&gt;&lt;data_defined_properties&gt;&lt;Option type=&quot;Map&quot;&gt;&lt;Option name=&quot;name&quot; type=&quot;QString&quot; value=&quot;&quot;/&gt;&lt;Option name=&quot;properties&quot;/&gt;&lt;Option name=&quot;type&quot; type=&quot;QString&quot; value=&quot;collection&quot;/&gt;&lt;/Option&gt;&lt;/data_defined_properties&gt;&lt;layer class=&quot;SimpleLine&quot; enabled=&quot;1&quot; pass=&quot;0&quot; locked=&quot;0&quot;&gt;&lt;Option type=&quot;Map&quot;&gt;&lt;Option name=&quot;line_color&quot; type=&quot;QString&quot; value=&quot;60,60,60,255&quot;/&gt;&lt;Option name=&quot;line_width&quot; type=&quot;QString&quot; value=&quot;0.3&quot;/&gt;&lt;/Option&gt;&lt;/layer&gt;&lt;/symbol&gt;"/>'
        '<Option name="minLength" type="double" value="0"/>'
        '<Option name="minLengthMapUnitScale" type="QString" value="3x:0,0,0,0,0,0"/>'
        '<Option name="minLengthUnit" type="QString" value="MM"/>'
        '<Option name="offsetFromAnchor" type="double" value="0"/>'
        '<Option name="offsetFromAnchorMapUnitScale" type="QString" value="3x:0,0,0,0,0,0"/>'
        '<Option name="offsetFromAnchorUnit" type="QString" value="MM"/>'
        '<Option name="offsetFromLabel" type="double" value="0"/>'
        '<Option name="offsetFromLabelMapUnitScale" type="QString" value="3x:0,0,0,0,0,0"/>'
        '<Option name="offsetFromLabelUnit" type="QString" value="MM"/>'
        '</Option>\n'
        '      </callout>\n'
        '    </settings>\n'
        '  </labeling>\n'
        '  <blendMode>0</blendMode>\n'
        '  <featureBlendMode>0</featureBlendMode>\n'
        '  <layerOpacity>1</layerOpacity>\n'
        '</qgis>\n'
    )
    qml_path.write_text(content, encoding="utf-8")


def parse_orig_arg(s: str):
    parts = s.split(",", 2)
    if len(parts) < 2:
        raise argparse.ArgumentTypeError(f"--orig は lat,lon または lat,lon,name の形式で指定してください: {s!r}")
    lat  = float(parts[0])
    lon  = float(parts[1])
    name = parts[2].strip() if len(parts) == 3 else f"{lat:.4f},{lon:.4f}"
    return lat, lon, name


def load_orig_csv(path: str):
    origins = []
    with open(path, encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            lat  = float(row["lat"])
            lon  = float(row["lon"])
            name = row.get("name", f"{lat:.4f},{lon:.4f}").strip()
            origins.append((lat, lon, name))
    return origins


def build_multisource_graph(G: csr_matrix, orig_snap_infos: list):
    """
    スーパーソースノード（インデックス = G.shape[0]）を追加し、
    各始点ノードへのエッジ（重み = snap_time_min）を張った拡張グラフを返す。
    orig_snap_infos: [(node_idx, snap_time_min), ...]
    """
    from scipy.sparse import coo_matrix, vstack as sp_vstack

    n_v = G.shape[0]
    super_src_idx = n_v

    super_row_arr = np.zeros(len(orig_snap_infos), dtype=np.int32)  # shape=(1,*) なので行インデックスは0
    super_col_arr = np.array([idx for idx, _ in orig_snap_infos], dtype=np.int32)
    super_w_arr   = np.array([t for _, t in orig_snap_infos], dtype=np.float64)

    extra = coo_matrix(
        (super_w_arr, (super_row_arr, super_col_arr)),
        shape=(1, n_v + 1),
    ).tocsr()

    G_padded = csr_matrix(
        (G.data, G.indices, G.indptr),
        shape=(n_v, n_v + 1),
    )

    G_ext = sp_vstack([G_padded, extra], format="csr")
    return G_ext, super_src_idx


def main():
    ap = argparse.ArgumentParser(description="到達圏分析（複数始点対応）")
    ap.add_argument("--links",    default=DEFAULT_LINKS,  help="道路リンク parquet")
    ap.add_argument("--nodes",    default=DEFAULT_NODES,  help="道路ノード parquet")
    ap.add_argument("--access",   default=DEFAULT_ACCESS, help="L6アクセスリンク parquet")
    ap.add_argument("--orig",     action="append", default=[],
                    metavar="lat,lon[,name]",
                    help="始点（複数回指定可）。例: 35.8578,139.6490,埼玉県庁")
    ap.add_argument("--orig-csv", default=None,
                    help="始点CSVファイル（lat,lon,name 列）")
    ap.add_argument("--out-dir",  default=DEFAULT_OUT, help="出力ディレクトリ")
    args = ap.parse_args()

    origins = []
    if args.orig_csv:
        origins.extend(load_orig_csv(args.orig_csv))
    for s in args.orig:
        origins.append(parse_orig_arg(s))
    if not origins:
        origins = [(35.8578, 139.6490, "埼玉県庁")]

    t0      = time.time()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print("道路リンク読み込み中...")
    links = gpd.read_parquet(args.links)
    print(f"  {len(links):,} 本")

    print("道路ノード読み込み中...")
    nodes = gpd.read_parquet(args.nodes)
    print(f"  {len(nodes):,} 件")

    print("アクセスリンク読み込み中...")
    access = gpd.read_parquet(args.access)
    print(f"  {len(access):,} 件  ({time.time()-t0:.1f}s)")

    print("グラフ構築中...")
    _, n2i, G = build_sparse_graph(links)
    print(f"  ノード: {G.shape[0]:,}  エッジ: {G.nnz:,}  ({time.time()-t0:.1f}s)")

    snap_results = []
    for lat, lon, name in origins:
        nid, s_lat, s_lon, s_m = nearest_road_node(nodes, lat, lon)
        if nid not in n2i:
            raise ValueError(f"始点ノード {nid}（{name}）がグラフに存在しません")
        idx = n2i[nid]
        print(f"始点: {name}  ノード={nid}  ({s_lat:.5f},{s_lon:.5f})  snap={s_m:.0f}m")
        snap_results.append((lat, lon, name, nid, s_lat, s_lon, s_m, idx))

    print("\nメッシュ別到達時間計算中...")
    road_nids  = access["road_node"].astype(int).to_numpy()
    acc_time   = (access["time_001min"].astype(float) * 0.01).to_numpy()
    mesh_codes = access["mesh_code"].astype(str).to_numpy()
    graph_idxs = np.array([n2i.get(int(rn), -1) for rn in road_nids], dtype=np.int32)
    valid      = graph_idxs >= 0
    safe_idxs  = np.where(valid, graph_idxs, 0)

    # ── 始点ごとに個別 Dijkstra → 個別到達圏マップを出力 ──────
    print(f"\n到達時間マップ出力中...")
    for lat, lon, name, nid, s_lat, s_lon, s_m, orig_idx in snap_results:
        print(f"  Dijkstra: {name} ...")
        dist_road = sp_dijkstra(G, directed=True, indices=orig_idx)
        da_road   = np.where(valid, dist_road[safe_idxs], np.inf)
        dist_mesh = np.where(da_road < np.inf, da_road + acc_time, np.inf)
        save_arrival_map(
            mesh_codes, dist_mesh,
            out_dir / f"arrival_map_{name}.parquet",
            out_dir / f"arrival_map_{name}.qml",
        )
        print(f"    完了  ({time.time()-t0:.1f}s)")

    od_features = []
    for lat, lon, name, nid, s_lat, s_lon, s_m, _ in snap_results:
        od_features.append({
            "type": "Feature",
            "geometry": {"type": "Point", "coordinates": [lon, lat]},
            "properties": {
                "type": "origin",
                "name": name,
                "lat": lat,
                "lon": lon,
                "snap_lat": s_lat,
                "snap_lon": s_lon,
                "snap_m": round(s_m),
            },
        })

    origins_label = "_".join(name for _, _, name, *_ in snap_results)
    od_geojson = {"type": "FeatureCollection", "features": od_features}
    od_path = out_dir / f"origins_{origins_label}.geojson"
    with open(od_path, "w", encoding="utf-8") as f:
        json.dump(od_geojson, f, ensure_ascii=False, indent=2)
    od_qml_path = out_dir / f"origins_{origins_label}.qml"
    write_od_qml(od_qml_path)
    print(f"  {od_path.name}")
    print(f"  {od_qml_path.name}")

    print(f"\n総処理時間: {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
