# 到達圏分析マップ

国土数値情報（KSJ）道路データを使い、任意の始点からの到達時間を 125m メッシュで可視化する Web マップ。

`input/origins.csv` を編集してコミットするだけで GitHub Actions が自動的に分析を実行し、GitHub Pages に公開する。

**公開 URL**: https://shiwaku.github.io/ksj-reachability-analysis-map/

---

## デモ（埼玉県・3始点）

| 始点 | 概要 |
|---|---|
| 川越市役所 | 埼玉県中部 |
| 東松山市役所 | 埼玉県中部 |
| 大宮駅 | 埼玉県南部 |

---

## 仕組み

```
input/origins.csv を編集してコミット
    ↓ GitHub Actions がトリガー
    Dijkstra で各始点からの到達時間を計算（約 10 秒 / 始点）
    → GeoJSON → tippecanoe で PMTiles に変換
    → MapLibre GL JS ビューワーと一緒に GitHub Pages にデプロイ
```

---

## 使い方

### 始点を変更する

`input/origins.csv` を編集してコミットするだけで自動更新される。

```csv
lat,lon,name
35.9249,139.4858,川越市役所
36.0420,139.4006,東松山市役所
35.9063,139.6239,大宮駅
```

緯度・経度は Google マップで地点を右クリックするとコピーできる。

### 手動で Actions を実行する

GitHub の **Actions タブ → Compute & Deploy → Run workflow** から手動実行できる。

---

## ディレクトリ構成

```
.github/workflows/
  compute.yml           GitHub Actions ワークフロー

src/
  reachability_search.py  到達圏分析（Dijkstra）
  to_geojson.py           parquet → GeoJSON 変換（tippecanoe 入力用）
  ksj_to_network_csv.py   国土数値情報 GeoJSON → 道路ネットワーク変換
  make_access_links.py    L6 アクセスリンク生成

viewer/
  index.html            MapLibre GL JS ビューワー

data/
  prefecture.parquet    都道府県境界（--pref クリップ用）
  city.parquet          市区町村境界（--city クリップ用）

network/
  saitama_pref/         埼玉県サンプルネットワーク（リポジトリに同梱）

input/
  origins.csv           始点リスト（lat,lon,name）← ここを編集
```

---

## 新しいエリアで使う

### 1. 国土数値情報のダウンロード

[国土数値情報ダウンロードサービス](https://nlftp.mlit.go.jp/ksj/gml/datalist/KsjTmplt-N13-v2_1.html) から道路データ（N13-24）の GeoJSON をダウンロードし、`input/` に配置する。

```
input/
├── N13-24_5338/
│   └── N13-24_5338.geojson
└── N13-24_5339/
    └── N13-24_5339.geojson
```

### 2. ネットワーク生成

```bash
python3 src/ksj_to_network_csv.py \
  --meshes 5338,5339,5438,5439 \
  --case saitama_pref \
  --pref 埼玉県

python3 src/make_access_links.py \
  --meshes 5338,5339,5438,5439 \
  --case saitama_pref \
  --level 6 \
  --pref 埼玉県
```

### 3. 到達圏分析を実行（ローカル確認用）

```bash
python3 src/reachability_search.py --orig-csv input/origins.csv \
  --links network/saitama_pref/KSJ_N13-24_saitama_pref_道路リンク.parquet \
  --nodes network/saitama_pref/KSJ_N13-24_saitama_pref_道路ノード.parquet \
  --access network/saitama_pref/KSJ_N13-24_saitama_pref_アクセスリンク_L6.parquet
```

### 4. ネットワークをリポジトリに追加

`.gitignore` の除外対象になっているため、新しいエリアのネットワークを含める場合は `.gitignore` に例外を追記する。

```gitignore
!network/{case}/
!network/{case}/**
```

---

## 主要都道府県の 1 次メッシュコード

| 都道府県 | 1 次メッシュコード |
|---|---|
| 埼玉県 | 5338, 5339, 5438, 5439 |
| 東京都 | 5338, 5339, 5438, 5439 |
| 神奈川県 | 5238, 5239, 5338, 5339 |
| 千葉県 | 5239, 5240, 5339, 5340, 5439, 5440 |
| 愛知県 | 5236, 5237, 5336, 5337, 5436, 5437 |
| 大阪府 | 5135, 5235 |
| 福岡県 | 4930, 5030, 5031, 5032 |

---

## 必要環境（ローカル実行時）

```bash
pip install geopandas pyarrow scipy
```

---

## 制約・注意事項

- **一方通行未考慮**: 国土数値情報には一方通行フィールドが存在しないため、全道路を双方向リンクとして扱っている
- **メッシュレベル**: L6（125m）固定。変更する場合は `make_access_links.py --level` とネットワーク再生成が必要
- **速度モデル**: vehicle モード（道路種別・幅員による速度テーブル）

---

## データについて

- **出典**: 国土数値情報 道路データ（N13）/ 国土交通省
- **ダウンロード**: https://nlftp.mlit.go.jp/ksj/gml/datalist/KsjTmplt-N13-v2_1.html
- **複製承認**: 数値地図（国土基本情報）/ 国土地理院長承認（複製）R 6JHf503
