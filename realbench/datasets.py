"""Dataset loaders, inclusion rules, cluster-level dev/test split and blind sampling.

Every loader turns raw files into Curve objects. Nothing here looks at any method's output.
Raw files live under REALBENCH_DATA (default: <repo>/data/raw); see data/README.md.
Anything a loader drops is counted in SKIPS and reported, never silently discarded.
"""

from __future__ import annotations

import collections
import glob
import hashlib
import heapq
import json
import math
import os
import re
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Iterator

import numpy as np

from . import config as C
from .geometry import arclength, dedupe, diameter, is_closed

ROOT = Path(__file__).resolve().parents[1]
EARTH_KM = 6371.0088
SKIPS: collections.Counter = collections.Counter()


def data_root() -> Path:
    return Path(os.environ.get("REALBENCH_DATA", ROOT / "data" / "raw"))


@dataclass
class Curve:
    dataset: str
    curve_id: str
    cluster: str
    pts: np.ndarray
    closed: bool
    meta: dict = field(default_factory=dict)

    @property
    def n(self) -> int:
        return len(self.pts)


# --------------------------------------------------------------------------- helpers

def _project(lonlat, scale_km: float) -> np.ndarray:
    a = np.asarray(lonlat, dtype=np.float64)[:, :2]
    lon0, lat0 = a[:, 0].mean(), a[:, 1].mean()
    x = np.radians(a[:, 0] - lon0) * math.cos(math.radians(lat0)) * scale_km
    y = np.radians(a[:, 1] - lat0) * scale_km
    return np.column_stack([x, y])


def _rel(path: str) -> str:
    try:
        return str(Path(path).resolve().relative_to(data_root().resolve()))
    except ValueError:
        return os.path.basename(path)


def _cell(lon: float, lat: float) -> str:
    return f"cell{math.floor(lon / C.MAP_CELL_DEG)},{math.floor(lat / C.MAP_CELL_DEG)}"


# --------------------------------------------------------------------------- maps

def load_geojson(pattern: str, group_fields: tuple = ()):
    """GeoJSON (Multi)LineStrings, projected to km.

    Cluster (protocol 3.2): features sharing the same non-empty values of
    ``group_fields`` (e.g. a river's name) form one group; each group is assigned to
    the 10 x 10 degree cell of its centroid, and the cell is the cluster. Nearby and
    same-named features therefore always share a cluster and a split.
    """
    def gen():
        for fp in sorted(glob.glob(str(data_root() / pattern))):
            with open(fp) as f:
                gj = json.load(f)
            parts, groups = [], collections.defaultdict(list)
            for fi, ft in enumerate(gj["features"]):
                g = ft.get("geometry")
                if not g:
                    continue
                props = ft.get("properties") or {}
                vals = tuple(props.get(k) for k in group_fields)
                key = ("g",) + vals if group_fields and all(v not in (None, "") for v in vals) else ("f", fi)
                coords = [g["coordinates"]] if g["type"] == "LineString" else (
                    g["coordinates"] if g["type"] == "MultiLineString" else [])
                for pi, c in enumerate(coords):
                    if len(c) >= 2:
                        parts.append((fi, pi, c, key, props.get("name")))
                        groups[key].append(np.asarray(c, float)[:, :2])
            cell = {k: _cell(*np.vstack(v).mean(0)) for k, v in groups.items()}
            for fi, pi, c, key, name in parts:
                pts = _project(c, EARTH_KM)
                yield Curve("", f"{_rel(fp)}#{fi}.{pi}", f"{_rel(fp)}#{cell[key]}", pts, False, {"name": name})
    return gen


# --------------------------------------------------------------------------- GPS

def load_gpx(pattern: str = "gpx/**/*.gpx"):
    """GPX track segments only (routes and waypoints are plans, not recordings). Metres."""
    def gen():
        import gpxpy
        for fp in sorted(glob.glob(str(data_root() / pattern), recursive=True)):
            try:
                with open(fp, "rb") as f:
                    g = gpxpy.parse(f)
            except Exception:
                SKIPS["gpx parse failure"] += 1
                continue
            SKIPS["gpx routes ignored"] += len(g.routes)
            for ti, t in enumerate(g.tracks):
                for si, s in enumerate(t.segments):
                    if len(s.points) < 2:
                        continue
                    pts = _project([(q.longitude, q.latitude) for q in s.points], EARTH_KM * 1000.0)
                    yield Curve("", f"{_rel(fp)}#t{ti}s{si}", _rel(fp), pts, False)
    return gen


def load_geolife(pattern: str = "geolife/**/Trajectory/*.plt"):
    """Microsoft GeoLife Trajectories 1.3: one .plt file = one recorded trajectory. Metres.

    Format: 6 header lines, then lat,lon,0,alt,days,date,time. Cluster = user (the
    directory above Trajectory/), so all of one person's trips share a split.
    """
    def gen():
        for fp in sorted(glob.glob(str(data_root() / pattern), recursive=True)):
            user = Path(fp).parent.parent.name
            try:
                rows = [l.split(",") for l in open(fp).read().splitlines()[6:] if l.strip()]
                ll = np.array([(float(r[1]), float(r[0])) for r in rows])
            except Exception:
                SKIPS["geolife parse failure"] += 1
                continue
            if len(ll) >= 2:
                yield Curve("", f"{user}/{Path(fp).stem}", user, _project(ll, EARTH_KM * 1000.0), False)
    return gen


# --------------------------------------------------------------------------- pen strokes

def load_quickdraw_raw(pattern: str = "quickdraw/*.ndjson"):
    """Quick Draw RAW ndjson: drawing = [[x...], [y...], [t...]] per stroke.

    Only the pre-registered categories (config.QUICKDRAW_CATEGORIES) are read; other
    files are ignored and counted. Simplified files (2 arrays per stroke) were produced
    with RDP and are refused. Cluster = drawing key_id.
    """
    allowed = set(C.QUICKDRAW_CATEGORIES)

    def gen():
        for fp in sorted(glob.glob(str(data_root() / pattern))):
            cat = Path(fp).stem.replace("full_raw_", "").replace("_", " ")
            if cat not in allowed:
                SKIPS["quickdraw file not in pre-registered categories"] += 1
                continue
            with open(fp) as f:
                for line in f:
                    d = json.loads(line)
                    for si, st in enumerate(d["drawing"]):
                        if len(st) < 3:
                            raise ValueError(f"{fp}: stroke has no timing array; this is a "
                                             "'simplified' (RDP-processed) file, excluded by protocol 3.3")
                        pts = np.column_stack([st[0], st[1]]).astype(np.float64)
                        yield Curve("", f"{d['word']}/{d['key_id']}#{si}", f"{d['word']}/{d['key_id']}", pts, False,
                                    {"word": d["word"]})
    return gen


def load_omniglot(pattern: str = "omniglot/strokes_*/*/*/*.txt"):
    """Omniglot raw pen strokes (Lake et al. 2015, MIT licence): 1623 characters x 20 drawers.

    File = one drawing of one character: ``START``, then ``x,y,t`` rows, ``BREAK`` between
    pen strokes. Each pen stroke is one open polyline (time is ignored). Both the
    background and evaluation sets are read whole. Cluster = alphabet/character, so all 20
    drawings of a character share a split. Coordinates are snapped to whole pixels by the
    collection tool, so these strokes carry quantisation noise by construction.
    """
    def gen():
        for fp in sorted(glob.glob(str(data_root() / pattern))):
            p = Path(fp)
            char = f"{p.parent.parent.name}/{p.parent.name}"
            strokes, cur = [], []
            for line in p.read_text().splitlines():
                line = line.strip()
                if not line:
                    continue
                if line in ("START", "BREAK"):
                    if cur:
                        strokes.append(cur)
                    cur = []
                    continue
                v = line.split(",")
                try:
                    cur.append((float(v[0]), float(v[1])))
                except (ValueError, IndexError):
                    SKIPS["omniglot malformed line"] += 1
            if cur:
                strokes.append(cur)
            for si, st in enumerate(strokes):
                if len(st) >= 2:
                    yield Curve("", f"{char}/{p.stem}#{si}", char, np.array(st, dtype=np.float64), False)
    return gen


def load_inkml(pattern: str = "inkml/**/*.inkml"):
    """InkML traces (CROHME, MathWriting). X/Y channels read from <traceFormat> when given.

    Cluster = writer annotation when present (CROHME writers span many files), else the file.
    """
    def gen():
        for fp in sorted(glob.glob(str(data_root() / pattern), recursive=True)):
            try:
                root = ET.parse(fp).getroot()
            except Exception:
                SKIPS["inkml parse failure"] += 1
                continue
            ns = root.tag.split("}")[0] + "}" if root.tag.startswith("{") else ""
            writer = next((a.text.strip() for a in root.iter(f"{ns}annotation")
                           if a.get("type") == "writer" and a.text and a.text.strip()), None)
            chans = [c.get("name") for c in root.iter(f"{ns}channel")]
            ix = chans.index("X") if "X" in chans else 0
            iy = chans.index("Y") if "Y" in chans else 1
            cluster = f"writer:{writer}" if writer else _rel(fp)
            for ti, tr in enumerate(root.iter(f"{ns}trace")):
                pts = []
                for tok in (tr.text or "").strip().split(","):
                    v = tok.split()
                    if len(v) > max(ix, iy):
                        try:
                            pts.append((float(v[ix]), float(v[iy])))
                        except ValueError:
                            SKIPS["inkml non-numeric point"] += 1
                if len(pts) >= 2:
                    yield Curve("", f"{_rel(fp)}#{tr.get('id', ti)}", cluster, np.array(pts), False)
    return gen


def load_uji(pattern: str = "uji/*.txt"):
    """UJI Pen Characters v2 ("WORD <char> <session>", "NUMSTROKES n", "POINTS n # x y ...").

    Curve id = session/char/occurrence#stroke. Cluster = writer (W.. token of the session).
    """
    def gen():
        for fp in sorted(glob.glob(str(data_root() / pattern))):
            occ = collections.Counter()
            cid, k, writer = None, 0, None
            with open(fp, errors="replace") as f:
                for line in f:
                    t = line.strip()
                    if t.startswith("WORD"):
                        parts = t.split()
                        char, session = (parts[1], parts[2]) if len(parts) >= 3 else (parts[-1], parts[-1])
                        occ[(session, char)] += 1
                        cid = f"{session}/{char}/{occ[(session, char)]}"
                        m = re.search(r"W\d+", session)
                        writer = m.group(0) if m else session
                        k = 0
                    elif t.startswith("POINTS") and cid is not None:
                        nums = t.split("#", 1)[1].split()
                        xy = np.array(nums, dtype=np.float64).reshape(-1, 2)
                        yield Curve("", f"{cid}#{k}", writer, xy, False)
                        k += 1
    return gen


# --------------------------------------------------------------------------- CAD

def _uniform_with_corners(primitives) -> np.ndarray:
    """primitives: list of (f(t in [0,1]) -> (n,2), length).

    Samples each primitive with spacing <= loop length / CAD_SAMPLES_PER_LOOP and keeps
    every primitive endpoint (corners), per protocol 3.3.
    """
    total = sum(L for _, L in primitives)
    step = total / C.CAD_SAMPLES_PER_LOOP
    out = []
    for f, L in primitives:
        n = max(1, int(math.ceil(L / step)))
        seg = f(np.linspace(0.0, 1.0, n + 1))
        out.append(seg if not out else seg[1:])
    return np.vstack(out)


def _line(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    return (lambda t: a + np.outer(t, b - a)), float(np.hypot(*(b - a)))


def _arc(center, start, end, ccw: bool, full: bool = False):
    c = np.asarray(center, float)
    a0 = math.atan2(start[1] - c[1], start[0] - c[0])
    a1 = math.atan2(end[1] - c[1], end[0] - c[0])
    r = float(np.hypot(start[0] - c[0], start[1] - c[1]))
    if full:
        sweep = 2 * math.pi
    else:
        sweep = (a1 - a0) % (2 * math.pi)
        if not ccw:
            sweep = sweep - 2 * math.pi
        if abs(sweep) < 1e-12:
            sweep = 2 * math.pi if ccw else -2 * math.pi
    return (lambda t: c + r * np.column_stack([np.cos(a0 + sweep * t), np.sin(a0 + sweep * t)])), abs(sweep) * r


def _xy(p):
    return (float(p["x"]), float(p["y"]))


def _close(pts: np.ndarray) -> np.ndarray:
    pts = pts.copy()
    pts[-1] = pts[0]
    return pts


def load_fusion360(pattern: str = "fusion360/**/*.json"):
    """Fusion 360 Gallery reconstruction JSON: each closed profile loop is one curve.

    Verified on Line3D loops (repo test files). Arc3D / Circle3D follow the Fusion API
    field names (center_point, start_point, end_point, normal, radius) and are checked for
    consistency (both ends at the same radius). Loops containing any other curve type
    (splines, ellipses, ...) cannot be sampled and are skipped; skips are counted by type
    and reported (a stated threat to validity, protocol 7). Cluster = project id.
    """
    dist = lambda p, q: math.hypot(p[0] - q[0], p[1] - q[1])

    def gen():
        for fp in sorted(glob.glob(str(data_root() / pattern), recursive=True)):
            try:
                d = json.load(open(fp))
            except Exception:
                SKIPS["fusion360 parse failure"] += 1
                continue
            project = Path(fp).stem.split("_")[0]
            for eid, e in (d.get("entities") or {}).items():
                if e.get("type") != "Sketch":
                    continue
                for pid, pr in (e.get("profiles") or {}).items():
                    for li, lp in enumerate(pr.get("loops", [])):
                        prims, bad = [], None
                        for c in lp.get("profile_curves", []):
                            t = c.get("type")
                            try:
                                if t == "Line3D":
                                    prims.append(("seg", _xy(c["start_point"]), _xy(c["end_point"]), c))
                                elif t == "Arc3D":
                                    a, b, ctr = _xy(c["start_point"]), _xy(c["end_point"]), _xy(c["center_point"])
                                    if abs(dist(a, ctr) - dist(b, ctr)) > 1e-6 * max(dist(a, ctr), 1e-12):
                                        bad = "Arc3D inconsistent radius"
                                        break
                                    prims.append(("arc", a, b, c))
                                elif t == "Circle3D":
                                    prims.append(("circle", None, None, c))
                                else:
                                    bad = f"loop with {t}"
                                    break
                            except (KeyError, TypeError):
                                bad = f"{t} missing fields"
                                break
                        if bad or not prims:
                            SKIPS[f"fusion360 skipped: {bad or 'empty loop'}"] += 1
                            continue
                        built = []
                        if len(prims) == 1 and prims[0][0] == "circle":
                            c = prims[0][3]
                            ctr = _xy(c["center_point"])
                            r = float(c["radius"])
                            built.append(_arc(ctr, (ctr[0] + r, ctr[1]), (ctr[0] + r, ctr[1]), True, full=True))
                        elif any(k == "circle" for k, *_ in prims):
                            SKIPS["fusion360 skipped: circle inside multi-curve loop"] += 1
                            continue
                        else:
                            # orient the first primitive toward the second, then chain
                            pos = None
                            if len(prims) > 1:
                                _, a0, b0, _ = prims[0]
                                _, a1, b1, _ = prims[1]
                                b0_joins = min(dist(b0, a1), dist(b0, b1)) <= min(dist(a0, a1), dist(a0, b1))
                                pos = a0 if b0_joins else b0
                            start, gap = None, 0.0
                            for kind, a, b, c in prims:
                                rev = pos is not None and dist(a, pos) > dist(b, pos)
                                if rev:
                                    a, b = b, a
                                if start is None:
                                    start = a
                                else:
                                    gap = max(gap, dist(a, pos))
                                if kind == "seg":
                                    built.append(_line(a, b))
                                else:
                                    nz = float((c.get("normal") or {}).get("z", 1.0))
                                    built.append(_arc(_xy(c["center_point"]), a, b, (nz >= 0) != rev))
                                pos = b
                            scale = max(L for _, L in built)
                            if gap > C.CLOSED_TOL * scale or dist(pos, start) > C.CLOSED_TOL * scale:
                                SKIPS["fusion360 skipped: loop does not close"] += 1
                                continue
                        pts = _close(_uniform_with_corners(built))
                        yield Curve("", f"{Path(fp).stem}:{eid}:{pid}:{li}", project, pts, True)
    return gen


def load_dxf(pattern: str = "dxf/**/*.dxf"):
    """DXF via ezdxf. Polylines and block INSERTs are exploded into their LINE/ARC/... pieces;
    every piece is flattened finely on its own and pieces are chained end to end (growing at
    both ends). Each chain is resampled at uniform arc length (length / 1000) with every
    piece boundary kept as a corner. UNVERIFIED on a real corpus. Cluster = file.
    """
    def gen():
        import ezdxf
        from ezdxf import path as ezpath

        def explode(e, depth=0):
            t = e.dxftype()
            if t in ("LWPOLYLINE", "POLYLINE", "INSERT") and depth < 8:
                try:
                    for v in e.virtual_entities():
                        yield from explode(v, depth + 1)
                except Exception:
                    SKIPS["dxf explode failure"] += 1
            elif t in ("LINE", "ARC", "SPLINE", "CIRCLE", "ELLIPSE"):
                yield e

        for fp in sorted(glob.glob(str(data_root() / pattern), recursive=True)):
            try:
                doc = ezdxf.readfile(fp)
            except Exception:
                SKIPS["dxf read failure"] += 1
                continue
            pieces = []
            for e in doc.modelspace():
                for v in explode(e):
                    try:
                        p = ezpath.make_path(v)
                        dense = np.array([(q.x, q.y) for q in p.flattening(1e-9 + 1e-6 * max(
                            1.0, *(abs(c) for q in p.control_vertices() for c in (q.x, q.y))))])
                    except Exception:
                        SKIPS["dxf flatten failure"] += 1
                        continue
                    dense = dedupe(dense)
                    if len(dense) >= 2:
                        pieces.append(dense)
            loose = [d for d in pieces if not is_closed(d)]
            chains = [[d] for d in pieces if is_closed(d)]
            while loose:
                ch = [loose.pop(0)]
                tol = lambda: 1e-6 * max(diameter(np.vstack(ch)), 1e-9)
                grown = True
                while grown:
                    grown = False
                    head, tail = ch[0][0], ch[-1][-1]
                    for i, d in enumerate(loose):
                        if np.hypot(*(d[0] - tail)) <= tol():
                            ch.append(loose.pop(i)); grown = True; break
                        if np.hypot(*(d[-1] - tail)) <= tol():
                            ch.append(loose.pop(i)[::-1]); grown = True; break
                        if np.hypot(*(d[-1] - head)) <= tol():
                            ch.insert(0, loose.pop(i)); grown = True; break
                        if np.hypot(*(d[0] - head)) <= tol():
                            ch.insert(0, loose.pop(i)[::-1]); grown = True; break
                chains.append(ch)
            for ci, parts in enumerate(chains):
                prims = []
                for d in parts:
                    s = arclength(d)
                    if s[-1] > 0:
                        prims.append(((lambda t, d=d, s=s: np.column_stack([
                            np.interp(t * s[-1], s, d[:, 0]), np.interp(t * s[-1], s, d[:, 1])])), float(s[-1])))
                if prims:
                    pts = _uniform_with_corners(prims)
                    yield Curve("", f"{_rel(fp)}#{ci}", _rel(fp), pts, is_closed(pts, C.CLOSED_TOL))
    return gen


# --------------------------------------------------------------------------- registry

@dataclass
class Dataset:
    id: str
    category: str
    loader: Callable[[], Iterator[Curve]]
    description: str
    files: str        # glob (relative to data root) of the raw files this dataset reads


DATASETS = {d.id: d for d in [
    Dataset("ne_coastline", "maps", load_geojson("natural_earth/ne_10m_coastline.geojson"),
            "Natural Earth 1:10m coastline", "natural_earth/ne_10m_coastline.geojson"),
    Dataset("ne_rivers", "maps", load_geojson("natural_earth/ne_10m_rivers_lake_centerlines.geojson", ("name",)),
            "Natural Earth 1:10m river centerlines", "natural_earth/ne_10m_rivers_lake_centerlines.geojson"),
    Dataset("ne_roads", "maps", load_geojson("natural_earth/ne_10m_roads.geojson", ("sov_a3", "name")),
            "Natural Earth 1:10m roads", "natural_earth/ne_10m_roads.geojson"),
    Dataset("geolife", "GPS", load_geolife(), "Microsoft GeoLife Trajectories 1.3 (metres)",
            "geolife/**/Trajectory/*.plt"),
    Dataset("gpx", "GPS", load_gpx(), "GPX track segments (metres); exploratory only", "gpx/**/*.gpx"),
    Dataset("quickdraw_raw", "pen strokes", load_quickdraw_raw(), "Google Quick Draw raw strokes, 10 fixed categories",
            "quickdraw/*.ndjson"),
    Dataset("omniglot", "pen strokes", load_omniglot(), "Omniglot raw pen strokes, all 50 alphabets",
            "omniglot/strokes_*/*/*/*.txt"),
    Dataset("inkml", "pen strokes", load_inkml(), "InkML handwriting traces (CROHME / MathWriting)",
            "inkml/**/*.inkml"),
    Dataset("uji", "pen strokes", load_uji(), "UJI Pen Characters v2 strokes", "uji/*.txt"),
    Dataset("fusion360", "CAD", load_fusion360(), "Fusion 360 Gallery sketch profile loops",
            "fusion360/**/*.json"),
    Dataset("dxf", "CAD", load_dxf(), "DXF drawings (lines, arcs, polylines, splines)", "dxf/**/*.dxf"),
]}


# --------------------------------------------------------------------------- split and sampling

def _unit_hash(*parts) -> float:
    h = hashlib.sha256(":".join(str(p) for p in parts).encode()).digest()
    return int.from_bytes(h, "big") / 2 ** 256


def split_of(dataset: str, cluster: str) -> str:
    return "dev" if _unit_hash(C.SPLIT_SALT, dataset, cluster) < C.DEV_FRACTION else "test"


def _content_key(pts: np.ndarray, closed: bool) -> str:
    """Duplicate key invariant to position, scale, traversal direction and (closed) start point."""
    d = max(diameter(pts), 1e-300)
    q = np.round((pts - pts.min(0)) / d * 1e9).astype(np.int64)
    cands = []
    for arr in (q, q[::-1]):
        if closed:
            ring = arr[:-1]
            k = np.lexsort((ring[:, 1], ring[:, 0]))[0]
            arr = np.vstack([np.roll(ring, -k, axis=0), ring[k:k + 1]])
        cands.append(hashlib.sha1(np.ascontiguousarray(arr).tobytes()).hexdigest())
    return min(cands)


def raw_files(dataset: str) -> list[str]:
    return sorted(glob.glob(str(data_root() / DATASETS[dataset].files), recursive=True))


def load_split(dataset: str, split: str, cap: int | None = None, stats: dict | None = None) -> list[Curve]:
    """Eligible curves of one split, as a blind simple random sample of at most ``cap``.

    The sample is the ``cap`` curves with the smallest sha256(seed, dataset, curve_id),
    which is a uniform random sample that can be drawn while streaming large files.
    Exact duplicates (any position, scale, direction or start point) are removed across
    the whole dataset before splitting (first in file order kept), so a duplicate can
    never sit in both dev and test. Closed curves are snapped exactly shut.
    """
    ds = DATASETS[dataset]
    cap = C.CAP[split] if cap is None else cap
    SKIPS.clear()
    st = {"raw_curves": 0, "too_short": 0, "duplicates": 0, "eligible": {"dev": 0, "test": 0}}
    clusters = {"dev": set(), "test": set()}
    seen: set = set()
    ids_seen: set = set()
    heap: list = []
    for cv in ds.loader():
        st["raw_curves"] += 1
        pts = dedupe(cv.pts)
        if len(pts) < C.MIN_VERTICES or not np.isfinite(pts).all() or arclength(pts)[-1] <= 0:
            st["too_short"] += 1
            continue
        closed = cv.closed or is_closed(pts, C.CLOSED_TOL)
        if closed:
            pts = _close(pts)
            if len(pts) >= 3 and np.array_equal(pts[-2], pts[-1]):
                pts = pts[:-1]
                pts[-1] = pts[0]
        if cv.curve_id in ids_seen:
            raise ValueError(f"{dataset}: duplicate curve id {cv.curve_id!r} from the loader")
        ids_seen.add(cv.curve_id)
        key = _content_key(pts, closed)
        if key in seen:
            st["duplicates"] += 1
            continue
        seen.add(key)
        sp = split_of(dataset, cv.cluster)
        st["eligible"][sp] += 1
        clusters[sp].add(cv.cluster)
        if sp != split or cap <= 0:
            continue
        cv = Curve(dataset, cv.curve_id, cv.cluster, pts, closed, cv.meta)
        rank = _unit_hash(C.SEED, dataset, cv.curve_id)
        if len(heap) < cap:
            heapq.heappush(heap, (-rank, cv.curve_id, cv))
        elif -heap[0][0] > rank:
            heapq.heapreplace(heap, (-rank, cv.curve_id, cv))
    st["clusters"] = {k: len(v) for k, v in clusters.items()}
    st["skips"] = dict(SKIPS)
    if stats is not None:
        stats.update(st)
        stats["_cluster_ids"] = clusters
    return sorted((c for _, _, c in heap), key=lambda c: c.curve_id)
