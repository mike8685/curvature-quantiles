"""Dataset loaders, inclusion rules, cluster-level dev/test split and blind sampling.

Every loader turns raw files into Curve objects. Nothing here looks at any method's output.
Raw files live under REALBENCH_DATA (default: <repo>/data/raw); see data/README.md.
"""

from __future__ import annotations

import glob
import hashlib
import heapq
import json
import math
import os
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Iterator

import numpy as np

from . import config as C
from .geometry import arclength, dedupe, diameter, is_closed

ROOT = Path(__file__).resolve().parents[1]
EARTH_KM = 6371.0088


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


# --------------------------------------------------------------------------- maps

def load_geojson(pattern: str, units: str = "km") -> Callable[[], Iterator[Curve]]:
    def gen():
        for fp in sorted(glob.glob(str(data_root() / pattern))):
            with open(fp) as f:
                gj = json.load(f)
            for fi, ft in enumerate(gj["features"]):
                g = ft.get("geometry")
                if not g:
                    continue
                parts = [g["coordinates"]] if g["type"] == "LineString" else (
                    g["coordinates"] if g["type"] == "MultiLineString" else [])
                for pi, c in enumerate(parts):
                    if len(c) < 2:
                        continue
                    pts = _project(c, EARTH_KM * (1000.0 if units == "m" else 1.0))
                    yield Curve("", f"{_rel(fp)}#{fi}.{pi}", f"{_rel(fp)}#{fi}", pts, False,
                                {"name": (ft.get("properties") or {}).get("name")})
    return gen


# --------------------------------------------------------------------------- GPS

def load_gpx(pattern: str = "gpx/**/*.gpx"):
    def gen():
        import gpxpy
        for fp in sorted(glob.glob(str(data_root() / pattern), recursive=True)):
            try:
                with open(fp, "rb") as f:
                    g = gpxpy.parse(f)
            except Exception:
                continue
            segs = [(f"t{ti}s{si}", s.points) for ti, t in enumerate(g.tracks) for si, s in enumerate(t.segments)]
            segs += [(f"r{ri}", r.points) for ri, r in enumerate(g.routes)]
            for tag, p in segs:
                if len(p) < 2:
                    continue
                pts = _project([(q.longitude, q.latitude) for q in p], EARTH_KM * 1000.0)
                yield Curve("", f"{_rel(fp)}#{tag}", _rel(fp), pts, False)
    return gen


# --------------------------------------------------------------------------- pen strokes

def load_quickdraw_raw(pattern: str = "quickdraw/*.ndjson"):
    """Quick Draw RAW ndjson: drawing = [[x...], [y...], [t...]] per stroke.

    Simplified files (2 arrays per stroke) were produced with RDP and are refused.
    Cluster = drawing key_id.
    """
    def gen():
        for fp in sorted(glob.glob(str(data_root() / pattern))):
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


def load_inkml(pattern: str = "inkml/**/*.inkml"):
    """InkML traces (CROHME, MathWriting). Cluster = file (one expression by one writer)."""
    def gen():
        for fp in sorted(glob.glob(str(data_root() / pattern), recursive=True)):
            root = ET.parse(fp).getroot()
            ns = root.tag.split("}")[0] + "}" if root.tag.startswith("{") else ""
            for ti, tr in enumerate(root.iter(f"{ns}trace")):
                pts = []
                for tok in (tr.text or "").strip().split(","):
                    v = tok.split()
                    if len(v) >= 2:
                        pts.append((float(v[0]), float(v[1])))
                if len(pts) >= 2:
                    yield Curve("", f"{_rel(fp)}#{tr.get('id', ti)}", _rel(fp), np.array(pts), False)
    return gen


def load_uji(pattern: str = "uji/*.txt"):
    """UJI Pen Characters v2 text format (WORD / NUMSTROKES / POINTS n # x y ...).

    UNVERIFIED against the real file: parser follows the published format description.
    `python -m realbench datasets` reports how many strokes parsed. Cluster = writer.
    """
    def gen():
        for fp in sorted(glob.glob(str(data_root() / pattern))):
            word = None
            k = 0
            with open(fp, errors="replace") as f:
                for line in f:
                    t = line.strip()
                    if t.startswith("WORD"):
                        parts = t.split()
                        word = parts[2] if len(parts) > 2 else parts[-1]
                        k = 0
                    elif t.startswith("POINTS") and word is not None:
                        nums = t.split("#", 1)[1].split()
                        xy = np.array(nums, dtype=np.float64).reshape(-1, 2)
                        writer = next((p for p in word.replace("_", "-").split("-") if p.startswith("W")), word)
                        yield Curve("", f"{word}#{k}", writer, xy, False)
                        k += 1
    return gen


# --------------------------------------------------------------------------- CAD

def _uniform_with_corners(primitives) -> np.ndarray:
    """primitives: list of callables f(t in [0,1]) -> (n,2) and their lengths.

    Samples each primitive with spacing <= loop length / CAD_SAMPLES_PER_LOOP and keeps
    every primitive endpoint (corners), per protocol 3.3.
    """
    total = sum(L for _, L in primitives)
    step = total / C.CAD_SAMPLES_PER_LOOP
    out = []
    for f, L in primitives:
        n = max(1, int(math.ceil(L / step)))
        t = np.linspace(0.0, 1.0, n + 1)
        seg = f(t)
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


def load_fusion360(pattern: str = "fusion360/**/*.json"):
    """Fusion 360 Gallery reconstruction JSON: each closed profile loop is one curve.

    Verified on Line3D loops (repo test files). Arc3D / Circle3D follow the Fusion API
    field names (center_point, start_point, end_point, normal, radius); loops with any
    other curve type are skipped and counted under meta 'skipped_types'.
    Cluster = project id (first field of the file name).
    """
    def gen():
        for fp in sorted(glob.glob(str(data_root() / pattern), recursive=True)):
            try:
                d = json.load(open(fp))
            except Exception:
                continue
            ents = d.get("entities", {})
            project = Path(fp).stem.split("_")[0]
            for eid, e in ents.items():
                if e.get("type") != "Sketch":
                    continue
                for pid, pr in (e.get("profiles") or {}).items():
                    for li, lp in enumerate(pr.get("loops", [])):
                        prims, ok, cur = [], True, None
                        for c in lp.get("profile_curves", []):
                            t = c.get("type")
                            if t == "Line3D":
                                a, b = _xy(c["start_point"]), _xy(c["end_point"])
                                prims.append(("seg", a, b, c))
                            elif t == "Arc3D":
                                a, b = _xy(c["start_point"]), _xy(c["end_point"])
                                prims.append(("arc", a, b, c))
                            elif t == "Circle3D":
                                prims.append(("circle", None, None, c))
                            else:
                                ok = False
                                break
                        if not ok or not prims:
                            continue
                        built = []
                        if len(prims) == 1 and prims[0][0] == "circle":
                            c = prims[0][3]
                            ctr = _xy(c["center_point"])
                            r = float(c["radius"])
                            built.append(_arc(ctr, (ctr[0] + r, ctr[1]), (ctr[0] + r, ctr[1]), True, full=True))
                        else:
                            # chain: orient each primitive so its start meets the previous end;
                            # the first one is oriented toward the second
                            if any(k == "circle" for k, *_ in prims):
                                continue
                            dist = lambda p, q: math.hypot(p[0] - q[0], p[1] - q[1])
                            pos = None
                            if len(prims) > 1:
                                _, a0, b0, _ = prims[0]
                                _, a1, b1, _ = prims[1]
                                pos = b0 if min(dist(b0, a1), dist(b0, b1)) <= min(dist(a0, a1), dist(a0, b1)) else a0
                                pos = a0 if pos is b0 else b0          # start point of the first primitive
                            start, gap = None, 0.0
                            for kind, a, b, c in prims:
                                if pos is not None and dist(a, pos) > dist(b, pos):
                                    a, b, rev = b, a, True
                                else:
                                    rev = False
                                if start is None:
                                    start = a
                                elif pos is not None:
                                    gap = max(gap, dist(a, pos))
                                if kind == "seg":
                                    built.append(_line(a, b))
                                else:
                                    nz = float((c.get("normal") or {}).get("z", 1.0))
                                    ccw = (nz >= 0) != rev
                                    built.append(_arc(_xy(c["center_point"]), a, b, ccw))
                                pos = b
                            scale = max(L for _, L in built)
                            if gap > 1e-6 * scale or dist(pos, start) > 1e-6 * scale:
                                continue                       # not a connected closed loop
                        pts = _uniform_with_corners(built)
                        yield Curve("", f"{Path(fp).stem}:{eid}:{pid}:{li}", project, pts, True)
    return gen


def load_dxf(pattern: str = "dxf/**/*.dxf"):
    """DXF via ezdxf: each LWPOLYLINE/POLYLINE/SPLINE/CIRCLE/ELLIPSE entity is a curve;
    LINE/ARC entities are chained end to end. Entities are flattened finely, then
    resampled at uniform arc length (loop length / 1000) keeping entity boundaries and
    polyline vertices as corners. UNVERIFIED on a real corpus. Cluster = file.
    """
    def gen():
        import ezdxf
        from ezdxf import path as ezpath
        for fp in sorted(glob.glob(str(data_root() / pattern), recursive=True)):
            try:
                doc = ezdxf.readfile(fp)
            except Exception:
                continue
            msp = doc.modelspace()
            pieces = []      # (dense points, corner flags)
            for k, e in enumerate(msp):
                if e.dxftype() not in ("LINE", "ARC", "LWPOLYLINE", "POLYLINE", "SPLINE", "CIRCLE", "ELLIPSE"):
                    continue
                try:
                    p = ezpath.make_path(e)
                    bb = np.array([(v.x, v.y) for v in p.control_vertices()])
                    tol = 1e-5 * max(float(np.ptp(bb, axis=0).max()) if len(bb) > 1 else 1.0, 1e-9)
                    dense = np.array([(v.x, v.y) for v in p.flattening(tol)])
                except Exception:
                    continue
                if len(dense) >= 2:
                    pieces.append((e.dxftype(), dense))
            # chain LINE/ARC pieces greedily
            singles = [d for t, d in pieces if t not in ("LINE", "ARC")]
            loose = [d for t, d in pieces if t in ("LINE", "ARC")]
            chains = []
            while loose:
                ch = [loose.pop(0)]
                grown = True
                while grown:
                    grown = False
                    end = ch[-1][-1]
                    for i, d in enumerate(loose):
                        tol = 1e-6 * max(diameter(np.vstack(ch)), 1e-9)
                        if np.hypot(*(d[0] - end)) <= tol:
                            ch.append(loose.pop(i)); grown = True; break
                        if np.hypot(*(d[-1] - end)) <= tol:
                            ch.append(loose.pop(i)[::-1]); grown = True; break
                chains.append(ch)
            curves = [[d] for d in singles] + chains
            for ci, parts in enumerate(curves):
                prims = []
                for d in parts:
                    s = arclength(d)
                    if s[-1] <= 0:
                        continue
                    prims.append(((lambda t, d=d, s=s: np.column_stack([
                        np.interp(t * s[-1], s, d[:, 0]), np.interp(t * s[-1], s, d[:, 1])])), float(s[-1])))
                if not prims:
                    continue
                pts = _uniform_with_corners(prims)
                yield Curve("", f"{_rel(fp)}#{ci}", _rel(fp), pts, is_closed(pts))
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
    Dataset("ne_rivers", "maps", load_geojson("natural_earth/ne_10m_rivers_lake_centerlines.geojson"),
            "Natural Earth 1:10m river centerlines", "natural_earth/ne_10m_rivers_lake_centerlines.geojson"),
    Dataset("ne_roads", "maps", load_geojson("natural_earth/ne_10m_roads.geojson"),
            "Natural Earth 1:10m roads", "natural_earth/ne_10m_roads.geojson"),
    Dataset("gpx", "GPS", load_gpx(), "GPX track segments and routes (metres)", "gpx/**/*.gpx"),
    Dataset("quickdraw_raw", "pen strokes", load_quickdraw_raw(), "Google Quick Draw raw strokes",
            "quickdraw/*.ndjson"),
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


def _content_key(pts: np.ndarray) -> str:
    d = max(diameter(pts), 1e-300)
    q = np.round((pts - pts.min(0)) / d * 1e9).astype(np.int64)
    return hashlib.sha1(q.tobytes()).hexdigest()


def raw_files(dataset: str) -> list[str]:
    return sorted(glob.glob(str(data_root() / DATASETS[dataset].files), recursive=True))


def load_split(dataset: str, split: str, cap: int | None = None, stats: dict | None = None) -> list[Curve]:
    """Eligible curves of one split, as a blind simple random sample of at most ``cap``.

    The sample is the ``cap`` curves with the smallest sha256(seed, dataset, curve_id),
    which is a uniform random sample that can be drawn while streaming large files.
    Exact duplicate curves are removed globally (first curve id kept) before splitting,
    so a duplicate can never sit in both dev and test.
    """
    ds = DATASETS[dataset]
    cap = C.CAP[split] if cap is None else cap
    st = {"raw_curves": 0, "too_short": 0, "duplicates": 0, "eligible": {"dev": 0, "test": 0},
          "clusters": {"dev": set(), "test": set()}}
    seen: dict[str, str] = {}
    heap: list = []   # max-heap on rank via negative
    for cv in ds.loader():
        st["raw_curves"] += 1
        pts = dedupe(cv.pts)
        if len(pts) < C.MIN_VERTICES or not np.isfinite(pts).all() or arclength(pts)[-1] <= 0:
            st["too_short"] += 1
            continue
        key = _content_key(pts)
        if key in seen:
            st["duplicates"] += 1
            continue
        seen[key] = cv.curve_id
        sp = split_of(dataset, cv.cluster)
        st["eligible"][sp] += 1
        st["clusters"][sp].add(cv.cluster)
        if sp != split or cap <= 0:
            continue
        cv = Curve(dataset, cv.curve_id, cv.cluster, pts, cv.closed or is_closed(pts), cv.meta)
        rank = _unit_hash(C.SEED, dataset, cv.curve_id)
        if len(heap) < cap:
            heapq.heappush(heap, (-rank, cv.curve_id, cv))
        elif -heap[0][0] > rank:
            heapq.heapreplace(heap, (-rank, cv.curve_id, cv))
    st["clusters"] = {k: len(v) for k, v in st["clusters"].items()}
    if stats is not None:
        stats.update(st)
    return sorted((c for _, _, c in heap), key=lambda c: c.curve_id)
