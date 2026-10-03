#!/usr/bin/env python3
"""
erode_stl.py - shrink an STL in X/Y (keeping Z unchanged) so it fits into a
hole of the same shape, with an optional bevel that tapers the part toward
one face.

Typical use: lettering (an extruded outline) that has to drop into a pocket
cut with the exact same outline. Printed at the same size it won't fit, so
this script insets every wall by a fixed clearance.

How it works
------------
A horizontal cross-section of the model is offset inward with a true 2D
erosion (manifold3d / Clipper2 polygon offset) and re-extruded. Z is never
scaled, so overall height and the top/bottom planes are preserved exactly.

  * No bevel: the model is treated as a straight extrusion (lettering
    usually is), so one cross-section is taken at mid-height and re-extruded
    over the full height. This is exact.
  * --bevel: the inset grows linearly from --offset at the "full" face to
    --offset + --bevel at the beveled face (bottom by default). The default
    'smooth' method cuts true sloped walls by subtracting the Minkowski sum
    of the outside region with a cone. The 'stepped' method instead stacks
    thin layers (--step, default 0.04 mm - below the A1's minimum layer
    height); it is faster but produces much larger files.
  * --reslice: re-slice the original at every --step layer instead of
    assuming a straight extrusion. Use this for models that aren't simple
    prisms (e.g. lettering that already has chamfers or a base plate).
    Bevels are then always built stepped.

Examples
--------
  # 0.15 mm clearance on every side
  python erode_stl.py letters.stl letters_fit.stl --offset 0.15

  # 0.1 mm clearance at the top face, 0.1 + 0.3 = 0.4 mm at the bottom face
  python erode_stl.py letters.stl letters_fit.stl --offset 0.1 --bevel 0.3

  # Taper only the bottom 1 mm of the part
  python erode_stl.py letters.stl out.stl --offset 0.1 --bevel 0.3 --bevel-height 1

Requirements:  pip install trimesh manifold3d numpy
"""

from __future__ import annotations

import argparse
import math
import sys

import numpy as np
import trimesh
from manifold3d import CrossSection, JoinType, Manifold, Mesh, OpType

JOIN_TYPES = {"round": JoinType.Round, "miter": JoinType.Miter, "square": JoinType.Square}


def trimesh_to_manifold(mesh: trimesh.Trimesh) -> Manifold:
    m = Manifold(
        Mesh(
            vert_properties=np.ascontiguousarray(mesh.vertices, dtype=np.float32),
            tri_verts=np.ascontiguousarray(mesh.faces, dtype=np.uint32),
        )
    )
    if m.is_empty():
        raise SystemExit(
            f"Input mesh is not a closed, watertight solid ({m.status()}). "
            "Repair it first (e.g. Bambu Studio/Orca 'Fix model', or Meshmixer)."
        )
    return m


def manifold_to_trimesh(m: Manifold) -> trimesh.Trimesh:
    mesh = m.to_mesh()
    verts = np.asarray(mesh.vert_properties)[:, :3]
    faces = np.asarray(mesh.tri_verts)
    return trimesh.Trimesh(vertices=verts, faces=faces, process=False)


def build(solid: Manifold, args) -> Manifold:
    (_, _, zmin), (_, _, zmax) = np.reshape(solid.bounding_box(), (2, 3))
    height = zmax - zmin
    if height <= 0:
        raise SystemExit("Model has zero height.")

    bevel = args.bevel
    bevel_height = height if args.bevel_height is None else min(args.bevel_height, height)
    if bevel > 0 and bevel_height <= 0:
        raise SystemExit("--bevel-height must be > 0.")

    def offset_at(z: float) -> float:
        """Total inset at height z."""
        if bevel <= 0:
            return args.offset
        # distance from the beveled face
        d = (z - zmin) if args.bevel_side == "bottom" else (zmax - z)
        if d >= bevel_height:
            return args.offset
        return args.offset + bevel * (1.0 - d / bevel_height)

    def slice_at(z: float) -> CrossSection:
        # Slicing a wall made of triangles leaves extra collinear points whose
        # position depends on z (where the plane crosses each diagonal). Drop
        # them and snap to a 0.1 micron grid so layers cut from the same
        # straight wall come out identical; otherwise near-coincident vertices
        # from neighbouring layers become non-manifold edges in a float32 STL.
        polys = [np.round(p, 4) for p in solid.slice(z).simplify(1e-4).to_polygons()]
        return CrossSection(polys) if polys else CrossSection()

    base_section = None
    if not args.reslice:
        base_section = slice_at(zmin + height / 2.0)
        if base_section.is_empty():
            raise SystemExit("Mid-height cross-section is empty; try --reslice.")

    def section_at(z: float) -> CrossSection:
        return base_section if base_section is not None else slice_at(z)

    def erode(cs: CrossSection, distance: float) -> CrossSection:
        if distance > 0:
            cs = cs.offset(-distance, JOIN_TYPES[args.join], args.miter_limit, args.segments)
        # Drop near-duplicate points, which would otherwise collapse into
        # degenerate triangles once the STL is written in float32.
        return cs.simplify(1e-4)

    if base_section is not None and (bevel <= 0 or args.bevel_method == "smooth"):
        # Straight extrusion: one slab with the plain offset is exact.
        cs = erode(base_section, args.offset)
        if cs.is_empty():
            raise SystemExit("Nothing left after erosion - the offset is larger than the features.")
        prism = Manifold.extrude(cs, height).translate((0.0, 0.0, zmin))
        if bevel <= 0:
            return prism
        return smooth_bevel(prism, cs, zmin, zmax, bevel, bevel_height, args)

    # Layered build: used for --reslice and --bevel-method stepped.
    if bevel <= 0 and not args.reslice:
        bounds = [zmin, zmax]
    else:
        n = max(1, math.ceil(height / args.step))
        bounds = list(np.linspace(zmin, zmax, n + 1))

    slabs = []
    for z0, z1 in zip(bounds[:-1], bounds[1:]):
        zc = 0.5 * (z0 + z1)
        cs = erode(section_at(zc), offset_at(zc))
        if not cs.is_empty():
            slabs.append(Manifold.extrude(cs, z1 - z0).translate((0.0, 0.0, z0)))

    if not slabs:
        raise SystemExit("Nothing left after erosion - the offset is larger than the features.")
    # Collapse sliver edges left where adjacent layers nearly line up.
    return Manifold.batch_boolean(slabs, OpType.Add).simplify(1e-4)


def smooth_bevel(prism: Manifold, cs: CrossSection, zmin: float, zmax: float,
                 bevel: float, bevel_height: float, args) -> Manifold:
    """Cut a linear bevel into a straight prism.

    Every point within `bevel * (1 - d / bevel_height)` (horizontally) of the
    prism's outside, where d is the distance from the beveled face, is
    removed. That region is the Minkowski sum of the outside, taken as a slab
    at the beveled face, with a cone of radius `bevel` tapering to a point at
    height `bevel_height`. Subtracting it gives clean conical side walls that
    stay parallel to the original outline (inside corners get rounded).
    """
    x0, y0, x1, y1 = cs.bounds()
    margin = bevel + 1.0
    frame = CrossSection.square((x1 - x0 + 2 * margin, y1 - y0 + 2 * margin)).translate(
        (x0 - margin, y0 - margin)
    )
    # 1 mm thick slab of "outside" just beyond the beveled face.
    outside = Manifold.extrude(frame - cs, 1.0)
    # The Minkowski sum adds positions, so the cone stays at the origin:
    # radius `bevel` at z=0, tapering to a point `bevel_height` into the part.
    cone = Manifold.cylinder(bevel_height, bevel, 0.0, args.cone_segments)
    if args.bevel_side == "bottom":
        outside = outside.translate((0.0, 0.0, zmin - 1.0))
    else:
        outside = outside.translate((0.0, 0.0, zmax))
        cone = cone.mirror((0.0, 0.0, 1.0))
    print("computing smooth bevel (Minkowski sum; can take a minute for long text)...")
    cutter = outside.minkowski_sum(cone)
    result = prism - cutter
    if result.is_empty():
        raise SystemExit("Nothing left after beveling - the bevel is larger than the features.")
    return result


def report(label: str, mesh: trimesh.Trimesh):
    lo, hi = mesh.bounds
    size = hi - lo
    print(
        f"{label}: {len(mesh.faces)} triangles, "
        f"size X {size[0]:.3f}  Y {size[1]:.3f}  Z {size[2]:.3f} mm, "
        f"watertight={mesh.is_watertight}"
    )


def main(argv=None):
    p = argparse.ArgumentParser(
        description="Erode an STL in X/Y by a fixed clearance (Z unchanged), "
        "with an optional linear bevel toward one face.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("input", help="input STL")
    p.add_argument("output", help="output STL")
    p.add_argument("-o", "--offset", type=float, required=True,
                   help="inset applied to every side wall, in mm (per side)")
    p.add_argument("-b", "--bevel", type=float, default=0.0,
                   help="extra inset at the beveled face, in mm; tapers linearly "
                        "to 0 at the other face (or at --bevel-height)")
    p.add_argument("--bevel-side", choices=["bottom", "top"], default="bottom",
                   help="which face gets the extra inset (the other face keeps just --offset)")
    p.add_argument("--bevel-height", type=float, default=None,
                   help="height in mm over which the bevel tapers (default: full model height)")
    p.add_argument("--bevel-method", choices=["smooth", "stepped"], default="smooth",
                   help="'smooth' cuts true sloped walls (Minkowski sum, slower); "
                        "'stepped' stacks thin layers of increasing inset (fast, larger file)")
    p.add_argument("--cone-segments", type=int, default=32,
                   help="facets around the bevel cone for --bevel-method smooth "
                        "(more = rounder inside corners, slower)")
    p.add_argument("--step", type=float, default=0.04,
                   help="layer thickness in mm for --bevel-method stepped and --reslice")
    p.add_argument("--reslice", action="store_true",
                   help="slice the original at every layer instead of assuming a straight extrusion")
    p.add_argument("--join", choices=sorted(JOIN_TYPES), default="round",
                   help="corner style for the offset; 'round' is a true erosion "
                        "(inside corners rounded, outside corners stay sharp)")
    p.add_argument("--miter-limit", type=float, default=2.0,
                   help="miter limit (multiples of the offset) when --join miter")
    p.add_argument("--segments", type=int, default=64,
                   help="segments per full circle for rounded corners")
    args = p.parse_args(argv)

    if args.offset < 0 or args.bevel < 0:
        p.error("--offset and --bevel must be >= 0")
    if args.step <= 0:
        p.error("--step must be > 0")

    mesh = trimesh.load(args.input, force="mesh")
    if not isinstance(mesh, trimesh.Trimesh) or len(mesh.faces) == 0:
        p.error(f"could not read a mesh from {args.input}")
    mesh.merge_vertices()
    report("input ", mesh)

    out = manifold_to_trimesh(build(trimesh_to_manifold(mesh), args))
    out.export(args.output)
    report("output", out)
    print(f"wrote {args.output}")


if __name__ == "__main__":
    sys.exit(main())
