# STL tools

## erode_stl.py

Shrinks an STL in X/Y by a fixed clearance (Z unchanged) so a part fits into a
hole of the exact same shape, e.g. lettering inlaid into a matching pocket.
Optionally adds a bevel, so the part keeps its size at one face and gets
linearly smaller toward the other. That makes it easy to start in the hole and
press in snugly.

```sh
pip install -r requirements.txt

# 0.15 mm clearance on every side wall
python erode_stl.py letters.stl letters_fit.stl --offset 0.15

# 0.1 mm clearance at the top, growing linearly to 0.4 mm at the bottom
python erode_stl.py letters.stl letters_fit.stl --offset 0.1 --bevel 0.3
```

| Option | Meaning |
| --- | --- |
| `-o/--offset MM` | Inset applied to every side wall (per side). |
| `-b/--bevel MM` | Extra inset at the beveled face, tapering linearly to 0. |
| `--bevel-side bottom\|top` | Which face gets the bevel (default `bottom`). |
| `--bevel-height MM` | Taper only over this much height (default: full height). |
| `--bevel-method smooth\|stepped` | `smooth` cuts true sloped walls; `stepped` stacks thin layers. |
| `--reslice` | Re-slice every layer, for models that aren't straight extrusions. |
| `--join round\|miter\|square` | Corner treatment. `round` is a true erosion: outside corners stay sharp, inside corners are rounded. |

Run `python erode_stl.py --help` for all options.

### Tips for the Bambu A1

* Start with about 0.1–0.15 mm per side for PLA and print a short test word
  first. Elephant's foot on the first layers often needs a little more
  clearance at the bottom, and a small `--bevel` covers that too.
* The bevel is measured relative to the model's orientation in the STL, so
  "bottom" means the face that sits on the build plate as exported.
* The input must be a closed (watertight) mesh. Separate overlapping letter
  bodies are fine; they are merged.
