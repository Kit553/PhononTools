# Style decisions and invariants

## Hard invariants

1. A single exported panel has a square PDF page.
2. Export never uses a tight bounding box.
3. Physical page dimensions are explicit and testable.
4. PDF text remains embedded text, not paths or bitmap glyphs.
5. Raster layers use at least 600 dpi in a final PDF.
6. Plot recipes may select a named style profile but may not override arbitrary
   Matplotlib rcParams. This prevents plot files from becoming hidden style
   sheets.

## Canvas square versus axes square

The current invariant applies to the **PDF artboard**, not the inner data box.
Axis labels, tick labels, legends, and colour bars occupy part of that artboard,
so forcing the inner data box to be square as well can waste space or cause
clipping. Individual plot types may request a square data box when its geometry
requires it, for example equal-scale x/y maps.

## Multipanel consequence

A linked grid will be sized from square cells, but the complete exported block
need not itself be square. A vertical stack of three square cells is naturally
taller than it is wide. This keeps single-panel geometry uniform without
distorting multipanel compositions.

## Typography

The default font request is Arial, followed by Helvetica, Nimbus Sans,
Liberation Sans, and DejaVu Sans. This gives a familiar scientific sans-serif
appearance on Windows and Linux without making a proprietary font a runtime
requirement. MathText uses DejaVu Sans and external LaTeX is disabled so clean
reproduction does not depend on a local TeX installation.

## Resolution

"600 dpi PDF" only has meaning for content that must be rasterized. Lines,
markers, patches, and text remain vector objects. Dense heatmaps and very large
scatter collections may be selectively rasterized at 600 dpi to control PDF
size without degrading labels or axes.

## Affinity boundary

Affinity may place, align, label, and combine complete quantitative panels with
VESTA/VMD renderings or schematics. Data geometry, ticks, axes labels, line
styles, and quantitative annotations should be regenerated from the recipe
rather than edited inside Affinity.

