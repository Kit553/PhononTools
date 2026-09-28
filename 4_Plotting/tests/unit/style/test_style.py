"""Geometry, palette, and export-policy tests."""

from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import matplotlib as mpl
from matplotlib.colors import is_color_like
from pypdf import PdfReader

from catplot import (
    DEFAULT_EDGE_MM,
    FINAL_DPI,
    apply_style,
    load_theme,
    mm_to_inches,
    new_panel,
    save_panel,
)


class StyleTests(unittest.TestCase):
    def test_default_geometry_is_exact_square(self) -> None:
        fig, _, _ = new_panel()
        width, height = fig.get_size_inches()
        self.assertLess(abs(width - height), 1e-12)
        self.assertLess(abs(width - mm_to_inches(DEFAULT_EDGE_MM)), 1e-12)

    def test_base_export_policy(self) -> None:
        apply_style()
        # Matplotlib normalizes the style value ``standard`` to None.
        self.assertIsNone(mpl.rcParams["savefig.bbox"])
        self.assertEqual(mpl.rcParams["savefig.dpi"], FINAL_DPI)
        self.assertEqual(mpl.rcParams["pdf.fonttype"], 42)
        self.assertFalse(mpl.rcParams["text.usetex"])

    def test_theme_colours_are_valid_and_indigo_led(self) -> None:
        theme = load_theme("default")
        self.assertEqual(theme["series"][0], "#332288")
        self.assertTrue(all(is_color_like(color) for color in theme["series"]))
        self.assertTrue(
            all(is_color_like(color) for color in theme["semantic"].values())
        )

    def test_exported_pdf_page_is_square(self) -> None:
        with TemporaryDirectory() as directory:
            fig, ax, _ = new_panel()
            ax.plot([0, 1], [0, 1])
            output = save_panel(fig, Path(directory) / "panel.pdf")
            page = PdfReader(output).pages[0]
            width_pt = float(page.mediabox.width)
            height_pt = float(page.mediabox.height)
            self.assertLess(abs(width_pt - height_pt), 0.01)
            expected_pt = mm_to_inches(DEFAULT_EDGE_MM) * 72
            self.assertLess(abs(width_pt - expected_pt), 0.02)

    def test_non_square_export_is_rejected(self) -> None:
        with TemporaryDirectory() as directory:
            fig, _, _ = new_panel()
            fig.set_size_inches(4, 3)
            with self.assertRaisesRegex(ValueError, "non-square"):
                save_panel(fig, Path(directory) / "bad.pdf")


if __name__ == "__main__":
    unittest.main()

