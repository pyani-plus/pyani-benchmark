#!/bin/env python
# -*- coding: utf-8 -*-
"""utils.py
Defines helper functions and variables.
"""

from matplotlib.colors import LinearSegmentedColormap

# Colormap for pairwise identity heatmap
my_cmap = LinearSegmentedColormap.from_list(
    "my_gradient",
    (
        # Edit this gradient at https://eltos.github.io/gradient/#80:2369BC-95:FFFFFF-100:A8373B
        (0.000, (0.137, 0.412, 0.737)),
        (0.750, (1.000, 1.000, 1.000)),
        (1.000, (0.659, 0.216, 0.231)),
    ),
)
