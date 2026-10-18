import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from nmcode import DataRepository, arguments, finish_figure, result_directory, set_style, write_table
from nmcode.plotting import PALETTE, heatmap
from nmcode.statistics import benjamini_hochberg, donor_composition, kruskal_by_feature, mann_whitney_by_feature, partial_spearman_table, spearman_table, welch_differential


def setup(description, name):
    args = arguments(description)
    set_style()
    return args, DataRepository(args.data_root), result_directory(args.output_root, name)
