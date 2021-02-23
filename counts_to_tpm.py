import os
import sys
import argparse
import subprocess
import shutil
import re
import glob
import shlex
import operator
from collections import defaultdict
import time
import ast
import pandas as pd
import numpy as np
from pathlib import Path
import textwrap
from sklearn import preprocessing
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import ListedColormap, LinearSegmentedColormap, Normalize
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1 import make_axes_locatable
from sklearn.linear_model import LinearRegression
from scipy import stats
from Bio import Phylo
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform
from scipy.cluster import hierarchy
import dendropy
from natsort import natsorted

from io import StringIO
from skbio import read
from skbio.tree import TreeNode

from skbio.diversity import alpha_diversity
from skbio.diversity import beta_diversity
from skbio.stats.distance import mantel
from skbio.stats.ordination import pcoa
from skbio.stats.distance import anosim
from bioinfokit.analys import norm, get_data
# load sugarcane RNA-seq expression dataset (Published in Bedre et al., 2019)
df = pd.read_csv("all_OM252_counts.csv", sep="\t", header=0, index_col=0)

# make gene column as index column
#df = df.set_index('ACC')


# now, normalize raw counts using TPM method
# gene length must be in bp
nm = norm()
nm.tpm(df=df, gl='genome_length')
# get TPM normalized dataframe
tpm_df = nm.tpm_norm

tpm_df = tpm_df.T

tpm_df.to_csv("test_tpm.tsv", sep="\t")
