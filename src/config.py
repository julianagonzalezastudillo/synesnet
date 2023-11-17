import os.path
from pathlib import Path

# File paths or directories
DATA_DIR = Path("/Users/juliana.gonzalez/ownCloud/graph_analysis/")

NET_DIR = DATA_DIR / "net_metrics/new/new"  # net_path
NET_DIR.mkdir(parents=True, exist_ok=True)

RAND_DIR = DATA_DIR / "rand_mat_strength(equal)"
RAND_DIR.mkdir(parents=True, exist_ok=True)

PLOT_DIR = Path(os.getcwd(), *["plots", "glb", "new"])
PLOT_DIR.mkdir(parents=True, exist_ok=True)


# Files
NODE_FILE = os.path.join(DATA_DIR, "BN_Atlas_246_LUT_reoriented.txt")
INFO_FILE = os.path.join(
    DATA_DIR, "resultsROI_Subject001_Condition001.mat"
)  # results_file


# Application settings
P_VAL = 0.05
CORR_TYPE = "_thr"
N_SUB = len(
    os.listdir(os.path.join(DATA_DIR, "symetrical_corr_mat"))
)  # number of subjects
