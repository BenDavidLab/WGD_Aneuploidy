import os

ROOT_FOLDER="/specific/elkon/kavyap/aneuploidy"
CODE_FOLDER=os.path.join(ROOT_FOLDER,"code")
CACHE_FOLDER=os.path.join(ROOT_FOLDER,"cache")
CACHE_FOLDER_TCGA=os.path.join(ROOT_FOLDER,"cache")
CACHE_FOLDER_ROW=os.path.join(ROOT_FOLDER,"cache_rowperm")
OUTPUT_FOLDER=os.path.join(ROOT_FOLDER,"output")
DATASETS_FOLDER=os.path.join(ROOT_FOLDER,"datasets")

INTRA=0
INTER=1

CHR_INTERACTION_NAMES={INTRA:"intra",
                   INTER: "inter"}
