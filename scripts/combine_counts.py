import pandas as pd
from glob import glob

# Recursively search for count files
count_files = sorted(glob("results/counts/*/*_counts.txt"))
combined = None

for file in count_files:
    sample = file.split("/")[-1].split("_counts.txt")[0]
    
    # Load only Geneid and count column, skipping all header lines starting with '#'
    df = pd.read_csv(
        file,
        sep="\t",
        comment="#",   # skip lines starting with #
        header=0       # use the actual header from featureCounts
    )

    df = df[["Geneid", df.columns[-1]]]  # Gene ID and the last column (counts)
    df.columns = ["Geneid", sample]
    df = df.set_index("Geneid")

    if combined is None:
        combined = df
    else:
        combined = combined.join(df, how="outer")

combined = combined.fillna(0).astype(int)
combined.to_csv("results/counts/all_counts.tsv", sep="\t")
