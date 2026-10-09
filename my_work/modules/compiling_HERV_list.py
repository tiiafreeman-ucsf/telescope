##For Telescope, first steps are to collect:
    #RNA seq FASTQs, sample metadata, Human reerence genome, HERV annotations
    #hg38 or GRCh28 

## should optimize aggregation if using again!! takes forever

import pandas as pd
from pathlib import Path

HERVd_inputs = r"my_work/HERV_datasets/HERVd human files.csv"
HervD_atlas_inputs = r"my_work/HERV_datasets/HervD Atlas Database(HervD).csv"
gEVE_inputs = r"my_work/HERV_datasets/gEVE_database_human.csv"

output_dir = Path(r"my_work/output")
output_dir.mkdir(parents=True, exist_ok=True)

def load_hervd(HERVd_inputs):
    df = pd.read_csv(HERVd_inputs)
    df.columns = df.columns.str.strip()
    df = df.loc[:, ~df.columns.str.startswith("Unnamed")].copy()

    columns = {
        "name": "HERVd_ID",
        "chrom": "Chromosome",
        "chromStart": "Start",
        "chromEnd": "End",
        "strand": "Strand",
        "thickStart": "HERVd_thickStart",
        "thickEnd": "HERVd_thickEnd",
        "blockCount": "HERVd_blockCount",
        "blockSizes": "HERVd_blockSizes",
        "blockStarts": "HERVd_blockStarts"
    }

    df = df[list(columns.keys())].rename(columns=columns).copy()
    df["Chromosome"] = df["Chromosome"].astype("string").str.strip().str.removeprefix("chr")

    parts = df["HERVd_ID"].str.split(":", n=2, expand=True)

    df["HERVd_num"] = parts[0]
    df["HERVd_repeat"] = parts[2]

    df["source_database"] = "HERVd"

    return df

def load_hervd_atlas(HervD_atlas_inputs):

    df = pd.read_csv(HervD_atlas_inputs)
    df.columns = df.columns.str.strip()
    df = df.loc[:, ~df.columns.str.startswith("Unnamed")].copy()

    columns = {
        "HERV ID": "HERVd_atlas_ID",
        "HERV Name": "HERVd_atlas_Name",
        "Strand": "Strand",
        "Position": "Position",
        "Length": "HERVd_atlas_Length",
        "Region Type": "HERVd_atlas_Region_Type",
        "Type": "HERVd_atlas_Type",
        "Group": "HERVd_atlas_Group"
    }

    df = df[list(columns.keys())].rename(columns=columns).copy()

    parts = df["Position"].str.split(":", n=1, expand=True)
    df["Chromosome"] = parts[0].astype("string").str.strip().str.removeprefix("chr")

    positions = parts[1].str.split("-", n=1, expand=True)
    df["Start"] = pd.to_numeric(positions[0], errors="raise")
    df["End"] = pd.to_numeric(positions[1], errors="raise")

    df = df.drop(columns=["Position"]).copy()

    df["source_database"] = "HERVd_atlas"

    return df


def load_geve(gEVE_inputs):
    
    df = pd.read_csv(gEVE_inputs)
    df.columns = df.columns.str.strip()
    df = df.loc[:, ~df.columns.str.startswith("Unnamed")].copy()

    columns = {
        "ID": "gEVE_ID",
        "chr": "Chromosome",
        "start": "Start",
        "end": "End",
        "strand": "Strand",
        "AA_length": "gEVE_AA_length",
        "method": "gEVE_method",
        "N_cnt": "gEVE_N_cnt",
        "MetID": "gEVE_MetID",
        "Met_AA_length": "gEVE_Met_AA_length",
        "HMM_profile": "gEVE_HMM_profile",
        "Viral_BLAST": "gEVE_Viral_BLAST",
        "NR_BLAST": "gEVE_NR_BLAST",
        "EVE_BLAST": "gEVE_EVE_BLAST",
        "RetroTector": "gEVE_RetroTector",
        "Repbase": "gEVE_Repbase"
    }

    df = df[list(columns.keys())].rename(columns=columns).copy()
    df["source_database"] = "gEVE"

    return df

def combine_herv_datasets(hervd_df, hervd_atlas_df, geve_df):
    combined_df = pd.concat(
        [geve_df, hervd_df, hervd_atlas_df],
        ignore_index=True
    )

    match_columns = ["Chromosome", "Start", "End", "Strand"]

    def concatenate_values(series):
        values = series.dropna().astype(str)
        values = values[values.str.strip() != ""]
        return ";".join(values) if not values.empty else pd.NA
    
    other_columns = [
        col for col in combined_df.columns if col not in match_columns
    ]

    combined_df = combined_df.groupby(match_columns, as_index=False).agg({col:concatenate_values for col in other_columns})

    return combined_df

hervd_df = load_hervd(HERVd_inputs)
hervd_atlas_df = load_hervd_atlas(HervD_atlas_inputs)
geve_df = load_geve(gEVE_inputs)

hervd_df.to_csv(output_dir / "HERVd_processed.csv", index=False)
hervd_atlas_df.to_csv(output_dir / "HERVd_atlas_processed.csv", index=False)
geve_df.to_csv(output_dir / "gEVE_processed.csv", index=False)

print(f"HERVd: {len(hervd_df):,} records exported")
print(f"HERVd Atlas: {len(hervd_atlas_df):,} records exported")
print(f"gEVE: {len(geve_df):,} records exported")
print(f"\nFiles saved to: {output_dir}")

combined_df = combine_herv_datasets(hervd_df, hervd_atlas_df, geve_df)

combined_df.to_csv(output_dir / "combined_herv_annotations.csv", index=False)

print(f"\nCombined file saved to: {output_dir}")