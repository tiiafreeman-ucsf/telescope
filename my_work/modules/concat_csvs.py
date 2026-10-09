import pandas as pd

def combine_csvs(file1, file2, match_column, output_file):
    df1 = pd.read_csv(file1)
    df2 = pd.read_csv(file2)

    # Match rows using the specified column
    combined_df = df1.merge(
        df2,
        on=match_column,
        how="left"
    )

    combined_df.to_csv(output_file, index=False)

    return combined_df


combined_df = combine_csvs(
    file1="my_work/HERV_datasets/HervD Atlas Database(HerD Elements).csv",
    file2="my_work/HERV_datasets/HervD Atlas Database(HervD).csv",
    match_column="HERV ID",
    output_file="combined.csv"
)

print(combined_df.head())