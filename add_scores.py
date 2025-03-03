import os
import pandas as pd
import numpy as np
import scipy.stats as stats

def compute_and_add_scores(file_paths):
    """
    Computes TARGET_VALUE, ENRICHMENT, SELECTIVE_ENRICHMENT, and PVALUE for a list of CSV files.
    
    Args:
        file_paths (list): List of file paths for separated CSV files.
    
    Returns:
        None (modifies files in-place).
    """
    if not file_paths:
        print("No files provided for score computation.")
        return

    # Load all CSV files into a dictionary
    dataframes = {f: pd.read_csv(f) for f in file_paths}

    # Ensure necessary columns exist
    required_columns = {"COMPOUND_ID", "POS_INT_REP1", "POS_INT_REP2", "POS_INT_REP3"}
    for filename, df in dataframes.items():
        if not required_columns.issubset(df.columns):
            raise ValueError(f"Missing required columns in {filename}: {required_columns - set(df.columns)}")

    # Step 1: Compute TARGET_VALUE for each file
    for df in dataframes.values():
        df["TARGET_VALUE"] = df[["POS_INT_REP1", "POS_INT_REP2", "POS_INT_REP3"]].mean(axis=1, skipna=True)

    # Step 2: Compute SELECTIVE_VALUE and NTC_VALUE across all files
    merged_df = pd.concat(dataframes.values(), ignore_index=True)

    # Compute max and min values for each COMPOUND_ID
    max_values = merged_df.groupby("COMPOUND_ID")["TARGET_VALUE"].max()
    min_values = merged_df.groupby("COMPOUND_ID")["TARGET_VALUE"].min()

    # Apply values to each file
    for df in dataframes.values():
        df["SELECTIVE_VALUE"] = df["COMPOUND_ID"].map(max_values)
        df["NTC_VALUE"] = df["COMPOUND_ID"].map(min_values)

        # Calculate enrichment scores
        df["ENRICHMENT"] = df["TARGET_VALUE"] / df["NTC_VALUE"]
        df["SELECTIVE_ENRICHMENT"] = df["TARGET_VALUE"] / df["SELECTIVE_VALUE"]

        # Compute PVALUE with variance check
        def calculate_p_value(row):
            """Computes p-value using t-test, avoiding precision loss by checking variance."""
            # Extract replicas for the current compound
            protein_interest = row[["POS_INT_REP1", "POS_INT_REP2", "POS_INT_REP3"]].dropna().tolist()

            # Extract replicas for other compounds with the same COMPOUND_ID
            protein_other = merged_df[merged_df["COMPOUND_ID"] == row["COMPOUND_ID"]]
            protein_other_values = protein_other[["POS_INT_REP1", "POS_INT_REP2", "POS_INT_REP3"]].stack().dropna().tolist()

            # Check if both lists have sufficient values
            if len(protein_interest) == 0 or len(protein_other_values) < 3:
                return None  # Not enough data for statistical testing

            # Convert to numpy arrays
            protein_interest = np.array(protein_interest)
            protein_other_values = np.array(protein_other_values)

            # Check variance to avoid precision loss
            if np.std(protein_interest) < 1e-8 or np.std(protein_other_values) < 1e-8:
                return 1.0  # Assign high p-value if there's no meaningful difference

            # Perform t-test
            _, p_value = stats.ttest_ind(protein_interest, protein_other_values, equal_var=False)
            
            return p_value

        df["PVALUE"] = df.apply(calculate_p_value, axis=1)

    # Step 3: Save the updated files
    for file_path, df in dataframes.items():
        df.to_csv(file_path, index=False)
        print(f"Updated scores saved: {file_path}")

    print("\nAll separated files have been updated with computed scores.")
