# -*- coding: utf-8 -*-
"""
Created on Fri Feb 28 11:32:57 2025

@author: shay
"""

import os
import pandas as pd

def filter_anomalous_data(df, sep_file_name):
    """
    Filters duplicate rows and processes SMILES with different ENRICHMENT values:
    - Removes fully duplicate rows.
    - If all rows for a SMILES have ENRICHMENT < 10, keeps only the row with the smallest ENRICHMENT.
    - If all rows for a SMILES have ENRICHMENT > 10, keeps all rows but logs them.

    Args:
        df (pd.DataFrame): The input DataFrame.
        sep_file_name (str): The name of the separated CSV file being processed.

    Returns:
        pd.DataFrame: Cleaned DataFrame with anomalies handled.
    """

    # Ensure the necessary columns exist
    required_columns = {"SMILES", "ENRICHMENT"}
    if not required_columns.issubset(df.columns):
        raise ValueError(f"Missing required columns: {required_columns - set(df.columns)}")

    # Step 1: Remove fully duplicate rows
    df_cleaned = df.drop_duplicates()

    # Step 2: Identify SMILES that have multiple ENRICHMENT values
    enrichment_groups = df_cleaned.groupby("SMILES")["ENRICHMENT"].nunique()
    conflicting_smiles = enrichment_groups[enrichment_groups > 1].index.tolist()

    # Prepare log dataframe
    log_df = pd.DataFrame(columns=["SMILES", "ENRICHMENT"])

    # Step 3: Process conflicting SMILES
    rows_to_keep = []
    for smiles in conflicting_smiles:
        subset = df_cleaned[df_cleaned["SMILES"] == smiles]  # Get all rows for this SMILES

        # Check if all ENRICHMENT values are < 10
        if subset["ENRICHMENT"].max() < 10:
            # Keep only the row with the lowest ENRICHMENT
            best_row = subset.loc[subset["ENRICHMENT"].idxmin()]
            rows_to_keep.append(best_row)
        elif subset["ENRICHMENT"].min() > 10:
            # Keep all rows but log them
            log_df = pd.concat([log_df, subset[["SMILES", "ENRICHMENT"]]])

            # Keep all rows (do nothing)
            rows_to_keep.extend(subset.to_dict(orient="records"))
        else:
            # Keep all rows if mix of >10 and <10 (no action needed)
            rows_to_keep.extend(subset.to_dict(orient="records"))

    # Convert filtered rows to DataFrame
    filtered_df = pd.DataFrame(rows_to_keep)

    # Step 4: Merge back with non-conflicting SMILES
    final_df = pd.concat([df_cleaned[~df_cleaned["SMILES"].isin(conflicting_smiles)], filtered_df])

    # Step 5: Save log file (if applicable)
    if not log_df.empty:
        log_file_path = os.path.join(os.getcwd(), f"DuplicatesLog_{sep_file_name}.csv")
        log_df.to_csv(log_file_path, index=False)
        print(f"Logged {len(log_df)} SMILES entries with ENRICHMENT > 10 in {log_file_path}")

    return final_df
