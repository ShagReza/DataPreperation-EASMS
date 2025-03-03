# -*- coding: utf-8 -*-
"""
Created on Fri Feb 28 12:26:45 2025

@author: shagh
"""

import pandas as pd

def generate_ml_labels(df):
    """
    Assigns BINARY_LABEL based on ENRICHMENT and PVALUE criteria.

    Args:
        df (pd.DataFrame): The input DataFrame.

    Returns:
        pd.DataFrame: Updated DataFrame with the BINARY_LABEL column.
    """

    # Ensure necessary columns exist
    required_columns = {"ENRICHMENT", "PVALUE"}
    if not required_columns.issubset(df.columns):
        raise ValueError(f"Missing required columns: {required_columns - set(df.columns)}")

    # Assign BINARY_LABEL based on conditions
    df["BINARY_LABEL"] = df.apply(lambda row: "Y" if row["ENRICHMENT"] > 10 and row["PVALUE"] < 0.05 else "N", axis=1)

    return df
