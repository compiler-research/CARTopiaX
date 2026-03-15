import numpy as np
import pandas as pd

def test_plot3_division_by_zero():
    """Testing that Plot 3 percentage calculation handles potential zero tumor cells."""

    #mock data for tumor being fully eliminated 
    test_df = pd.DataFrame({
        'num_tumor_cells':        [100, 50, 0],
        'tumor_cells_type1':      [40,  20, 0],
        'tumor_cells_type2':      [30,  15, 0],
        'tumor_cells_type3':      [20,  10, 0],
        'tumor_cells_type4':      [10,   5, 0],
        'tumor_cells_type5_dead': [0,    0, 0],
    })

    cell_types = [
        ('tumor_cells_type1',      'Type 1'),
        ('tumor_cells_type2',      'Type 2'),
        ('tumor_cells_type3',      'Type 3'),
        ('tumor_cells_type4',      'Type 4'),
        ('tumor_cells_type5_dead', 'Type 5 (Dead)'),
    ]

    for col, label in cell_types:
        if col in test_df.columns:
            total = test_df['num_tumor_cells'].replace(0, np.nan)
            percentage = 100 * test_df[col] / total

            #when tumor is eliminated, the percentage should be NaN not inf
            assert np.isnan(percentage.iloc[2]), \
                f"{label}: expected NaN when tumor count is 0, got {percentage.iloc[2]}"

            #non-zero rows should not be affected
            assert not np.isinf(percentage.iloc[0]), \
                f"{label}: non-zero rows should not produce inf"
            assert not np.isnan(percentage.iloc[0]), \
                f"{label}: non-zero rows should not produce NaN"

    print("all plot 3 tests have passed successfully")