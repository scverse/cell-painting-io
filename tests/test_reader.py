from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from cell_painting_io import read_profiles


@pytest.fixture
def profiles(tmp_path: Path) -> list[Path]:
    frame = pd.DataFrame({"Metadata_Well": ["A01", "A02"], "Cells_AreaShape_Area": [10.0, 20.0]})
    paths = []
    for name, rows in (("plate1", frame), ("plate2", frame.iloc[:0]), ("plate3", frame)):
        path = tmp_path / name / "profile.csv"
        path.parent.mkdir()
        rows.to_csv(path, index=False)
        paths.append(path)
    return paths


def test_a_file_with_no_rows_does_not_erase_the_features(profiles: list[Path]) -> None:
    adata = read_profiles(profiles)

    assert adata.shape == (4, 1)
    assert adata.var_names.tolist() == ["Cells_AreaShape_Area"]
    assert np.array_equal(adata.X.ravel(), [10.0, 20.0, 10.0, 20.0])


def test_path_columns_stay_aligned_when_a_file_is_dropped(profiles: list[Path]) -> None:
    adata = read_profiles(profiles, path_columns={"Metadata_Plate": 1})

    assert adata.obs["Plate"].tolist() == ["plate1", "plate1", "plate3", "plate3"]


def test_all_files_empty_yields_no_observations(profiles: list[Path]) -> None:
    # a header-only file carries no dtype information, so there is no way to tell which columns are features
    adata = read_profiles([profiles[1]])

    assert adata.n_obs == 0
