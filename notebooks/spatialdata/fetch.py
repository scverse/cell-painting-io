"""Download one well of one plate of a Cell Painting Gallery source, and nothing else.

The gallery is tens of terabytes, so the notebooks in this directory work from a slice of it: the well-level
profile of a whole plate, but the images and CellProfiler output of a single well.

    python fetch.py cpg0008-pki/broad 2021_04_07_Batch1 BR00122970 A01

Prints the profile it found, which is the `profile` argument `cell_painting_io.read_plate` then needs.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd
import s3fs

BUCKET = "cellpainting-gallery"
SKIP = ("Cytoplasm.csv", "Experiment.csv", "Image.csv", "cp.is.done")
PROFILE_ORDER = (
    "_normalized_feature_select_negcon_batch.csv.gz",
    "_normalized_feature_select_batch.csv.gz",
    "_normalized.csv.gz",
    ".csv.gz",
    ".parquet",
)

fs = s3fs.S3FileSystem(anon=True)


def download(keys: list[str], destinations: list[Path]) -> None:
    def one(pair: tuple[str, Path]) -> None:
        key, destination = pair
        if destination.exists():
            return
        destination.parent.mkdir(parents=True, exist_ok=True)
        fs.get(f"{BUCKET}/{key}", str(destination))

    with ThreadPoolExecutor(16) as pool:
        list(pool.map(one, zip(keys, destinations, strict=True)))


def find_profile(source: str, batch: str, plate: str) -> str | None:
    directory = f"{BUCKET}/{source}/workspace/profiles/{batch}/{plate}"
    if not fs.exists(directory):
        return None
    names = {Path(path).name for path in fs.ls(directory)}
    return next((f"{plate}{suffix}" for suffix in PROFILE_ORDER if f"{plate}{suffix}" in names), None)


def fetch(source: str, batch: str, plate: str, well: str, root: Path) -> str | None:
    load_data = root / "workspace/load_data_csv" / batch / plate / "load_data.csv"
    download([f"{source}/workspace/load_data_csv/{batch}/{plate}/load_data.csv"], [load_data])

    profile = find_profile(source, batch, plate)
    if profile is not None:
        download(
            [f"{source}/workspace/profiles/{batch}/{plate}/{profile}"],
            [root / "workspace/profiles" / batch / plate / profile],
        )

    analysis = f"{source}/workspace/analysis/{batch}/{plate}/analysis"
    keys, destinations = [], []
    for directory in fs.ls(f"{BUCKET}/{analysis}"):
        if not Path(directory).name.startswith(f"{plate}-{well}-"):
            continue
        for key in fs.find(directory):
            name = key.split(f"{analysis}/")[1]
            if not name.endswith(SKIP):
                keys.append(key.removeprefix(f"{BUCKET}/"))
                destinations.append(root / "workspace/analysis" / batch / plate / "analysis" / name)
    download(keys, destinations)

    frame = pd.read_csv(load_data)
    frame = frame[frame["Metadata_Well"].astype(str) == well]
    channels = [c for c in frame.columns if c.startswith("URL_Orig")] or [
        c for c in frame.columns if c.startswith("FileName_Orig")
    ]
    keys, destinations = [], []
    for _, row in frame.iterrows():
        for column in channels:
            if column.startswith("URL_"):
                location = str(row[column])
            else:
                location = f"{row['PathName_Orig' + column.removeprefix('FileName_Orig')]}/{row[column]}"
            tail = location[location.index(f"/{batch}/images/") + 1 :]
            keys.append(f"{source}/images/{tail}")
            destinations.append(root / "images" / tail)
    download(keys, destinations)
    return profile


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("source", help="accession and source, as in cpg0008-pki/broad")
    parser.add_argument("batch")
    parser.add_argument("plate")
    parser.add_argument("well")
    parser.add_argument("--root", type=Path, default=Path("~/data"), help="where the sources are kept")
    args = parser.parse_args()

    root = args.root.expanduser() / args.source
    profile = fetch(args.source, args.batch, args.plate, args.well, root)
    print(f"{root}\nprofile: {profile}")


if __name__ == "__main__":
    main()
