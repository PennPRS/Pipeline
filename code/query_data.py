"""
Query data utilities for PennPRS: resolve and download GWAS Catalog and FinnGen summary statistics.
"""

import os
import logging
from typing import Optional

import wget

logger = logging.getLogger(__name__)

# Default paths (override with environment variables if needed)
DEFAULT_HARMONISED_FILE = os.environ.get(
    "PENNPRS_HARMONISED_FILE",
    "/home/ubuntu/pennprs/frontend/public/harmonised.txt",
)
DEFAULT_GWAS_BASE_URL = "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
DEFAULT_GWAS_DATA_DIR = os.environ.get("PENNPRS_GWAS_DATA_DIR", "/home/ubuntu/data2/gwas_data/")
DEFAULT_FINNGEN_DATA_DIR = os.environ.get(
    "PENNPRS_FINNGEN_DATA_DIR", "/home/ubuntu/data2/finn12_data/"
)
FINNGEN_FILE_PREFIX = "EUR_finngen_R12_"
FINNGEN_FILE_SUFFIX = ".txt"


def get_query_path(
    trait_id: str,
    harmonised_file: Optional[str] = None,
) -> Optional[str]:
    """
    Resolve a GWAS Catalog trait ID to its summary statistics URL.

    Looks up the trait ID in the harmonised index file (EBI GWAS Catalog paths).
    Returns the full URL if found, None otherwise.

    Args:
        trait_id: GWAS Catalog study/trait identifier (e.g. GCST...).
        harmonised_file: Path to harmonised.txt. Defaults to PENNPRS_HARMONISED_FILE or built-in default.

    Returns:
        Full URL string, or None if trait_id is not in the index.
    """
    path = harmonised_file or DEFAULT_HARMONISED_FILE
    if not os.path.isfile(path):
        logger.warning("Harmonised file not found: %s", path)
        return None
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            # Format: prefix/TRACT_ID/suffix
            parts = line.split("/")
            if len(parts) >= 2 and trait_id == parts[1]:
                return DEFAULT_GWAS_BASE_URL + line
    return None


def get_finngen_trait_path(
    trait_id: str,
    data_dir: Optional[str] = None,
    prefix: Optional[str] = None,
    suffix: Optional[str] = None,
) -> str:
    """
    Build the local file path for a FinnGen trait (R12 EUR summary statistics).

    Does not check whether the file exists; use os.path.exists() if needed.

    Args:
        trait_id: FinnGen phenocode.
        data_dir: Base directory for FinnGen data. Defaults to PENNPRS_FINNGEN_DATA_DIR or built-in default.
        prefix: Filename prefix (default EUR_finngen_R12_).
        suffix: Filename suffix (default .txt).

    Returns:
        Absolute path to the expected summary statistics file.
    """
    base = data_dir or DEFAULT_FINNGEN_DATA_DIR
    pre = prefix if prefix is not None else FINNGEN_FILE_PREFIX
    suf = suffix if suffix is not None else FINNGEN_FILE_SUFFIX
    return os.path.join(base, f"{pre}{trait_id}{suf}")


def finngen_trait_exists(
    trait_id: str,
    data_dir: Optional[str] = None,
) -> bool:
    """
    Check whether a FinnGen trait file exists on disk.

    Args:
        trait_id: FinnGen phenocode.
        data_dir: Base directory for FinnGen data (optional).

    Returns:
        True if the file exists.
    """
    path = get_finngen_trait_path(trait_id, data_dir=data_dir)
    return os.path.isfile(path)


def download_gwas_catalog_file(
    trait_id: str,
    trait_url: str,
    data_dir: Optional[str] = None,
    output_dir: Optional[str] = None,
    log_file_path: Optional[str] = None,
) -> Optional[str]:
    """
    Download a GWAS Catalog summary statistics file by URL and save as trait_id + extension.

    If the file already exists at data_dir/{trait_id}{extension}, the download is skipped.
    On failure, a message is appended to query_data.log (in output_dir or at log_file_path).

    Args:
        trait_id: GWAS Catalog trait/study ID (used for the local filename).
        trait_url: Full URL to the summary statistics file.
        data_dir: Directory to download into. Defaults to PENNPRS_GWAS_DATA_DIR or built-in default.
        output_dir: Directory for log file if log_file_path not set.
        log_file_path: Optional path for query_data.log. If None, uses output_dir/query_data.log when output_dir is set.

    Returns:
        Path to the downloaded (or existing) file, or None on failure.
    """
    data_dir = data_dir or DEFAULT_GWAS_DATA_DIR
    os.makedirs(data_dir, exist_ok=True)

    file_name = os.path.basename(trait_url)
    ext = os.path.splitext(file_name)[1]
    new_filename = f"{trait_id}{ext}"
    new_file_path = os.path.join(data_dir, new_filename)

    if os.path.exists(new_file_path):
        logger.info("File %s already exists in %s, skipping download.", new_filename, data_dir)
        return new_file_path

    try:
        cwd = os.getcwd()
        downloaded = None
        try:
            os.chdir(os.path.dirname(data_dir) or "/")
            downloaded = wget.download(trait_url, out=data_dir)
        finally:
            os.chdir(cwd)
            if downloaded:
                temp_path = downloaded + ".tmp"
                if os.path.exists(temp_path):
                    try:
                        os.remove(temp_path)
                    except OSError:
                        pass

        if downloaded and os.path.exists(downloaded):
            if os.path.abspath(downloaded) != os.path.abspath(new_file_path):
                os.rename(downloaded, new_file_path)
            else:
                new_file_path = downloaded
            return new_file_path
    except Exception as e:
        logger.exception("Failed to download GWAS file: %s", e)
        log_path = log_file_path
        if log_path is None and output_dir:
            log_path = os.path.join(output_dir, "query_data.log")
        if log_path:
            os.makedirs(os.path.dirname(log_path), exist_ok=True)
            with open(log_path, "a") as log_file:
                log_file.write(
                    f"Failed to download the file from GWAS Catalog: {e}. "
                    f"Trait ID: {trait_id}. Please ensure the GWAS Catalog ID is in the "
                    "queryable data list on https://pennprs.org/data.\n"
                )
        return None

    return None


def resolve_gwas_trait_path(trait_id: str, harmonised_file: Optional[str] = None) -> Optional[str]:
    """
    Resolve a GWAS Catalog trait ID to its download URL. Convenience wrapper around get_query_path.

    Args:
        trait_id: GWAS Catalog study/trait identifier.
        harmonised_file: Optional path to harmonised index file.

    Returns:
        URL string or None if not found.
    """
    return get_query_path(trait_id, harmonised_file=harmonised_file)


def resolve_finngen_trait_path(
    trait_id: str,
    data_dir: Optional[str] = None,
    check_exists: bool = True,
) -> Optional[str]:
    """
    Resolve a FinnGen phenocode to its local file path.

    Args:
        trait_id: FinnGen phenocode.
        data_dir: Base directory for FinnGen data (optional).
        check_exists: If True, return None when the file does not exist.

    Returns:
        Absolute file path, or None if check_exists=True and file is missing.
    """
    path = get_finngen_trait_path(trait_id, data_dir=data_dir)
    if check_exists and not os.path.isfile(path):
        return None
    return path


# --- CLI for standalone use ---

def _main():
    import argparse
    parser = argparse.ArgumentParser(
        description="Resolve or download PennPRS queryable trait data (GWAS Catalog, FinnGen)."
    )
    parser.add_argument(
        "source",
        choices=["gwas", "finngen"],
        help="Data source: gwas (GWAS Catalog) or finngen (FinnGen R12).",
    )
    parser.add_argument(
        "trait_id",
        help="Trait identifier: GWAS Catalog ID (e.g. GCST...) or FinnGen phenocode.",
    )
    parser.add_argument(
        "--resolve-only",
        action="store_true",
        help="Only print the resolved path/URL; do not download (GWAS only).",
    )
    parser.add_argument(
        "--download",
        action="store_true",
        help="For GWAS: download the file into the default data dir. Ignored for FinnGen (local only).",
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Directory for logs (and optional output). Default: current directory.",
    )
    parser.add_argument(
        "--harmonised-file",
        default=None,
        help="Path to harmonised.txt (GWAS only). Overrides PENNPRS_HARMONISED_FILE.",
    )
    parser.add_argument(
        "--gwas-data-dir",
        default=None,
        help="Directory for GWAS downloads. Overrides PENNPRS_GWAS_DATA_DIR.",
    )
    parser.add_argument(
        "--finngen-data-dir",
        default=None,
        help="Base directory for FinnGen files. Overrides PENNPRS_FINNGEN_DATA_DIR.",
    )
    args = parser.parse_args()

    if args.source == "gwas":
        url = get_query_path(
            args.trait_id,
            harmonised_file=args.harmonised_file or DEFAULT_HARMONISED_FILE,
        )
        if url is None:
            print("NOT_FOUND")
            return 1
        print(url)
        if args.download and not args.resolve_only:
            out_dir = args.output_dir or os.getcwd()
            path = download_gwas_catalog_file(
                args.trait_id,
                url,
                data_dir=args.gwas_data_dir or DEFAULT_GWAS_DATA_DIR,
                output_dir=out_dir,
            )
            if path:
                print(f"DOWNLOADED: {path}")
                return 0
            return 1
        return 0

    if args.source == "finngen":
        path = get_finngen_trait_path(
            args.trait_id,
            data_dir=args.finngen_data_dir or DEFAULT_FINNGEN_DATA_DIR,
        )
        print(path)
        if not os.path.isfile(path):
            print("NOT_FOUND", file=__import__("sys").stderr)
            return 1
        return 0

    return 0


if __name__ == "__main__":
    raise SystemExit(_main())
