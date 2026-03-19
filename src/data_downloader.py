import csv
import requests
import os
import sys
import time
import logging
import hashlib

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.dirname(BASE_DIR)
LOG_DIR = os.path.join(PROJECT_DIR, "logs")
TSV_DIR = os.path.join(PROJECT_DIR, "data", "tsv")
RAW_DIR = os.path.join(PROJECT_DIR, "data", "raw")

for directory in [LOG_DIR, TSV_DIR, RAW_DIR]:
    os.makedirs(directory, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s | %(levelname)-7s | %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    handlers=[
        logging.FileHandler(os.path.join(LOG_DIR, "downloader.log")),
        logging.StreamHandler(sys.stdout)
    ]
)

class BaseDownloader:
    """Shared network and I/O utilities for all downloaders."""

    def request_data(self, url, headers=None):
        """Executes a GET request and returns a tuple (json_data, None) or (None, error_message)."""
        try:
            response = requests.get(url, headers=headers, timeout=(15, 60))
            response.raise_for_status()
            return response.json(), None
        except requests.exceptions.RequestException as e:
            return None, str(e)
        except ValueError as e:
            return None, f"JSON Decode Error: {str(e)}"

    def save_tsv(self, data, filename):
        """Saves a list of lists to a TSV file."""
        with open(filename, 'w', newline='') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerows(data)

    def calculate_md5(self, file_path):
        """Reads a file in 1MB chunks and returns its hexadecimal MD5 hash."""
        hash_md5 = hashlib.md5()
        with open(file_path, "rb") as f:
            for chunk in iter(lambda: f.read(1048576), b""):
                hash_md5.update(chunk)
        return hash_md5.hexdigest()

    def _format_size(self, size_in_bytes):
        """Converts raw bytes into a human-readable string (KB, MB, GB)."""
        for unit in ['B', 'KB', 'MB', 'GB']:
            if size_in_bytes < 1024.0:
                return f"{size_in_bytes:.2f} {unit}"
            size_in_bytes /= 1024.0
        return f"{size_in_bytes:.2f} TB"

    def _check_existing_file(self, dest_path, expected_size, expected_md5):
        """Validates a local file via size (short-circuit) and MD5. Returns a boolean."""
        if not os.path.exists(dest_path):
            return False

        local_size = os.path.getsize(dest_path)

        if expected_size > 0 and local_size != expected_size:
            os.remove(dest_path)
            return False

        if expected_md5:
            if self.calculate_md5(dest_path) == expected_md5:
                return True
            os.remove(dest_path)
            return False

        return True

    def _perform_download_and_verify(self, url, dest_path, expected_size, expected_md5):
        """Streams the file download and applies strict size and MD5 validation with UX feedback."""
        downloaded_bytes = 0
        file_name = os.path.basename(dest_path)
        start_time = time.time()

        with requests.get(url, stream=True, timeout=(15, 300)) as r:
            r.raise_for_status()
            with open(dest_path, 'wb') as f:
                for chunk in r.iter_content(chunk_size=1048576):
                    if chunk:
                        f.write(chunk)
                        downloaded_bytes += len(chunk)
                        
                        elapsed_time = time.time() - start_time
                        speed = (downloaded_bytes / 1048576) / elapsed_time if elapsed_time > 0 else 0
                        
                        if expected_size > 0:
                            sys.stdout.write(f"\rDownloading {file_name}: {self._format_size(downloaded_bytes)} / {self._format_size(expected_size)} @ {speed:.1f} MB/s")
                        else:
                            sys.stdout.write(f"\rDownloading {file_name}: {self._format_size(downloaded_bytes)} / Unknown @ {speed:.1f} MB/s")
                        sys.stdout.flush()
        
        sys.stdout.write("\n")

        final_size = os.path.getsize(dest_path)
        
        if expected_size > 0 and final_size != expected_size:
            raise ValueError(f"Size mismatch: {final_size} != {expected_size}")

        if expected_md5:
            sys.stdout.write(f"Verifying MD5 checksum for {file_name}... ")
            sys.stdout.flush()
            if self.calculate_md5(dest_path) != expected_md5:
                sys.stdout.write("FAILED!\n")
                raise ValueError("MD5 mismatch.")
            sys.stdout.write("OK!\n")

        return True

    def download_file(self, url, dest_path, expected_size, expected_md5, retries=3):
        """Manages the download cycle, validation, and retries upon failure."""
        exp_size = int(expected_size) if str(expected_size).isdigit() else 0
        file_name = os.path.basename(dest_path)

        for attempt in range(retries):
            try:
                if self._check_existing_file(dest_path, exp_size, expected_md5):
                    sys.stdout.write(f"{file_name}: Already downloaded and verified. Skipping.\n")
                    return True

                return self._perform_download_and_verify(url, dest_path, exp_size, expected_md5)

            except Exception as e:
                if os.path.exists(dest_path):
                    os.remove(dest_path)

                if attempt == retries - 1:
                    logging.error(f"Failed to download {file_name}: {e}")
                    return False

                time.sleep(2)

    def process_download_queue(self, tsv_path, output_dir):
        """Reads TSV metadata and queues the download for each validated link."""
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        with open(tsv_path, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                links = row['Download_Links'].split(';')
                sizes = row['File_Sizes'].split(';')
                md5s = row['MD5_Checksums'].split(';')
                is_paired = row['Paired_End'].upper() == 'TRUE'

                for i, link in enumerate(links):
                    suffix = f"_{i+1}" if (is_paired and len(links) > 1) else ""
                    file_name = f"{row['Accession']}{suffix}.fastq.gz"
                    full_dest = os.path.join(output_dir, file_name)

                    exp_size = sizes[i] if i < len(sizes) else "0"
                    exp_md5 = md5s[i] if i < len(md5s) else ""

                    self.download_file(link, full_dest, exp_size, exp_md5)

class EncodeDownloader(BaseDownloader):
    """Implements metadata extraction via the ENCODE Project REST API."""

    def _extract_fastq_info(self, data, condition):
        """Filters the ENCODE JSON to retain only released and paired FASTQ files."""
        extracted = []
        base_url = "https://www.encodeproject.org"
        for file in data.get('files', []):
            is_fastq = (file.get('file_type') == 'fastq' or file.get('file_format') == 'fastq')
            if is_fastq and file.get('status') == 'released':
                rep = file.get('replicate') or {}
                
                pair_val = str(file.get('paired_end'))
                acc = file.get('accession', 'N/A')
                if pair_val in ['1', '2']:
                    acc = f"{acc}_{pair_val}"
                
                is_paired = pair_val in ['1', '2']
                
                extracted.append([
                    condition,
                    acc,
                    is_paired,
                    str(rep.get('biological_replicate_number', 'N/A')),
                    base_url + file.get('href', ''),
                    str(file.get('file_size', '0')),
                    file.get('md5sum', '')
                ])
        return extracted

    def fetch(self, experiment_id):
        """Consumes the ENCODE endpoint and returns the standardized metadata."""
        url = f"https://www.encodeproject.org/experiments/{experiment_id}/?format=json"
        raw_data, error = self.request_data(url, headers={'accept': 'application/json'})

        if error or not isinstance(raw_data, dict):
            return [], error or "Invalid data received from API"

        assay = raw_data.get('assay_term_name', 'Unknown')
        target = raw_data.get('target', {}).get('label', '')
        condition = f"{assay}_{target}".strip('_')

        parsed_data = self._extract_fastq_info(raw_data, condition)
        if not parsed_data:
            return [], "No released FASTQs found"
        return parsed_data, None

class EnaDownloader(BaseDownloader):
    """Implements metadata extraction via the ENA portal, acting as a backend for GEO and SRA."""

    def _extract_fastq_info(self, data, condition):
        """Evaluates the priority of FTP links (FASTQ vs SRA) in the ENA API response."""
        extracted = []
        if not isinstance(data, list):
            return []
        for file_data in data:
            layout = file_data.get('library_layout', 'N/A')

            ftp_raw = file_data.get('fastq_ftp', '')
            bytes_raw = file_data.get('fastq_bytes', '')
            md5_raw = file_data.get('fastq_md5', '')

            if not ftp_raw:
                ftp_raw = file_data.get('sra_ftp', '')
                bytes_raw = file_data.get('sra_bytes', '')
                md5_raw = file_data.get('sra_md5', '')

            ftp_links = [
                f"http://{link}" if not link.startswith('http') else link
                for link in ftp_raw.split(';') if link
            ]

            extracted.append([
                condition,
                file_data.get('run_accession', 'N/A'),
                True if layout == "PAIRED" else False,
                file_data.get('sample_accession', 'N/A'),
                ";".join(ftp_links),
                bytes_raw,
                md5_raw
            ])
        return extracted

    def fetch(self, accession_id):
        """Consumes the ENA filereport endpoint to request read metadata."""
        url = (
            f"https://www.ebi.ac.uk/ena/portal/api/filereport"
            f"?accession={accession_id}&result=read_run"
            f"&fields=run_accession,sample_accession,library_layout,"
            f"fastq_ftp,fastq_bytes,fastq_md5,sra_ftp,sra_bytes,sra_md5"
            f"&format=json"
        )
        raw_data, error = self.request_data(url)
        if error:
            return [], error
        if not raw_data:
            return [], "Empty response"
        return self._extract_fastq_info(raw_data, "ENA_Data"), None

class GeoDownloader(EnaDownloader):
    """Download interface for GEO experiments (GSE/PRJ prefixes)."""
    pass

class SraDownloader(EnaDownloader):
    """Download interface for SRA runs (SRR/ERR prefixes)."""
    pass

class TcgaDownloader(BaseDownloader):
    """Placeholder interface for future integration with the TCGA GDC-Client."""

    def fetch(self, experiment_id):
        raise NotImplementedError("TCGA requires GDC-Client.")

if __name__ == "__main__":
    descriptions = {
        "ENCSR000EUN": "p53 ChIP-seq",
        "PRJNA148505": "p53 RNA-seq",
        "ENCSR000CPK": "p53 RNA-seq",
        "ENCSR042AWH": "ATAC-seq",
        "ENCSR000EUM": "H3K4me3",
        "ENCSR000EUL": "H3K27ac",
        "ENCSR000EUK": "H3K27me3"
    }

    target_ids = sys.argv[1:] if len(sys.argv) > 1 else list(descriptions.keys())

    downloaders_map = {
        "ENC": EncodeDownloader, "GSE": GeoDownloader, "PRJ": GeoDownloader,
        "SRR": SraDownloader, "ERR": SraDownloader, "TCG": TcgaDownloader
    }

    logging.info("Starting processing")

    for eid in target_ids:
        prefix = eid[:3].upper()
        downloader_class = downloaders_map.get(prefix)

        if not downloader_class:
            logging.warning(f"[SKIP] {eid}: Unknown prefix.")
            continue

        try:
            downloader = downloader_class()
            data, err = downloader.fetch(eid)

            if err:
                logging.error(f"[{eid}] Error: {err}")
                continue

            if data:
                tsv_path = os.path.join(TSV_DIR, f"{eid}_metadata.tsv")
                header = [['Condition', 'Accession', 'Paired_End', 'Replicate',
                           'Download_Links', 'File_Sizes', 'MD5_Checksums']]
                downloader.save_tsv(header + data, tsv_path)

                logging.info(f"[{eid}] Downloading {len(data)} file(s)...")
                downloader.process_download_queue(tsv_path, RAW_DIR)
                logging.info(f"[{eid}] Done")
            else:
                logging.warning(f"[{eid}] No valid data found.")

        except NotImplementedError as nie:
            logging.warning(f"[{eid}] Warning: {nie}")
        except Exception as e:
            logging.critical(f"[{eid}] Critical failure: {e}")

    logging.info("Processing complete")