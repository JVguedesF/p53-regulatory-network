import csv
import requests
import os
import sys

class BaseDownloader:
    def request_data(self, url, headers=None):
        try:
            response = requests.get(url, headers=headers)
            response.raise_for_status()
            return response.json()
        except requests.exceptions.RequestException:
            return {}

    def save_tsv(self, data, filename):
        with open(filename, 'w', newline='') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerows(data)

    def download_file(self, url, dest_path):
        try:
            response_head = requests.head(url, allow_redirects=True)
            remote_size = int(response_head.headers.get('content-length', 0))

            if os.path.exists(dest_path):
                local_size = os.path.getsize(dest_path)
                if local_size == remote_size and remote_size > 0:
                    return True
                else:
                    print(f"      - File corrupted or incomplete ({local_size}/{remote_size}). Retrying...")

            with requests.get(url, stream=True) as r:
                r.raise_for_status()
                with open(dest_path, 'wb') as f:
                    for chunk in r.iter_content(chunk_size=1048576):
                        f.write(chunk)
            return True
        except Exception:
            return False

    def process_download_queue(self, tsv_path, output_dir):
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        with open(tsv_path, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                accession = row['Accession']
                links = row['Download_Links'].split(';')
                is_paired = row['Paired_End'].upper() == 'TRUE'
                
                print(f"  > Processing: {accession}")
                
                for i, link in enumerate(links):
                    suffix = f"_{i+1}" if is_paired else ""
                    file_name = f"{accession}{suffix}.fastq.gz"
                    full_dest = os.path.join(output_dir, file_name)
                    
                    print(f"    - Downloading to: {file_name}")
                    self.download_file(link, full_dest)

class EncodeDownloader(BaseDownloader):
    def _extract_fastq_info(self, data, condition):
        extracted = []
        base_url = "https://www.encodeproject.org"
        for file in data.get('files', []):
            if file.get('file_type') == 'fastq':
                rep = file.get('replicate') or {}
                download_link = base_url + file.get('href', '')
                extracted.append([
                    condition,
                    file.get('accession', 'N/A'),
                    file.get('paired_end', 'N/A'),
                    str(rep.get('biological_replicate_number', 'N/A')),
                    download_link
                ])
        return extracted

    def fetch(self, experiment_id):
        url = f"https://www.encodeproject.org/experiments/{experiment_id}/?format=json"
        try:
            raw_data = self.request_data(url, headers={'accept': 'application/json'})
            if not raw_data: return []

            parsed_data = []
            parsed_data.extend(self._extract_fastq_info(raw_data, 'KD'))

            controls = set()
            for file in raw_data.get('files', []):
                rep = file.get('replicate') or {}
                exp = rep.get('experiment') or {}
                controls.update(exp.get('possible_controls', []))

            for ctrl in controls:
                ctrl_id = ctrl.strip('/')
                ctrl_url = f"https://www.encodeproject.org/{ctrl_id}/?format=json"
                ctrl_data = self.request_data(ctrl_url, headers={'accept': 'application/json'})
                if ctrl_data:
                    parsed_data.extend(self._extract_fastq_info(ctrl_data, 'Control'))

            return parsed_data
        except Exception:
            return []

class EnaDownloader(BaseDownloader):
    def _extract_fastq_info(self, data, condition):
        extracted = []
        if not isinstance(data, list): return []
        for idx, file_data in enumerate(data, start=1):
            layout = file_data.get('library_layout', 'N/A')
            paired = True if layout == "PAIRED" else False
            ftp_raw = file_data.get('fastq_ftp', '')
            ftp_links = [f"http://{link}" for link in ftp_raw.split(';') if link]
            extracted.append([
                condition,
                file_data.get('run_accession', 'N/A'),
                paired,
                str(idx),
                ";".join(ftp_links)
            ])
        return extracted

    def fetch(self, accession_id, condition="KD"):
        url = f"https://www.ebi.ac.uk/ena/portal/api/filereport?accession={accession_id}&result=read_run&fields=run_accession,library_layout,fastq_ftp&format=json"
        try:
            raw_data = self.request_data(url)
            return self._extract_fastq_info(raw_data, condition) if raw_data else []
        except Exception:
            return []
        
class GeoDownloader(EnaDownloader):
    pass

class SraDownloader(EnaDownloader):
    pass

class TcgaDownloader(BaseDownloader):
    pass

if __name__ == "__main__":
    descriptions = {
        "ENCSR000EUN": "p53 ChIP-seq",
        "GSE33256":    "p53 RNA-seq (Activation)",
        "ENCSR913HOH": "p53 RNA-seq (Knockdown)",
        "ENCSR453NUV": "ATAC-seq (Chromatin)",
        "ENCSR000EUM": "H3K4me3 (Promoters)",
        "ENCSR000EUL": "H3K27ac (Enhancers)",
        "ENCSR000EUK": "H3K27me3 (Repression)",
        "TCG-BRCA":    "TCGA Breast Cancer",
        "TCG-COAD":    "TCGA Colon Cancer"
    }

    target_ids = sys.argv[1:] if len(sys.argv) > 1 else list(descriptions.keys())
    
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_dir = os.path.dirname(script_dir)
    tsv_dir = os.path.join(project_dir, "data", "tsv")
    raw_dir = os.path.join(project_dir, "data", "raw")

    os.makedirs(tsv_dir, exist_ok=True)
    os.makedirs(raw_dir, exist_ok=True)

    downloaders_map = {
        "ENC": EncodeDownloader, "GSE": GeoDownloader, "PRJ": GeoDownloader,
        "SRR": SraDownloader, "ERR": SraDownloader, "TCG": TcgaDownloader
    }

    print("\n" + "="*100)
    print(f"{'Experiment ID':<15} | {'Type':<5} | {'Description':<25} | {'Status'}")
    print("-" * 100)

    for experiment_id in target_ids:
        prefix = experiment_id[:3].upper()
        desc = descriptions.get(experiment_id, "Unknown Dataset")
        downloader_class = downloaders_map.get(prefix)

        if not downloader_class:
            print(f"{experiment_id:<15} | {prefix:<5} | {desc:<25} | [SKIP] Invalid Prefix")
            continue

        downloader = downloader_class()
        print(f"{experiment_id:<15} | {prefix:<5} | {desc:<25} | [BUSY] Fetching...", end="\r")
        
        data = downloader.fetch(experiment_id)
        if data:
            tsv_path = os.path.join(tsv_dir, f"{experiment_id}_metadata.tsv")
            header = [['Condition', 'Accession', 'Paired_End', 'Replicate', 'Download_Links']]
            downloader.save_tsv(header + data, tsv_path)
            
            print(f"{experiment_id:<15} | {prefix:<5} | {desc:<25} | [BUSY] Downloading...", end="\r")
            downloader.process_download_queue(tsv_path, raw_dir)
            print(f"{experiment_id:<15} | {prefix:<5} | {desc:<25} | [DONE] Success")
        else:
            print(f"{experiment_id:<15} | {prefix:<5} | {desc:<25} | [FAIL] No Data")

    print("="*100 + "\n")