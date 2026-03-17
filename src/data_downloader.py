import csv
import requests
import os

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
            with requests.get(url, stream=True) as r:
                r.raise_for_status()
                with open(dest_path, 'wb') as f:
                    for chunk in r.iter_content(chunk_size=8192):
                        f.write(chunk)
            return True
        except requests.exceptions.RequestException:
            return False

class EncodeDownloader(BaseDownloader):
    def _extract_fastq_info(self, data, condition):
        extracted = []
        for file in data.get('files', []):
            if file.get('file_type') == 'fastq':
                rep = file.get('replicate') or {}
                extracted.append([
                    condition,
                    file.get('accession', 'N/A'),
                    file.get('paired_end', 'N/A'),
                    rep.get('biological_replicate_number', 'N/A')
                ])
        return extracted

    def fetch(self, experiment_id):
        url = f"https://www.encodeproject.org/experiments/{experiment_id}/?format=json"
        
        try:
            raw_data = self.request_data(url, headers={'accept': 'application/json'})
            if not raw_data:
                return []

            parsed_data = []
            parsed_data.extend(self._extract_fastq_info(raw_data, 'KD'))

            controles = set()
            for file in raw_data.get('files', []):
                rep = file.get('replicate') or {}
                exp = rep.get('experiment') or {}
                controles.update(exp.get('possible_controls', []))

            for ctrl in controles:
                ctrl_id = ctrl.strip('/')
                ctrl_url = f"https://www.encodeproject.org/{ctrl_id}/?format=json"
                ctrl_data = self.request_data(ctrl_url, headers={'accept': 'application/json'})
                if ctrl_data:
                    parsed_data.extend(self._extract_fastq_info(ctrl_data, 'Control'))

            return parsed_data
            
        except Exception:
            return []

class GeoDownloader(BaseDownloader):
    def fetch(self, geo_id):
        pass

class SraDownloader(BaseDownloader):
    def fetch(self, sra_id):
        pass

class TcgaDownloader(BaseDownloader):
    def fetch(self, tcga_id):
        pass

if __name__ == "__main__":
    downloader = EncodeDownloader()
    dados = downloader.fetch("ENCSR000EUN")
    if dados:
        cabecalho = [['Condicao', 'Accession', 'Paired_End', 'Replicata']]
        downloader.save_tsv(cabecalho + dados, "ENCSR000EUN_metadata.tsv")