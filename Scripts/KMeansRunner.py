from random import seed
from os.path import join
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
from typing import List, Dict, Tuple
import pandas as pd
import collections
from tqdm import tqdm

FOLDER_PATH = "C:/Users/Patrick/Documents/Github/P6/Dataset"
FILE_NAME = "features.tsv"
OUTPUT_NAME = "kmeans_out_"
OUTPUT_TYPE = ".tsv"

seed(2)

class TsvLoader:
    @staticmethod
    def load_tsv(file) -> List:
        result = []
        indexes = []

        with open(file, 'r') as f:
            lines = f.readlines()
            for line_idx in range(len(lines)):
                line = lines[line_idx]
                line_split = line.split('\t')
                line_data = []
                for data_idx in range(1, len(line_split)):
                    line_data.append(float(line_split[data_idx].replace('\n', '')))
                result.append(line_data)
                indexes.append(line_idx)
        return result, indexes

class KMeansRunner:

    def __init__(self, target_clusters: int = 350, variens: int = 10, tsv_data: Tuple[List, List] = None) -> None:
        self.min_clusters = target_clusters - variens
        self.max_clsuters = target_clusters + variens
        self.data, self.indexes = tsv_data

    def make_tsv(self, cluster_map: pd.DataFrame, nr_name: int) -> None:
        data_dct: Dict[int, List[int]] = {}

        file_path = join(FOLDER_PATH, OUTPUT_NAME + str(nr_name) + OUTPUT_TYPE)

        for idx_, row in cluster_map.iterrows():
            cluster_id = row['cluster']
            data_idx = row['data_index']

            if cluster_id not in data_dct:
                data_dct[cluster_id] = [data_idx]
            else:
                data_dct[cluster_id].append(data_idx)
        
        with open(file_path, 'w') as f:
            for cluster_idx, item_lst in collections.OrderedDict(sorted(data_dct.items())).items():
                for item in item_lst:
                    line = f"{cluster_idx}\tedge_{item+1}\n"
                    f.write(line)
    
    def run(self, target_clusters: int) -> None:
        kmeans = KMeans(n_clusters=target_clusters)
        kmeans.fit(self.data)

        data_map = pd.DataFrame()
        data_map['data_index'] = self.indexes
        data_map['cluster'] = kmeans.labels_

        self.make_tsv(data_map, target_clusters)

    def run_range(self) -> None:
        for x in tqdm(range(self.min_clusters, self.max_clsuters)):
            self.run(x)



if __name__ == "__main__":    
    file_path = join(FOLDER_PATH, FILE_NAME)
    data = TsvLoader.load_tsv(file_path)

    runner = KMeansRunner(target_clusters=389, variens=15, tsv_data=data)
    runner.run_range()

    # kmeans = KMeans(n_clusters=389)
    # kmeans.fit(data)
    # y_kmeans = kmeans.predict(data)


    # make_tsv(data_map, len(indexes))

    print("Done")