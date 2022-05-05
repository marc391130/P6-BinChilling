from typing import Dict, List, Tuple

class SCGReader:
    def __init__(self, filepath: str) -> None:
        self.filepath = filepath

    def read_scg(self) -> Dict[str, List[str]]:

        temp: Dict[str, List[str]] = {}

        with open(self.filepath, 'r') as f:
            lines = f.readlines()

            for line in lines:
                split_line = line.split('\t')
                name, scg = split_line[0], split_line[1]
                contig_name = name.split('_')[0]
                scg_name = scg.replace('\n', '')

                # print(contig_name, scg_name)

                _ = temp.get(contig_name, [])
                _.append(scg_name)
                temp[contig_name] = _

        result_dct = {item: set(value) for item, value in temp.items()}
        return result_dct

if __name__ == '__main__':
    scg_reader = SCGReader('../Dataset/_proteins.faa.bacteria.scg')
    data = scg_reader.read_scg()

    for item, value in data.items():
        print(item,value)
