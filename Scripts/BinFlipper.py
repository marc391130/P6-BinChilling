from os import listdir
from os.path import join
from typing import Callable, List
import sys



class bin_flipper:
    def __init__(self, input_folder_path: str, output_folder_path: str, file_predicate: Callable[[str], bool], cluster_renamer: Callable[[str, str], str]) -> None:
        self.input_folder_path = input_folder_path
        self.output_folder_path = output_folder_path
        self.file_predicate = file_predicate
        self.cluster_renamer = cluster_renamer
        
    def __get_folder_content(self, folder_path: str) -> List[str]:
        return [f for f in listdir(self.input_folder_path) if self.file_predicate(f)]

    def work(self):
        files = self.__get_folder_content(self.input_folder_path)
        
        for file in files:
            result = []
            with open(join(self.input_folder_path, file), 'r') as f:
                lines = f.readlines()
                for line in lines:
                    temp = line.replace('\n', '').split('\t')
                    temp.reverse()
                    temp[1] = self.cluster_renamer(temp[1], file)
                    result.append('\t'.join(temp) + '\n')
            with open(join(self.output_folder_path, file), 'x') as f:
                f.writelines(result)
                
                
if __name__ == '__main__':
    print(sys.argv[1], sys.argv[2])
    flipper = bin_flipper(sys.argv[1], sys.argv[2], \
        lambda x: x.startswith('cluster') , lambda x,y : y.replace('.tsv', '') + ':' + x)
    flipper.work()