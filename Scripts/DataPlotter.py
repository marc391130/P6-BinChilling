import matplotlib.pyplot as plt
import argparse
from typing import List, Dict, Tuple

def parse_args() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
            prog='DataPlotter',
            description="""Data Plotting""",
            formatter_class=argparse.RawDescriptionHelpFormatter,
            usage="%(prog)s WRITE THIS LATER PLZ",
            add_help=True
            )
        
    p_args = parser.add_argument_group(title='Input Arguments (required)', description=None)
    p_args.add_argument('--data' ,'-d', metavar='', required=True,\
        dest='data', help='path to data file')
    p_args.add_argument('--ignore' ,'-ig', metavar='', nargs='+', type=int, required=False,\
        dest='ignore', help='Collumns to ignore!', default = [])

    return parser.parse_args()

def get_data(args) -> Dict[str, List[float]]:

    def __get_line_data__(line: str, ignore_set, convert_to_float = True) -> List[str]:
        result = []
        split_line = line.split('\t')

        for idx in range(len(split_line)):
            if idx in ignore_set: continue
            segment = split_line[idx].replace('\n', '')

            try:
                if convert_to_float: result.append(float(segment))
                else: result.append(segment)
            except Exception as e:
                print(e)
                print(f"Data: {segment} is misformatted!")
        return result

    ignore_set = args.ignore
    ignore_set = [ignore_set] if isinstance(ignore_set, int) else ignore_set

    print(f"Parsing data from file: {args.data}")
    with open(args.data) as f:
        lines = f.readlines()
        titles = __get_line_data__(lines[0], ignore_set, convert_to_float=False)
        result: Dict[str, List[float]] = {title: [] for title in titles}
        title_map: Dict[int, str] = {i: titles[i] for i in range(len(titles))}

        for idx in range(1, len(lines)):
            line = lines[idx]
            data = __get_line_data__(line, ignore_set)

            for idx_d in range(len(data)):
                data_seg = data[idx_d]
                result[title_map[idx_d]] += [data_seg]
    
    return result

def plot(data: Dict[str, List[float]]) -> None:

    for title, data_lst in data.items():
        x_data = [i for i in range(len(data_lst))]
        y_data = data_lst

        plt.plot(x_data, y_data, label=title)
    plt.legend()
    plt.show()


if __name__ == '__main__':
    args = parse_args()
    data = get_data(args)
    plot(data)