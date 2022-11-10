!Bin Chilling: A Metagenomic Ensemble Method

Docker file is also avaiable, and is the prefered way to run Bin Chilling.
only tested on linux.

Quickstart guide:

1. python3 setup.py
3. python3 BinChilling -h
4. follow guide


!Data requirements:


The data requirements of BinChillign are the following:
- Fasta file
- Depth/abundance file (when standable mode)
- marker gene stats
    - file MUST end in .tsv
    - alternativly, a .scg file can be used, with format "EDGE_NAME\tSCG"
- input partitions (unless in standalone mode)


!Docker Example:
Data used in this example is the strong100 dataset, which can be found at:
https://zenodo.org/record/6122610
Or can be directly downloaded from:
https://zenodo.org/record/6122610/files/strong100.zip?download=1

Example assumes the contents of the dataset is located in the Data folder. 

Add your dataset files to 
run following command in folder with dockerfile:
```
sudo docker build -t bin_chilling --no-cache .
```

Run docker container as using:
```
sudo docker run -it --rm -v /YOURPATH/P6/Data/:/P6/Data/ bin_chilling
```
This should start a docker container in interactive mode. 
Navigate to inside the BinChilling folder (Must be inside this folder, as the internal code expects the origin to be there.)
```
cd ./BinChilling
```
A help menu can be displayed using:
```
python3 . -h
```
Running in Ensemble mode:
```
python3 . -M Ensemble -f ../Data/edges.fasta -p ../Data/Partitions/ -g ../Data/marker_gene_stats.tsv -o ../Data/ens_output.tsv -j ../Data/edges_depth.txt --log ../Data/log.txt -ch 50 -am 0.90
```

Running in standalone mode:
```
python3 . -M Bin -f ../Data/edges.fasta -p ../Data/ClusterData/ -g ../Data/marker_gene_stats.tsv -o ../Data/bin_output.tsv -j ../Data/edges_depth.txt --log ../Data/log.txt -ch 50 -am 0.90
```

A file matching the output partition path file alongside one appended with .ref.tsv, will be outputted.
The .ref.tsv file is the refined version, which is the optimized output.

```
-f = fasta file
-p = partition path
-g = marker genes
-o = output partition path
-j = depth file
--log = log file path (can be omitted)
-ch = chunksize (performance impact, does not change result)
-am = a1min value
```


!Evaluation
Evaluation is done using marker_gene_stats the formulas from DAS_TOOls paper.
Navigate to the scripts folder.
python3 EvaluateBins.py -f ../Data/edges.fasta -g ../Data/marker_gene_stats.tsv  -p ../Data/Bin_output.tsv.ref.tsv -d ../Data/edges_depth.txt -o ../Data/evaluation.tsv