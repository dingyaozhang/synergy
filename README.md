# synergy
Our software, synergy, calculates the interaction among different genes in the cancer background. Theoretically, it could also be used to understand the synergy in the context of other complex diseases.

When calculating the synergy based on entropy, there are two different modes:
Mode 1: Continuous variable model. 
Mode 2: Discrete variable model.
The software also has a mode to calculate synergy based on PPI network:
Mode 3: PPI network model.

## Mode 1
### Quick Start
The input file should be matrix (tab-separated values file, tsv) file. In this file, the columns are features (like gene exoression level), and the rows are samples. 
The last column is the result feature (like cancer/normal) we want to study the influence of features on the result feature. In the Mode 1, all features except result feature 
are continuous variable. A name (tab-separated values file, tsv) file with two colunms are needed, too. The first column includes gene name, and the second row includes gene IDs. The name file is corresponding to each feature in matrix file. 

The matrix file looks like (only the numeric part is in the file.): 
|          | *feature 1*|  *feature 2*|*...* | *Feature N* | *Outcome* |
|  ----    |    ----  | ---- | ----|----  | ---- |
| *Sample 1* | 1.12     | 7.89 |  ...| 3.21 | 1    |
| *Sample 2* | 3.45     | 1.12 |  ...| 2.31 | 0    |
| ...      | ...      | ...  | ... | ...  |...   |
| *Sample M* | 13.45    | 21.12| ... | 9.87 | 0    |

(An example of this file is in.)

The name file looks like (The first row is not in the file.): 

|  *Gene name*   | *Gene id*  |
| ----  | ---- |
| ABCA13  |	ENSP00000411096 |
| ABCA3 	| ENSP00000301732 |
| ...     | ...             |
| ABCC9   |	ENSP00000261200 |

(An example of this file is in.)

There are two functions in mode 1: synergycon and synergycon2nd.
The difference is that synergycon2nd also considers another output file of synergy pipeline, which enables an combination of all three modes.
```
from synergycon import synergycon
synergycon.synergycon('matrix.txt',  'output.txt', 'matrixname.txt')
synergycon.synergycon2nd('matrix.txt','temp_output.txt','second_output.txt','matrixname.txt', digitnumber=3,iternum=0):

```
### Options
```

digitnumber:
How many die

iternum:
The number of iteration used to calculate p-value.

limitinter:
How many features could be included in one interaction.

```

## Mode 2
### Quick Start
The input file should be matrix (tab-separated values file, tsv) file. In this file, the columns are features (like mutation), and the rows are samples. 
The last column is the result feature (like cancer/normal) we want to study the influence of features on the result feature. In the Mode 2, all features except result feature 
are discrete variable. A name (tab-separated values file, tsv) file with two colunms are needed, too. The first column includes gene name, and the second row includes gene IDs.
The name file is corresponding to each feature in matrix file. 

The matrix file looks like (only the numeric part is in the file.): 
|          | *feature 1*|  *feature 2*|*...* | *Feature N* | *Outcome* |
|  ----    |    ----  | ---- | ----|----  | ---- |
| *Sample 1* | 11     | 7 |  ...| 3 | 1    |
| *Sample 2* | 3     | 1 |  ...| 23 | 0    |
| ...      | ...      | ...  | ... | ...  |...   |
| *Sample M* | 13    | 21| ... | 9 | 0    |

(An example of this file is in.)

The name file looks like (The first row is not in the file.): 

|  *Gene name*   | *Gene id*  |
|     ----       | ----       |
|ABCA13 |	ENSP00000411096 |
|ABCA3	 | ENSP00000301732 |
| ...    | ...             |
|ABCC9   |	ENSP00000261200 |

(An example of this file is in.)

There are two functions in mode 1: synergycon and synergycon2nd.
The difference is that synergycon2nd also considers another output file of synergy pipeline, which enables an combination of all three modes.
```
from synergydrv import synergydrv
synergydrv.synergydrv('matrix.txt',  'output.txt', 'matrixname.txt')
synergydrv.synergydrv2nd('matrix.txt','temp_output.txt','second_output.txt','matrixname.txt', digitnumber=3,iternum=0):

```
### Options
```
iternum:
The number of iteration used to calculate p-value.

limitinter:
How many features could be included in one interaction.
```


## Mode 3
### Quick Start
This mode needs several input files:
ppi_matrix:
A tsv (tab-separated values file) file with nrows and ncols， it should be a symmetric matrix.
(An example of this file is in.)

A name file which is corresponding to each feature in ppi_matrix file.
The name file looks like (The first row is not in the file.): 
| *Gene id*  |
| ---- |
|	ENSP00000411096 |
| ENSP00000301732 |
| ...             |
|	ENSP00000261200 |
(An example of this file is in.)

A geneexp file including gene expression level which is corresponding to each feature in ppi_matrix/name file.
The geneexp file looks like (The first row is not in the file.): 
| *Gene id*  |
| ---- |
|	1.23 |
| 4.56 |
| ...  |
|	7.89 |
(An example of this file is in.)





There are two functions in mode 1: GCana and GCjudge.
The difference is that GCjudge also considers another output file of synergy pipeline, which enables an combination of all three modes.
```
from synergy import sygyPPI
sygyPPI.GCana('ppi_matrix.txt','name.txt','geneexp.txt', 'output.txt')
sygyPPI.GCjudge('ppi_matrix.txt','name.txt','geneexp.txt','output.txt','temp_output.txt')

```
### Options
```
dataselect:
the genes which are selected to be calculated. If this value is not input, all the genes in the ppi_matrix will be calculated.

iternum:
The number of iteration used to calculate p-value.

limitinter:
How many features could be included in one interaction.

threshold:
the cutoff of ppi giant center.

repnum: 
the maximun number of random walk

changethreshold:
If the change is smaller than this value, the random walk will stop.

alpha:
the alpha parameter in the random walk

```

