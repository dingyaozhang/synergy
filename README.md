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
are continuous variable. A name (tab-separated values file, tsv) file with two colunms are needed, too. The first column includes gene name, and the second row includes gene IDs. 

The matrix file looks like (only the numeric part is in the file.): 
|          | *feature 1*|  *feature 2*|*...* | *Feature N* | *Outcome* |
|  ----    |    ----  | ---- | ----|----  | ---- |
| *Sample 1* | 1.12     | 7.89 |  ...| 3.21 | 1    |
| *Sample 2* | 3.45     | 1.12 |  ...| 2.31 | 0    |
| ...      | ...      | ...  | ... | ...  |...   |
| *Sample M* | 13.45    | 21.12| ... | 9.87 | 0    |

(An example of this file is in.)

The name file looks like (The first row is not in the file.): 

|  Gene name   | Gene id  |
| ----  | ---- |
|ABCA13 |	ENSP00000411096 |
|ABCA3	| ENSP00000301732 |
| ...   | ...             |
|ABCC9  |	ENSP00000261200 |

(An example of this file is in.)

There are two functions in mode 1: synergycon and synergycon2nd.
The difference is that synergycon2nd also considers another output file of synergy pipeline, which enables an combination of all three modes.
```
from synergycon import synergycon
synergycon('matrix.txt',  'output.txt', 'matrixname.txt')
from synergycon import synergycon2nd
synergycon2nd('matrix.txt','temp_output.txt','second_output.txt','matrixname.txt', digitnumber=3,iternum=0):

```
### Options
```
synergycon(datain,outputfile,dataanno, digitnumber=2,iternum=20, limitinter=5)
synergycon2nd(datain,dataresult,outputfile,dataanno, digitnumber=2,iternum=0):

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
The matrix file looks like (only the numeric part is in the file.): 
|          | *feature 1*|  *feature 2*|*...* | *Feature N* | *Outcome* |
|  ----    |    ----  | ---- | ----|----  | ---- |
| *Sample 1* | 11     | 7 |  ...| 3 | 1    |
| *Sample 2* | 3     | 1 |  ...| 23 | 0    |
| ...      | ...      | ...  | ... | ...  |...   |
| *Sample M* | 13    | 21| ... | 9 | 0    |

(An example of this file is in.)

The name file looks like (The first row is not in the file.): 

|  Gene name   | Gene id  |
| ----  | ---- |
|ABCA13 |	ENSP00000411096 |
|ABCA3	| ENSP00000301732 |
| ...   | ...             |
|ABCC9  |	ENSP00000261200 |

(An example of this file is in.)

There are two functions in mode 1: synergycon and synergycon2nd.
The difference is that synergycon2nd also considers another output file of synergy pipeline, which enables an combination of all three modes.
```
from synergydrv import synergydrv
synergydrv('matrix.txt',  'output.txt', 'matrixname.txt')
from synergydrv import synergydrv2nd
synergydrv2nd('matrix.txt','temp_output.txt','second_output.txt','matrixname.txt', digitnumber=3,iternum=0):

```
### Options
```
synergycon(datain,outputfile,dataanno,iternum=20, limitinter=5)
synergycon2nd(datain,dataresult,outputfile,dataanno, digitnumber=2,iternum=0):

iternum:
The number of iteration used to calculate p-value.

limitinter:
How many features could be included in one interaction.


```
