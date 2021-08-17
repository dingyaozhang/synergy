# synergy
Our software, synergy, calculates the interaction among different genes in the cancer background. Theoretically, it could also be used to understand the synergy in the context of other complex diseases.

When calculating the synergy based on entropy, there are two different modes:
Mode 1: Continuous variable model. 
Mode 2: Discrete variable model.
The software also has a mode to calculate synergy based on PPI network:
Mode 3: PPI network model.

## Mode 1
### Quick Start
The input file should be matrix (tab-separated values file, tsv) file. In this file, the columns are features (like mutation), and the rows are samples. 
The last column is the result feature (like cancer/normal) we want to study the influence of features on the result feature. In the Mode 1, all features except result feature 
are continuous variable. A name (tab-separated values file, tsv) file with two colunms are needed, too. The first column is gene name, and the second row is the No. of genes.  
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

```
from synergycon import synergycon
synergycon('matrix.txt',  'output.txt', 'matrixname.txt')
synergycon2nd('matrix.txt','temp_output.txt','second_output.txt','matrixname.txt', digitnumber=2,iternum=0):

```
### Options
```
synergycon(datain,outputfile,dataanno, digitnumber=2,iternum=20, limitinter=5)
synergycon2nd(datain,dataresult,outputfile,dataanno, digitnumber=2,iternum=0):

digitnumber:

iternum:

limitinter:

```
