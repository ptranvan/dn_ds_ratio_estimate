## Introduction

dN/dS ratio estimates of mites using [Codeml (PAML)](http://abacus.gene.ucl.ac.uk/software/paml.html).

The input data don't work with Codeml.

Alignments are in [phylip](http://evolution.genetics.washington.edu/phylip.html) format.

Trees have been generated with [RAxML](http://sco.h-its.org/exelixis/software.html) but have to be reformatted depending on the model.

```mites2codeml.py``` is useful to reformat the input data to make it compatible with [Codeml (PAML)](http://abacus.gene.ucl.ac.uk/software/paml.html).
 
## Table of contents

1. [Required files](#1_input)
2. [dN/dS ratio estimates](#2_paml)

## <a name="1_input"></a>1) Required input

1. Need 2 directory and 1 control file:

	* alignments/ directory: contain alignments in phylip format.
	* branch_lengths/ directory: contain the RAxML's trees in newick format.
	* control template (available here codeml_[model].ctl)

## <a name="2_paml"></a>2) dN/dS ratio estimates

1) Create a directory for your species (mites).

```
mkdir mites
```

Copy the 3 required input into this directory.

```
ll mites

total 9
drwxr-xr-x 2 ptranvan dee_schwan 3545  2 mar 11:00 alignments
drwxr-xr-x 2 ptranvan dee_schwan 3545  2 mar 11:01 branch_lengths
-rw-r--r-- 1 ptranvan dee_schwan 2931  2 mar 11:00 codeml_[model].ctl
```

2.1) For SAI and TI models:

```
python mites2codeml.py -s sai_ti -i1 alignments -i2 branch_lengths -i3 codeml_sai_ti.ctl -e <your_email> -o sai_ti.sh
```

2.2) For Hydrophobicity model:

```
python mites2codeml.py -s hydrophobicity -i1 alignments -i2 branch_lengths -i3 codeml_hydrophobicity.ctl -e <your_email> -o hydrophobicity.sh
```

3) A script (*.sh) will be created. It has to be run on Vital-IT server (usr@prd.vital-it.ch):

```
chmod +x *.sh
./*.sh
```

Once done, a directory output will be created.

4) Options:

```
-i1: PATH of alignments directory.
-i2: PATH of branch_lengths directory.
-i3: Control file.
-e: Email for job confirmation.
-o: Name of the script file to be executed.
```

