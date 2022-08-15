# Code for Finding Atomic Cliques in Temporal Graphs


This code accompanies the paper "On Atomic Cliques in Temporal Graphs" and is written in C++. If you wish to use or cite this code, please cite the paper: 

    @article{LMSB2022Atomic-clique, 
  	title = {On Atomic Cliques in Temporal Graphs}, 
  	author = {Lu, Yajun and Miao, Zhuqi and Sahraeian, Parisa and Balasundaram, Balabhaskar}, 
	note = {Under Review},
  	year = {2022}
	}



## Compiling the code
The following steps show how to compile and run the code to find atomic cliques in a Linux environment using a makefile (you can also run the code in Mac or Windows environment by configuring your IDE appropriately). 

- "parameter.txt" file is used to configure testbed folder and type of solver
- "data" folder includes a subset of graph instances that are used in this study. Complete graph instances can be downloaded from [Mendeley Data](https://data.mendeley.com/datasets/sj2jzztv4z/1).


### Steps to run the code to find atomic cliques solver in Linux environment:
1. Download or clone the repository to your machine.
2. Open the "Makefile" and set GUROBI_HOME to the directory of your Gurobi installation, e.g.: /opt/gurobi/9.5.1/linux64.
3. From the terminal, go to the folder containing "Makefile".
4. Type "make" and hit enter to compile. 
5. Solver selection through "parameter.txt" file:
- If you would like to run Enhanced Formulation (EF) solver, set solver type =1.
- If you would like to run Edge Peeling (EP) + [WB solver](https://github.com/jwalteros/dOmega), set solver type =2. After Step 6, the code will run EP algorithm to generate auxiliary graphs as WB solver's input. Then run the [WB solver](https://github.com/jwalteros/dOmega) on these auxiliary graphs.
6. Type "./main" to run the code.


## Acknowledgments
We would like to thank Jose Walteros and Austin Buchanan for making their [solver](https://github.com/jwalteros/dOmega) publicly available.


## Terms and Use:

MIT License

Copyright (c) 2022 Yajun Lu, Zhuqi Miao, Parisa Sahraeian, and Balabhaskar Balasundaram.

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
