# Comparison Results

The results of using `PowerModelsADA`v0.1.1 on 9 test cases from **[PGLib-OPF](https://github.com/power-grid-lib/pglib-opf)** is shown here. We benchmark three distributed algorithms with 5 power flow formulations.

### Simulation Setup

We run the three distributed algorithms on a high-performance computing service with a 16-core CPU and 16GB of RAM. We produced the results shown here using `Julia` v1.8, and [`Ipopt`](https://github.com/jump-dev/Ipopt.jl) solver.

We report the results that achieved the $l_2$-norm of the mismatches less than 0.01 (radians and per unit) within 10,000 iterations and the absolute value of the relative error less than 1% of the central solution from `PowerModels.jl`.

We tune the ADAs parameters by selecting large values and gradually reducing the values until reaching a good solution. For the ADMM and APP, we started with $\alpha = 10^6$ and divided by 10 for the next run, while for the ATC we started with $\alpha =1.2$ and subtracted 0.099 for the next run.

### Polar Form ACOPF

| **Algorithm** |      	| **ADMM**  |      	| **ATC**   |      	| **APP**   |      	|
|-------------	|:-----:|---------:	|------:|--------:	|------:|---------:	|------:|
| **Case name** |**Area**| **Time** |**Itr.**| **Time** |**Itr.**| **Time** |**Itr.**|
| 14_ieee     	| 2    	| 1.13    	| 14   	| 1.70   	| 28   	| 5.36    	| 31   	|
| 24\_ieee_rts 	| 4    	| 21.02    	| 97  	| 7.82   	| 67   	| 37.84   	| 207  	|
| 30_ieee     	| 3    	| 2.18    	| 24   	| 3.89   	| 43   	| 2.33    	| 25   	|
| 30pwl       	| 3    	| 7.40    	| 24   	| 4.08   	| 36   	| 10.73    	| 49   	|
| 39_epri     	| 3    	| 20.14    	| 89   	| 239.80 	| 1261  | 179.63   	| 873  	|
| 73\_ieee_rts 	| 3    	| 14.37   	| 61   	| 18.61  	| 58   	| 23.02   	| 83   	|
| 179_goc     	| 3    	| 31.42   	| 66   	| 62.06 	| 81  	| 76.82   	| 166  	|
| 300_ieee    	| 4    	| 21.51   	| 66   	| 651.26  	| 920  	| 28.16   	| 77   	|
| 588_sdet    	| 8    	| 295.05 	| 871 	| 3133.82 	| 1971 	| 437.17 	| 1283 	|

### Rectangular Form ACOPF

| **Algorithm** |      	| **ADMM**  |      	| **ATC**   |      	| **APP**   |      	|
|-------------	|:-----:|---------:	|------:|--------:	|------:|---------:	|------:|
| **Case name** |**Area**| **Time** |**Itr.**| **Time** |**Itr.**| **Time** |**Itr.**|
| 14_ieee     	| 2     | 1.04    	| 14   	| 1.75  	| 28   	| 2.80    	| 33   	|
| 24\_ieee_rts 	| 4     | 11.45   	| 95  	| 8.19  	| 66   	| 23.70   	| 206  	|
| 30_ieee     	| 3     | 3.20    	| 22   	| 3.53  	| 40   	| 3.36    	| 24   	|
| 30pwl       	| 3     | 1.07    	| 10   	| 6.82  	| 59   	| 6.13    	| 48   	|
| 39_epri     	| 3     | 11.40   	| 86  	| 323.25   	| 1243  | 12.12   	| 95   	|
| 73\_ieee_rts 	| 3     | 10.38   	| 60  	| 21.85 	| 58   	| 19.76   	| 116  	|
| 179_goc     	| 3     | 47.25   	| 122  	| 63.58 	| 83   	| 64.74   	| 170  	|
| 300_ieee    	| 4     | 17.89   	| 44   	| 1088.83 	| 900  	| 33.34   	| 94   	|
| 588_sdet    	| 8     | 401.70 	| 1031 	| 3838.72  	| 1977 	| 473.00 	| 1181 	|

### DC Approximation

| **Algorithm** |      	| **ADMM**  |      	| **ATC**   |      	| **APP**   |      	|
|-------------	|:-----:|---------:	|------:|--------:	|------:|---------:	|------:|
| **Case name** |**Area**| **Time** |**Itr.**| **Time** |**Itr.**| **Time** |**Itr.**|
| 14_ieee     	| 2     |   1.24 	|   15 	|   0.79 	|   45 	|  1.24 	|   21 	|
| 24\_ieee_rts 	| 4     |   12.02 	|  174 	|   2.61 	|   60 	|  14.05 	|  199 	|
| 30_ieee     	| 3     |   1.41 	|   21 	|   0.35 	|   45 	|  1.44 	|   20 	|
| 30pwl       	| 3     |   1.53 	|   18 	|   0.34 	|   42 	|  1.38 	|   20 	|
| 39_epri     	| 3     |   4.90 	|   69 	|  6.181 	|   69 	|  4.51 	|   64 	|
| 73\_ieee_rts 	| 3     |   6.81 	|   75 	|   4.88 	|   55 	|  6.77 	|   79 	|
| 179_goc     	| 3     |   4.56 	|   37 	|   3.24 	|   27 	|  5.21 	|   44 	|
| 300_ieee    	| 4     |   3.60 	|   26 	|   11.31 	|   58 	|  6.03 	|   36 	|
| 588_sdet    	| 8     | 106.01 	| 656 	| 18.535 	|  655 	| 185.20 	| 1156 	|

### SOCP Relaxation

| **Algorithm** |      	| **ADMM**  |      	| **ATC**   |      	| **APP**   |      	|
|-------------	|:-----:|---------:	|------:|--------:	|------:|---------:	|------:|
| **Case name** |**Area**| **Time** |**Itr.**| **Time** |**Itr.**| **Time** |**Itr.**|
| 14_ieee     	| 2      | 1.29    	| 12   	| 2.16  	| 39   	| 0.98    	| 12   	|
| 24\_ieee_rts 	| 4      | 3.85    	| 39   	| 4.18  	| 32   	| 4.94    	| 51   	|
| 30_ieee     	| 3      | 1.32    	| 12   	| 3.96  	| 33   	| 1.46    	| 13   	|
| 30pwl       	| 3      | 1.30    	| 11   	| 4.53  	| 29   	| 1.42    	| 13   	|
| 39_epri     	| 3      | 5.34    	| 41  	| 8.53  	| 47   	| 16.12    	| 119  	|
| 73\_ieee_rts 	| 3      | 0.51    	| 3    	| 8.03  	| 23   	| 0.45    	| 3    	|
| 179_goc     	| 3      | 7.56    	| 16   	| 20.87  	| 23   	| 3.54    	| 8   	|
| 300_ieee    	| 4      | 5.64   	| 11   	| 51.75 	| 42   	| 7.18    	| 14   	|
| 588_sdet    	| 8      | 87.89 	| 131 	| 145.32 	| 53   	| 85.91 	| 130 	|

### QC Relaxation

| **Algorithm** |      	| **ADMM**  |      	| **ATC**   |      	| **APP**   |      	|
|-------------	|:-----:|---------:	|------:|--------:	|------:|---------:	|------:|
| **Case name** |**Area**| **Time** |**Itr.**| **Time** |**Itr.**| **Time** |**Itr.**|
| 14_ieee     	| 2     | 2.09   	| 15   	| 4.37   	| 24   	| 4.47   	| 15   	|
| 24\_ieee_rts 	| 4     | 15.53   	| 55   	| 22.07   	| 50   	| 19.10   	| 66   	|
| 30_ieee     	| 3     | 11.05   	| 32   	| 12.66   	| 33   	| 9.89   	| 30   	|
| 30pwl       	| 3     | 1.73   	| 9   	| 11.25   	| 29   	| 2.68   	| 14   	|
| 39_epri     	| 3     | 31.80  	| 85  	| 82.02    	| 104 	| 36.31  	| 101  	|
| 73\_ieee_rts 	| 3     | 2.59   	| 7    	| 1268.75  	| 1009 	| 2.81   	| 7    	|
| 179_goc     	| 3     | 49.63  	| 24   	| 118.34  	| 24   	| 77.16  	| 39   	|
| 300_ieee    	| 4     | 18.52  	| 11   	| 602.90   	| 60 	| 57.31  	| 27   	|
| 588_sdet    	| 8     | 321.28 	| 177  	| 500.23 	| 56  	| 316.85 	| 177  	|
