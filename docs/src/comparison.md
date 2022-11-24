# Comparison Results 

The results of using `PowerModelsADA`v0.1 on 9 selected test cases with three distributed algorithms and 5 power flow formulations are shown here. We run the three distributed algorithms on a high-performance computing service with a 16-core CPU and 16GB of RAM. The results here used PowerModelsADA v0.1 in Julia v1.6. The results use Ipopt solver for the polar and rectangular ACOPF and Gurobi for the DC approximation, and SOCP and QC relaxations of the OPF.


### Polar AC OPF

| Algorithm   	|      	| ADMM    	|      	| ATC    	|      	| APP     	|      	|
|-------------	|:-----:|---------:	|------:	|--------:	|------:	|---------:	|------:	|
| Case name   	| Area 	|   Time  	| Itr. 	|  Time  	| Itr. 	|   Time  	| Itr. 	|
| 14_ieee     	| 2    	| 0.45    	| 16   	| 0.78   	| 27   	| 0.82    	| 31   	|
| 24_ieee_rts 	| 4    	| 8.09    	| 100  	| 4.56   	| 66   	| 10.80   	| 142  	|
| 30_ieee     	| 3    	| 1.17    	| 23   	| 1.39   	| 27   	| 1.56    	| 31   	|
| 30pwl       	| 3    	| 1.40    	| 23   	| 3.42   	| 53   	| 1.73    	| 26   	|
| 39_epri     	| 3    	| 8.09    	| 89   	| 625.31 	| 198  	| 7.81    	| 94   	|
| 73_ieee_rts 	| 3    	| 11.39   	| 61   	| 10.74  	| 53   	| 13.52   	| 70   	|
| 179_goc     	| 3    	| 33.44   	| 66   	| 363.90 	| 110  	| 72.02   	| 163  	|
| 300_ieee    	| 4    	| 40.33   	| 99   	| 40.25  	| 82   	| 36.15   	| 76   	|
| 588_sdet    	| 8    	| 1818.97 	| 1877 	| 660.10 	| 625  	| 1280.41 	| 1278 	|


### Rectangular AC OPF

| Algorithm   	|       | ADMM    	|      	| ATC   	|      	| APP     	|      	|
|-------------	|:-----: |    |--------:	|------:	|--------:	|------:	|-------:	|------:	|
| Case name   	| Area  |   Time  	| Itr. 	|  Time 	| Itr. 	|   Time  	| Itr. 	|
| 14_ieee     	| 2     | 0.61    	| 19   	| 1.22  	| 28   	| 1.06    	| 32   	|
| 24_ieee_rts 	| 4     | 15.29   	| 133  	| 9.10  	| 77   	| 20.36   	| 205  	|
| 30_ieee     	| 3     | 1.91    	| 27   	| 2.40  	| 33   	| 1.55    	| 23   	|
| 30pwl       	| 3     | 3.28    	| 42   	| 5.27  	| 62   	| 4.57    	| 58   	|
| 39_epri     	| 3     | 50.77   	| 397  	| NC    	| NC   	| 10.17   	| 94   	|
| 73_ieee_rts 	| 3     | 24.94   	| 111  	| 18.92 	| 62   	| 25.23   	| 115  	|
| 179_goc     	| 3     | 93.58   	| 165  	| 32.17 	| 45   	| 86.78   	| 169  	|
| 300_ieee    	| 4     | 40.20   	| 73   	| 82.60 	| 94   	| 47.26   	| 92   	|
| 588_sdet    	| 8     | 1538.00 	| 1337 	| NC    	| NC   	| 1366.12 	| 1180 	|


### DC Approximation 

| Algorithm:   	|       | ADMM   	|      	| ATC    	|      	| APP   	|      	|
|-------------	|:------: |--------:	|------:	|--------:	|------:	|-------:	|------:	|
| Case name   	| Area  |  Time  	| Itr. 	|  Time  	| Itr. 	|  Time 	| Itr. 	|
| 14_ieee     	| 2     |   0.10 	|   14 	|   0.23 	|   34 	|  1.13 	|   20 	|
| 24\_ieee_rts 	| 4     |   3.79 	|  204 	|   1.67 	|   62 	|  3.71 	|  198 	|
| 30_ieee     	| 3     |   0.18 	|   21 	|   0.35 	|   45 	|  0.27 	|   35 	|
| 30pwl       	| 3     |   0.18 	|   22 	|   0.34 	|   42 	|  0.16 	|   19 	|
| 39_epri     	| 3     |   1.21 	|   71 	|  6.181 	|   69 	|  1.38 	|   63 	|
| 73\_ieee_rts 	| 3     |   2.16 	|   74 	|   2.26 	|   68 	|  2.53 	|   78 	|
| 179_goc     	| 3     |   2.85 	|   61 	|   1.01 	|   26 	|  1.94 	|   43 	|
| 300_ieee    	| 4     |   2.38 	|   37 	|   4.73 	|   59 	|  2.43 	|   35 	|
| 588_sdet    	| 8     | 102.48 	| 1168 	| 18.535 	|  655 	| 89.58 	| 1155 	|


### SOCP Relaxation 

| Algorithm   	|        | ADMM    	|      	| ATC   	|      	| APP     	|      	|
|-------------	|:-------: |--------:	|------:	|--------:	|------:	|-------:	|------:	|
| Case name   	| Area   |   Time  	| Itr. 	|  Time 	| Itr. 	|   Time  	| Itr. 	|
| 14_ieee     	| 2      | 0.49    	| 15   	| 0.76  	| 31   	| 0.75    	| 17   	|
| 24_ieee_rts 	| 4      | 2.81    	| 52   	| 2.55  	| 35   	| 3.15    	| 51   	|
| 30_ieee     	| 3      | 0.47    	| 11   	| 1.35  	| 33   	| 1.46    	| 25   	|
| 30pwl       	| 3      | 0.55    	| 11   	| 2.61  	| 46   	| 1.27    	| 15   	|
| 39_epri     	| 3      | 9.53    	| 105  	| 3.92  	| 72   	| 7.42    	| 117  	|
| 73_ieee_rts 	| 3      | 0.16    	| 2    	| 1.52  	| 22   	| 0.19    	| 2    	|
| 179_goc     	| 3      | 3.37    	| 19   	| 3.78  	| 23   	| 5.44    	| 35   	|
| 300_ieee    	| 4      | 11.04   	| 25   	| 20.11 	| 52   	| 8.93    	| 18   	|
| 588_sdet    	| 8      | 3381.35 	| 5000 	| 44.55 	| 79   	| 3177.05 	| 5000 	|


### QC Relaxation 

| Algorithm   	|       | ADMM   	|      	| ATC    	|      	| APP    	|      	|
|-------------	|:----: |--------:|------:	|--------:	|------:	|-------:	|------:	|
| Case name   	| Area  |  Time  	| Itr. 	|  Time  	| Itr. 	|  Time  	| Itr. 	|
| 14_ieee     	| 2     | 0.64   	| 18   	| 0.75   	| 23   	| 0.46   	| 14   	|
| 24_ieee_rts 	| 4     | 5.47   	| 58   	| 3.62   	| 53   	| 5.57   	| 65   	|
| 30_ieee     	| 3     | 0.88   	| 14   	| 1.43   	| 33   	| 1.99   	| 26   	|
| 30pwl       	| 3     | 1.96   	| 39   	| 5.36   	| 79   	| 1.22   	| 12   	|
| 39_epri     	| 3     | 15.74  	| 195  	| NC     	| 1000 	| 10.63  	| 102  	|
| 73_ieee_rts 	| 3     | 2.04   	| 6    	| NC     	| 1000 	| 1.17   	| 6    	|
| 179_goc     	| 3     | 30.46  	| 44   	| 60.45  	| 63   	| 20.83  	| 39   	|
| 300_ieee    	| 4     | 63.74  	| 34   	| NC     	| 1000 	| 54.44  	| 26   	|
| 588_sdet    	| 8     | 478.99 	| 184  	| 220.44 	| 102  	| 455.58 	| 177  	|