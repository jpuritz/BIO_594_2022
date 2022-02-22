###Exome capture

##### Capture
```
fastp --in1 /home/BIO594/DATA/Week4/realdata/exome_capture/Capture1.F.fq --in2 /home/BIO594/DATA/Week4/realdata/exome_capture/Capture1.R.fq --out1 trim.Capture1.F.fq --out2 trim.Capture1.R.fq --cut_by_quality5 20 --cut_by_quality3 20 --cut_window_size 5 --cut_mean_quality 15 -q 15 -u 50 -j Capture.json -h Capture.html 

Read1 before filtering:
total reads: 100000
total bases: 15000000
Q20 bases: 14229580(94.8639%)
Q30 bases: 13142251(87.615%)

Read1 after filtering:
total reads: 98251
total bases: 13814610
Q20 bases: 13200706(95.5561%)
Q30 bases: 12276251(88.8643%)

Read2 before filtering:
total reads: 100000
total bases: 15000000
Q20 bases: 13842504(92.2834%)
Q30 bases: 12603463(84.0231%)

Read2 aftering filtering:
total reads: 98251
total bases: 13791012
Q20 bases: 12875724(93.3632%)
Q30 bases: 11826266(85.7534%)

Filtering result:
reads passed filter: 196502
reads failed due to low quality: 3412
reads failed due to too many N: 86
reads failed due to too short: 0
reads with adapter trimmed: 59180
bases trimmed due to adapters: 1818389

JSON report: Capture.json
HTML report: Capture.html

fastp --in1 /home/BIO594/DATA/Week4/realdata/exome_capture/Capture1.F.fq --in2 /home/BIO594/DATA/Week4/realdata/exome_capture/Capture1.R.fq --out1 trim.Capture1.F.fq --out2 trim.Capture1.R.fq --cut_by_quality5 20 --cut_by_quality3 20 --cut_window_size 5 --cut_mean_quality 15 -q 15 -u 50 -j Capture.json -h Capture.html 
fastp v0.12.4, time used: 0 seconds
```

##### CASE_J03

```
fastp --in1 /home/BIO594/DATA/Week4/realdata/exome_capture/CASE_J03.F.fq.gz --in2 /home/BIO594/DATA/Week4/realdata/exome_capture/CASE_J03.R.fq.gz --out1 trim.CASE_J03.F.fq.gz --out2 trim.CASE_J03.R.fq.gz --cut_by_quality5 20 --cut_by_quality3 20 --cut_window_size 5 --cut_mean_quality 15 -q 15 -u 50 -j CASE_J03.json -h CASE_J03.html 

Read1 before filtering:
total reads: 100000
total bases: 14200000
Q20 bases: 13547660(95.4061%)
Q30 bases: 12612857(88.8229%)

Read1 after filtering:
total reads: 98087
total bases: 12809020
Q20 bases: 12317045(96.1592%)
Q30 bases: 11549216(90.1647%)

Read2 before filtering:
total reads: 100000
total bases: 14200000
Q20 bases: 13431086(94.5851%)
Q30 bases: 12488351(87.9461%)

Read2 aftering filtering:
total reads: 98087
total bases: 12801873
Q20 bases: 12221061(95.4631%)
Q30 bases: 11452288(89.4579%)

Filtering result:
reads passed filter: 196174
reads failed due to low quality: 3636
reads failed due to too many N: 190
reads failed due to too short: 0
reads with adapter trimmed: 68108
bases trimmed due to adapters: 2271633

JSON report: CASE_J03.json
HTML report: CASE_J03.html

fastp --in1 /home/BIO594/DATA/Week4/realdata/exome_capture/CASE_J03.F.fq.gz --in2 /home/BIO594/DATA/Week4/realdata/exome_capture/CASE_J03.R.fq.gz --out1 trim.CASE_J03.F.fq.gz --out2 trim.CASE_J03.R.fq.gz --cut_by_quality5 20 --cut_by_quality3 20 --cut_window_size 5 --cut_mean_quality 15 -q 15 -u 50 -j CASE_J03.json -h CASE_J03.html 
fastp v0.12.4, time used: 2 seconds
```

###Rad-seq

#### FCLib1

```
fastp --in1 /home/BIO594/DATA/Week4/realdata/rad_seq/FCLib1.F.fastq.gz --in2 /home/BIO594/DATA/Week4/realdata/rad_seq/FCLib1.R.fastq.gz --out1 trim.FCLib1.F.fastq.gz --out2 trim.FCLib1.R.fastq.gz --cut_by_quality5 20 --cut_by_quality3 20 --cut_window_size 5 --cut_mean_quality 15 -q 15 -u 50 -j FCLib1.json -h FCLib1.html 

Read1 before filtering:
total reads: 100000
total bases: 15100000
Q20 bases: 14315329(94.8035%)
Q30 bases: 13250625(87.7525%)

Read1 after filtering:
total reads: 98073
total bases: 14758373
Q20 bases: 14080999(95.4102%)
Q30 bases: 13077139(88.6083%)

Read2 before filtering:
total reads: 100000
total bases: 15100000
Q20 bases: 13779401(91.2543%)
Q30 bases: 12201153(80.8023%)

Read2 aftering filtering:
total reads: 98073
total bases: 14728322
Q20 bases: 13604367(92.3688%)
Q30 bases: 12100766(82.1598%)

Filtering result:
reads passed filter: 196146
reads failed due to low quality: 3854
reads failed due to too many N: 0
reads failed due to too short: 0
reads with adapter trimmed: 938
bases trimmed due to adapters: 32855

JSON report: FCLib1.json
HTML report: FCLib1.html

fastp --in1 /home/BIO594/DATA/Week4/realdata/rad_seq/FCLib1.F.fastq.gz --in2 /home/BIO594/DATA/Week4/realdata/rad_seq/FCLib1.R.fastq.gz --out1 trim.FCLib1.F.fastq.gz --out2 trim.FCLib1.R.fastq.gz --cut_by_quality5 20 --cut_by_quality3 20 --cut_window_size 5 --cut_mean_quality 15 -q 15 -u 50 -j FCLib1.json -h FCLib1.html 
fastp v0.12.4, time used: 2 seconds
```

#### WC1_407

```
fastp --in1 /home/BIO594/DATA/Week4/realdata/rad_seq/WC1_407.F.fq.gz --in2 /home/BIO594/DATA/Week4/realdata/rad_seq/WC1_407.R.fq.gz --out1 trim.WC1_407.F.fq.gz --out2 trim.WC1_407.R.fq.gz --cut_by_quality5 20 --cut_by_quality3 20 --cut_window_size 5 --cut_mean_quality 15 -q 15 -u 50 -j WC1_407.json -h WC1_407.html 

Read1 before filtering:
total reads: 100000
total bases: 14600000
Q20 bases: 13989497(95.8185%)
Q30 bases: 13179453(90.2702%)

Read1 after filtering:
total reads: 98651
total bases: 14365378
Q20 bases: 13827118(96.2531%)
Q30 bases: 13054332(90.8736%)

Read2 before filtering:
total reads: 100000
total bases: 15100000
Q20 bases: 14126787(93.5549%)
Q30 bases: 12951507(85.7716%)

Read2 aftering filtering:
total reads: 98651
total bases: 14829691
Q20 bases: 14001419(94.4148%)
Q30 bases: 12876295(86.8278%)

Filtering result:
reads passed filter: 197302
reads failed due to low quality: 2698
reads failed due to too many N: 0
reads failed due to too short: 0
reads with adapter trimmed: 552
bases trimmed due to adapters: 17250

JSON report: WC1_407.json
HTML report: WC1_407.html

fastp --in1 /home/BIO594/DATA/Week4/realdata/rad_seq/WC1_407.F.fq.gz --in2 /home/BIO594/DATA/Week4/realdata/rad_seq/WC1_407.R.fq.gz --out1 trim.WC1_407.F.fq.gz --out2 trim.WC1_407.R.fq.gz --cut_by_quality5 20 --cut_by_quality3 20 --cut_window_size 5 --cut_mean_quality 15 -q 15 -u 50 -j WC1_407.json -h WC1_407.html 
fastp v0.12.4, time used: 1 seconds
```

### RNA-seq

```
fastp --in1 /home/BIO594/DATA/Week4/realdata/rna_seq/rna1.F.fq.gz --in2 /home/BIO594/DATA/Week4/realdata/rna_seq/rna1.R.fq.gz --out1 trim.rna1.F.fq.gz --out2 trim.rna1.R.fq.gz --cut_by_quality5 20 --cut_by_quality3 20 --cut_window_size 5 --cut_mean_quality 15 -q 15 -u 50 -j rna1.json -h rna1.html 

Read1 before filtering:
total reads: 1000
total bases: 101000
Q20 bases: 92490(91.5743%)
Q30 bases: 84200(83.3663%)

Read1 after filtering:
total reads: 773
total bases: 75285
Q20 bases: 73879(98.1324%)
Q30 bases: 67620(89.8187%)

Read2 before filtering:
total reads: 1000
total bases: 101000
Q20 bases: 83420(82.5941%)
Q30 bases: 74022(73.2891%)

Read2 aftering filtering:
total reads: 773
total bases: 74404
Q20 bases: 72587(97.5579%)
Q30 bases: 65042(87.4173%)

Filtering result:
reads passed filter: 1546
reads failed due to low quality: 336
reads failed due to too many N: 0
reads failed due to too short: 118
reads with adapter trimmed: 0
bases trimmed due to adapters: 0

JSON report: rna1.json
HTML report: rna1.html

fastp --in1 /home/BIO594/DATA/Week4/realdata/rna_seq/rna1.F.fq.gz --in2 /home/BIO594/DATA/Week4/realdata/rna_seq/rna1.R.fq.gz --out1 trim.rna1.F.fq.gz --out2 trim.rna1.R.fq.gz --cut_by_quality5 20 --cut_by_quality3 20 --cut_window_size 5 --cut_mean_quality 15 -q 15 -u 50 -j rna1.json -h rna1.html 
fastp v0.12.4, time used: 0 seconds
```

### WGS

```
fastp --in1 /home/BIO594/DATA/Week4/realdata/wgs/HC_1.F.fq.gz --in2 /home/BIO594/DATA/Week4/realdata/wgs/HC_1.R.fq.gz --out1 trim.HC_1.F.fq.gz --out2 trim.HC_1.R.fq.gz --cut_by_quality5 20 --cut_by_quality3 20 --cut_window_size 5 --cut_mean_quality 15 -q 15 -u 50 -j HC_1.json -h HC_1.html 

Read1 before filtering:
total reads: 100000
total bases: 15100000
Q20 bases: 14454860(95.7275%)
Q30 bases: 13562956(89.8209%)

Read1 after filtering:
total reads: 97534
total bases: 14633102
Q20 bases: 14064283(96.1128%)
Q30 bases: 13236705(90.4573%)

Read2 before filtering:
total reads: 100000
total bases: 15100000
Q20 bases: 13873553(91.8778%)
Q30 bases: 12588578(83.3681%)

Read2 aftering filtering:
total reads: 97534
total bases: 14608236
Q20 bases: 13615657(93.2053%)
Q30 bases: 12429426(85.0851%)

Filtering result:
reads passed filter: 195068
reads failed due to low quality: 4932
reads failed due to too many N: 0
reads failed due to too short: 0
reads with adapter trimmed: 8180
bases trimmed due to adapters: 203785

JSON report: HC_1.json
HTML report: HC_1.html

fastp --in1 /home/BIO594/DATA/Week4/realdata/wgs/HC_1.F.fq.gz --in2 /home/BIO594/DATA/Week4/realdata/wgs/HC_1.R.fq.gz --out1 trim.HC_1.F.fq.gz --out2 trim.HC_1.R.fq.gz --cut_by_quality5 20 --cut_by_quality3 20 --cut_window_size 5 --cut_mean_quality 15 -q 15 -u 50 -j HC_1.json -h HC_1.html 
fastp v0.12.4, time used: 2 seconds
```


