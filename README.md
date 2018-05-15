# countingBarcodeInFastq
Counting multiple DNA barcodes in gzipped fastq 

Scanning a fastq.gz file to count DNA barcodes/tag pre-defined in an Excel file.            
                                                                                
This module aims to count occurrence of prefined short DNA sequences            
(barcode/tag) and its combinations in a fastq.gz file from Next Generation      
Sequencing.                                                                     
                                                                                
The input files include a .fastq.gz file and barcode sets in Excel file with    
specific format (-b). Other required arguments are library length (insert_size, 
-s), direction (-d, 0 - forward, 1 - backward, 2 - both), type of barcode being 
put at column in output file (-c), type of barcode being put at row (-r),       
number of processes (-p) and number of reads to be scanned (-n).                
                                                                                
The output file is an Excel file with 4 sheets: summary, barcode_count,         
barcode_Read-Per-Million, and BarPlot.
