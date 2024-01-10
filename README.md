# cfTFBS-DL
CNN plus LSTM for cfDNA sequences TFBS prediction
We trained the CNN + LSTM model for TFBS prediction of cfDNA sequences based on the short(40bp-70bp) TFBS sequences detected by CHIP-seq experients. Then the cfDNA fragment characteristics related to the TFBS was analyzed, includiong the motif information and fragment coverage depth around potential TFBSs and TFBSs distribution on the genome.
![image](https://github.com/QiT-lab/cfTFBS-DL/assets/104767835/e5be5604-46d5-444d-acb7-40d79cb1b679)
All the operation was performed in R language script, this is a good attempt to run deep learning models on the R language. 
First of all, keras and tensorflow R library should be installed successfully.
```
install.packages('tensorflow')
install.packages('keras')
```

Then, the python environment should be configured to connected to R language. The conda environment is recommended.
```
library(reticulate)
Sys.setenv(RETICULATE_PYTHON='your_python_path') 
repl_python()
```

