import rpy2.robjects as robjects

r_code_file_name = 'LDA_label.R'
robjects.r.source(r_code_file_name, encoding="utf-8")