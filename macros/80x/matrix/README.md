# matrix



Run testAnalysis_test_new.C to split some of the histograms, the resulting histograms are stored in input

run matrix.C which takes root files from input to produce the covariance matrix, root files and pdfs, root files go to output_root and the pdfs go to plots

run convert.C which takes the root files from output_root directly, outputs correlation matrices(PDF only) to plots 
