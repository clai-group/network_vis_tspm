# Outputs from TSPM+
### corseq
8 Columns: patient_num, sequence, endphenx, durationBucket, count, value.var, startphen, startphen_dur\
*This table is used to create the freq table later.
### corrs_final
8 Columns: endphenx, startphen_dur, sequence, startphen, rho, p.value, p.adjust, rho.abs\
*This table is used to create the data_graph table later.

# Dataframes for plots
### Hierarchy
2 Columns: from, to

### Connect
3 Columns: from, to, value, freq\
*The value here is mean(rho).\
*The freq here is the frequency of each sequency.

### Vertices
3 Columns: name, value, group\
*The value here is the frequency of each ICD code.
