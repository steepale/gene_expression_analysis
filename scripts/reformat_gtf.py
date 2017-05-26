# import sys

# infile=sys.argv[1]

# outfile=open(sys.argv[2], 'w')


# for line in open(infile):
# 	col = line.split('\t')
# 	transcript_id = col[8].split(';')[0]
# 	gene_id = col[8].split(';')[2].replace('gene_name', 'gene_id')
# 	outfile.write(col[0] + '\t' + col[1] + '\t' + col[2] + '\t' + col[3] + '\t' + col[4] + '\t' + col[5] + '\t' + col[6] + '\t' + col[7] + '\t' + transcript_id + '; ' + gene_id + '\n')
# print('Finished')
