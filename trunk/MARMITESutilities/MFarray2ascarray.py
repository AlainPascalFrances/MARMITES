import os

arrays_PEST_fn = [r'E:\00code_ws\LAMATA_new\MF_ws\00Sy_l1.ref', r'E:\00code_ws\LAMATA_new\MF_ws\00Sy_l2.ref', r'E:\00code_ws\LAMATA_new\MF_ws\00Ss_l1.ref', r'E:\00code_ws\LAMATA_new\MF_ws\00Ss_l2.ref', r'E:\00code_ws\LAMATA_new\MF_ws\00hk_l1.ref', r'E:\00code_ws\LAMATA_new\MF_ws\00hk_l2.ref']
ncol = 60

header_fn = r'E:\00code_ws\LAMATA_new\MF_ws\00header_asc.txt'
if os.path.exists(header_fn):
    header = open(header_fn, 'r')
    header_lst = []
    for line in header:
        header_lst.append(line)
    header.close()
else:
    raise SystemExit("File %s does not exist." % header_fn)   

print ''

for l, fin_fn in enumerate(arrays_PEST_fn):

    if os.path.exists(fin_fn):
        fin = open(fin_fn, 'r')
    else:
        raise SystemExit("File %s does not exist." % fin_fn) 
    fout_fn = '%s.asc' % fin_fn.split('.')[0]
    fout = open(fout_fn, 'w')

    for line in header_lst:
        fout.write(line)

    fout.write('\n')

    fin.readline()
    ncol_tmp = 0
    for line in fin:
        line_tmp = line.split()
        for e in line_tmp:
            fout.write(e + ' ')
            ncol_tmp += 1
        if ncol_tmp == ncol:
            ncol_tmp = 0
            fout.write('\n')

    fin.close()
    fout.close()
    print 'File %s done!' % fin_fn
    
print '\nDone!\n'
