import os

ws = r'H:\00code_ws\LaMata_new\MF_ws'
fn = ['h_MF4PEST.smp','h_MF4PESTdiff']
for e in fn:
    inputFile_fn = os.path.join(ws, e)
    if os.path.exists(inputFile_fn):
        fin = open(inputFile_fn, 'r')
        textin = fin.read()
        fin.close()
        textout = textin.replace('NaN', '0.0')
        fout = open(inputFile_fn, 'w')
        fout.write(textout)
        fout.close()
    else:
        self.ErrorExit(msg = "File [%s] doesn't exist, verify name and path!"%inputFile_fn)


