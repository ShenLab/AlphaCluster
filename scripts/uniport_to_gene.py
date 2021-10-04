from os import walk
import os
import mygene

filenames = next(walk("../AlphaFold2_homo_sapiens"), (None, None, []))[2]
mg = mygene.MyGeneInfo()

for pdb_name in filenames:
    print(pdb_name.split('.')[0])
    symbol = mg.query("uniprot:"+pdb_name.split('.')[0], field=["symbol"])
    if not len(symbol["hits"]) == 0 and "symbol" in symbol["hits"][0]:
        print symbol["hits"][0]["symbol"]
        new_name = symbol["hits"][0]["symbol"]
        os.rename("../AlphaFold2_homo_sapiens/" + pdb_name, "../proteins/" + new_name+".pdb.gz")
