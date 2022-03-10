from os import walk
import os
import mygene

filenames = next(walk("../coevol"), (None, None, []))[2]
mg = mygene.MyGeneInfo()

for transcript_name in filenames:
    print(transcript_name.split('.')[0])
    if not transcript_name.split('.')[0].startswith('ENST'):
       continue
   
    symbol = mg.query("ensembl.transcript:"+transcript_name.split('.')[0], field=["symbol"])
    if not len(symbol["hits"]) == 0 and "symbol" in symbol["hits"][0]:
        print(symbol["hits"][0]["symbol"])
        new_name = symbol["hits"][0]["symbol"]
        os.rename("../coevol/" + transcript_name, "../coevol/" + new_name +".coevol.txt")
