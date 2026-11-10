import os as _os; _os.chdir(_os.path.dirname(_os.path.abspath(__file__)))  # paths below are relative to scripts/
import sys
import utils as utl

if __name__ == "__main__":
  inpath = sys.argv[1]
  outpath = sys.argv[2]
  print('inpath:', inpath)
  print('outpath:', outpath)
  utl.MakeDatabase(inpath, outpath)

