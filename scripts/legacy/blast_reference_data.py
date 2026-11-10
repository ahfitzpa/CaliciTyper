import os as _os; _os.chdir(_os.path.dirname(_os.path.abspath(__file__)))  # paths below are relative to scripts/

import DataPreparation as dp

dp.DownloadFastaFromGenbak()

dp.annotate_refseq()
