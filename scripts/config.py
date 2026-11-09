import os, shutil
_D = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # repository root
paths = {
    # annotation file for blast ref seq
    "annotation": _D+"/norotyper/data/database/refseq_annotation.csv",

    # file prefix for fastas/orf/annotation/amino acids
    "prefix":
    {
      'GI':[_D+"/norotyper/data/reference_sequences_RdRp/GI", _D+"/norotyper/data/reference_sequences_VP1/GI"],
      'GII':[_D+"/norotyper/data/reference_sequences_RdRp/GII", _D+"/norotyper/data/reference_sequences_VP1/GII"],
      'GIII':[_D+"/norotyper/data/reference_sequences_RdRp/GIII", _D+"/norotyper/data/reference_sequences_VP1/GIII"],
      'GIV':[_D+"/norotyper/data/reference_sequences_RdRp/GIV", _D+"/norotyper/data/reference_sequences_VP1/GIV"],
      'GV':[_D+"/norotyper/data/reference_sequences_RdRp/GV", _D+"/norotyper/data/reference_sequences_VP1/GV"],
      'GVI':[_D+"/norotyper/data/reference_sequences_RdRp/GVI", _D+"/norotyper/data/reference_sequences_VP1/GVI"],
      'GVII':[_D+"/norotyper/data/reference_sequences_RdRp/GVII", _D+"/norotyper/data/reference_sequences_VP1/GVII"],
      'GVIII':["", _D+"/norotyper/data/reference_sequences_VP1/GVIII"],
      'GNA1':[_D+"/norotyper/data/reference_sequences_RdRp/GNA1", _D+"/norotyper/data/reference_sequences_VP1/GNA1"],
      'GNA2':[_D+"/norotyper/data/reference_sequences_RdRp/GNA2", _D+"/norotyper/data/reference_sequences_VP1/GNA2"],
      'GIX':["", _D+"/norotyper/data/reference_sequences_VP1/GIX"],
      'GX':[_D+"/norotyper/data/reference_sequences_RdRp/GX", _D+"/norotyper/data/reference_sequences_VP1/GX"],
    },
    "aln":
      [_D+"/norotyper/data/reference_sequences_RdRp/all.aln", _D+"/norotyper/data/reference_sequences_VP1/all.aln"]
    ,
    # extracted from the published paper of noronet
    "taxonomy":
    {
      "rdrp":_D+'/validation/rdrp_taxonomy.csv',
      "vp1":_D+'/validation/vp1_taxonomy.csv'
    },
    # BLAST dataset
    "database": _D+"/norotyper/data/database/refseqdb_2",
    "exe":
    {
        "makeblastdb": os.environ.get("NOROTYPER_MAKEBLASTDB") or shutil.which("makeblastdb"),
        "blastn": os.environ.get("NOROTYPER_BLASTN") or shutil.which("blastn"),
        "clustal": os.environ.get("NOROTYPER_CLUSTALW") or shutil.which("clustalw2") or shutil.which("clustalw"),
    },
}
# reference sequences for Blast
numbers = {
    'GI':['NC_001959', 'NC_039897.1', 'NC_044853','NC_044854.1', 'NC_044856.1'],
    'GII':['NC_029646', 'NC_039475.1', 'NC_039476.1', 'NC_039477.1','NC_040876', 'NC_044045.1','NC_044046.1', 'NC_044932'], 
    'GIII':['NC_029645'],
    'GIV':['NC_029647', 'NC_044855.1'],
    'GV':['NC_008311', 'MW174170.1'],
    'GVI':['NC_044047','MW662289.1','MN908340.1','MW945229.1'],
    'GVII':['FJ692500', 'OL757872.1'],
    'GVIII':['AB985418.2'],
    'GIX':['OR050586.1','OR050585.1','OR050584.1'],
    'GX':['KJ790198.1','MF373609.1'],
    'Caliciviridae_Sapovirus_GI':['NC_006269'],
    'Caliciviridae_Sapovirus_GII':['NC_006554'],
    'Caliciviridae_Sapovirus_GIII':['NC_000940'],
    'Caliciviridae_Sapovirus_GIV':['DQ366346'],
    'Caliciviridae_Sapovirus_GV':['NC_027026'],
    'Caliciviridae_Sapovirus_GVI':['KJ508818'],
    'Caliciviridae_Sapovirus_GVII':['MF766258'],
    'Caliciviridae_Sapovirus_GVIII':['LC215895'],
    'Caliciviridae_Sapovirus_GX':['LC215896'],
    'Caliciviridae_Sapovirus_GXI':['LC215899'],
    'Caliciviridae_Sapovirus_GXII':['KX000385'],
    'Caliciviridae_Recovirus_unclassified':['EU391643'],
    'Caliciviridae_Chicken_calicivirus_unclassified':['HQ010042'],
    'Caliciviridae_Vesivirus_Feline_calicivirus':['NC_001481'],
    'Caliciviridae_Nacovirus_Turkey_calicivirus':['JQ347522'],
    'Caliciviridae_Lagovirus_Rabbit_hemorrhagic_disease_virus':['NC_001543'],
    'Caliciviridae_Lagovirus_unclassified':['NC_011704'],
    'Caliciviridae_unclassified_Bovine_calicivirus':['KT119483'],
    'Caliciviridae_Vesivirus_Vesicular_exanthema_of_swine_virus':['NC_002551'],
    'Caliciviridae_Lagovirus_European_brown_hare_syndrome_virus':['NC_002615'],
    'Caliciviridae_Vesivirus_Canine_calicivirus':['NC_004542'],
    'Caliciviridae_Nebovirus_unclassified':['NC_006875'],
    'Caliciviridae_Vesivirus_Rabbit_vesivirus':['NC_008580'],
    'Caliciviridae_Vesivirus_Steller_sea_lion_vesivirus':['NC_011050'],
    'Caliciviridae_Valovirus_unclassified':['NC_012699'],
    'Caliciviridae_Sapovirus_Bat_sapovirus':['NC_017936'],
    'Caliciviridae_Salmonid_calicivirus_unclassified':['NC_024031'],
    'Caliciviridae_Nacovirus_Goose_calicivirus':['NC_024078'],
    'Caliciviridae_Vesivirus_Ferret_badger_vesivirus':['NC_027122'],
    'Caliciviridae_Vesivirus_San_Miguel_sea_lion_virus':['U15301']
}

# All the accession number of norovirus for BLAST analysis
# RIVM has only one refseq for each genogroup, in this project, more refseqs are included.
noronumbers = {
    'GI':['NC_001959', 'NC_039897.1', 'NC_044853','NC_044854.1', 'NC_044856.1'],
    'GII':['NC_029646', 'NC_039475.1', 'NC_039476.1', 'NC_039477.1','NC_040876', 'NC_044045.1','NC_044046.1', 'NC_044932'], 
    'GIII':['NC_029645'],
    'GIV':['NC_029647', 'NC_044855.1'],
    'GV':['NC_008311', 'MW174170.1'],
    'GVI':['NC_044047','MW662289.1','MN908340.1','MW945229.1'],
    'GVII':['FJ692500', 'OL757872.1'],
    'GVIII':['AB985418.2'],
    'GIX':['OR050586.1','OR050585.1','OR050584.1'],
    'GX':['KJ790198.1','MF373609.1']
}
