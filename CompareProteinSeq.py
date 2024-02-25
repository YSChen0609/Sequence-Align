"""
This is a module for experiments, can be adjusted for different needs.
"""

from string import Template
import requests
import re
import SeqAlign as sa

class CompareProteinSeq():
    
    def __init__(self, compare_set=("P02189", "P68082", "P02144")):
        self.UNIPROT_URL_TEMPLATE = Template("https://rest.uniprot.org/uniprotkb/${seq_id}.fasta")
        self.compare_set = compare_set
        self.proteinSeqs = {_id: (k,v) for _id,k,v in (self.getProteinSeq(_id) for _id in self.compare_set)}

    def getProteinSeq(self, seq_id):
        """
        Extract the Protein Sequence from UNIPROT, given the sequence id (name).
        """
        proteinseq_url = self.UNIPROT_URL_TEMPLATE.substitute(seq_id=seq_id)
        r = requests.get(proteinseq_url)
        seq_name = re.search(r"\|\w+\|(\w+)", r.text).group(1)
        seq = r.text.split("\n",1)[1].replace("\n", "")
        
        return seq_id, seq_name, seq
    
    def main(self):
        for seq_id in self.compare_set[1:]:
            s = sa.SeqAlign(self.proteinSeqs['P02189'],self.proteinSeqs[seq_id])
            
            print(f"{self.proteinSeqs['P02189'][0]}(P02189) vs. {self.proteinSeqs[seq_id][0]}({seq_id}):\n")
            print(self.proteinSeqs['P02189'][1])
            print(self.proteinSeqs[seq_id][1])
            print("before align")

            s.main()
            print("="*25)