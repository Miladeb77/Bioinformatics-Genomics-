# retreiving mRNA sequence 
from Bio import SeqIO 
for seq_record in SeqIO.parse("insuline.fasta","fasta"): 
    print(seq_record.id)
    print(seq_record.description)
    print(str(seq_record.seq))
    print(len(seq_record))
    
#alignment 
from Bio.Align.Applications import ClustalwCommandline
from Bio import AlignIO
align = AlignIO.read("opuntia.aln", "clustal")
print(align)

#phylogenic tree
from Bio import Phylo
tree = Phylo.read("opuntia.dnd", "newick")
Phylo.draw_ascii(tree)

#BLAST with NCBI 
from Bio.Blast import NCBIWWW
fasta_string= open('myseq.fa').read() # myseq.fa is the patient's sample 
result_handle= NCBIWWW.qblast('blastn','nt', fasta_string)

#using XML package
from Bio.Blast import NCBIXML
blast_record= NCBIXML.read(result_handle)

#parse the output 
len(blast_record.alignments)
#filtering the matches based on the e value 
E_VALUE_THRESH= 0.01
for alignment in blast_record.alignments:
    for hsp in alignment.hsps: # hsps are the subrecords of the alignments 
        if hsp.expect < E_VALUE_THRESH: 
            print('****Alignment****')
            print('sequence:', alignment.title)
            print('length:', alignment.length)
            print('e value:', hsp.expect)
            print(hsp.query)
            print(hsp.match)
            print(hsp.sbjct)
