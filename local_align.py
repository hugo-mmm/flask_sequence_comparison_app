from Bio.Align import substitution_matrices as subs
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.Align import MultipleSeqAlignment

import sys
import os


# Initialize aligner globally

blosum = subs.load("BLOSUM62")
aligner = PairwiseAligner()
aligner.mode="global"
aligner.substitution_matrix = blosum
aligner.open_gap_score = -11
aligner.extend_gap_score = -1

def parse_sequence_data(data):
    lines = data.strip().split('\n')
    name_line = next(line for line in lines if line.startswith(">"))
    name = name_line[1:].strip()
    sequence = ''.join(''.join(line.split()[1:]) for line in lines[1:])
    sequence += lines[-1]
    return name, sequence


def print_similarity_matrix():
    print(str(blosum))    
    
def seq_similarity(seq1: str, seq2: str):
    # compute alignment
    alignments = aligner.align(seq1, seq2)
    
    # select the best alignment
    best_alignment = alignments[0]
    
    # compute aligned sequences
    aligned_seq1, aligned_seq2 = get_aligned_seqs(best_alignment)

    # compute positive matches (accounting for amino acid substitutions)
    positives, len_alignment, percent_positives = compute_positives(aligned_seq1, aligned_seq2)

    # compute amino acid substitutions
    substitutions = compute_substitutions(aligned_seq1, aligned_seq2)

    # compute identities
    identity, percent_ident = calc_identity(aligned_seq1, aligned_seq2)
    
    return positives, percent_positives, identity, percent_ident, substitutions, len_alignment, best_alignment


def compute_substitutions(seq1:str, seq2:str):
    substitutions = []
    for i, (a, b) in enumerate(zip(seq1, seq2), start=1):
        if a != b:
            substitutions.append(f"{a}{i}{b}")
    
    return ' '.join(substitutions) if substitutions else 'NONE'

def compute_positives(aligned_seq1: str, aligned_seq2: str):

    positives = 0
    
    for a,b in zip(aligned_seq1,aligned_seq2):
        sub_score = blosum.get((a,b),-4)
        iden = 1 if sub_score>0 else 0
        positives +=iden

    percent_positives = round((positives/len(aligned_seq1)*100),2)
    return positives, len(aligned_seq1), percent_positives

def calc_identity(aligned1:str, aligned2:str):

    ident = 0
    for a, b in zip(aligned1,aligned2):
        scr = 1 if a==b else 0
        ident+=scr

    return ident, round((ident/len(aligned1))*100,2)


def get_aligned_seqs(alignment: str):
    split_seq = str(alignment).split("\n")

    aligned_seq1 = ""
    aligned_seq2 = ""
    pivot_seq1 = 0
    pivot_seq2 = 2

    num_rows = int(len(split_seq) / 4)

    for i in range(num_rows):
        seq1_to_add = split_seq[pivot_seq1]
        seq2_to_add = split_seq[pivot_seq2]

        max_index = len(seq1_to_add) - 4 if i == num_rows - 1 else len(seq1_to_add)

        aligned_seq1 += seq1_to_add[20:max_index]
        aligned_seq2 += seq2_to_add[20:max_index]

        pivot_seq1 += 4
        pivot_seq2 += 4
    
    return aligned_seq1,aligned_seq2

def get_max_score(seq1:str, seq2:str) :
    score_seq1seq1 = aligner.align(seq1,seq1).score
    score_seq2seq2 = aligner.align(seq2,seq2).score

    return score_seq1seq1, score_seq2seq2


def main():
    # Example input sequence data
    #seq1 = "MKKVIAVLFVALIISPEAAPQGPEDSRAKPLAPAFDDLLAAMGALMVGFFLFVWGRQAP"
    #seq2 = "MKKVIAVLFVALIISPEAAPQGPEDSRAKPLAPAFDDLLAAMGALLVGFFLFVWGRQAP"

    # [DPA5-ITQB]-TYR
    seq1 = "MASLYPSPFSMAINASSPTGTYLSTSPCIPLLRNQNSNALKVAVWRVSCKATKDGDQKDREEESIINNGTFERRKILVGVGGLYGSLSNKALLSLAAPISPPDVTKCGPPVLPSGAKPTNCCPPISSKNIIDLTLPNNPQVKVRPAAHLVDSTYIQNYKEALRRMKALPLDDPRNFTQQANIHCAYCDGAYHQLGFPSLDFQVHNSWLFFPFHRWYLYFYEKILGSLIKDLDPNFAIPFWNWDSPNGMPIPTMYADPKSPLYDPLRNANHKPPKPVDLDYNGIEDQTPTQQQISTNLNTMYRQLVSSSKTSTLFFGSAYRAGEDSDPGGGVVENIPHGPVHVWTGDNTQPNFEDMGTLYSAARDPIFFSHHANVDRMWTIWKTLGGKRSDIKDPDWLESGFLFYDENKNLVRVKVKDCLDTRSLGYVYQDVDVPWLDSKPTPGIRRAVVKSFGAGAALAAETSKESTKWPVVLDSSVSVVVKRPRKGRNEREKKEEEEVLVIEGIEFERDLGVKFDVYINDEDDVHGGPTKTEFAGSFVSVAHKHKHKHSKMKTNLRLGITELLEDLGAEDDEHVVVTLVPRFGKGHVIVGGIKIEFHK"

    #Patent WO2022013407A1 SEQ08
    seq2 = "MKSAREDKVPWFPRKVSELDKCHHLVTKFDPDLDLDHPGFSDQVYRQRRKLIAEIAFQYKHGEPIPHVEYTAEEIATWKEVYVTLKGLYATHACREHLEGFQLLERYCGYREDSIPQLEDVSRFLKERTGFQLRPVAGLLSARDFLASLAFRVFQCTQYIRHASSPMHSPEPDCCHELLGHVPMLADRTFAQFSQDIGLASLGASDEEIEKLSTVYWFTVEFGLCKQNGELKAYGAGLLSSYGELLHSLSEEPEVRAFDPDTAAVQPYQDQTYQPVYFVSESFNDAKDKLRNYASRIQRPFSVKFDPYTLAIDVLDSPHTIQRSLEGVQDELHTLAHALSAIS"

    #>US 6242221 B1 SEQ26
    seq3 = "MKSAREDKVPWFPRKVSELDKCHHLVTKFDPDLDLDHPGFSDQVYRQRRKLIAEIAFQYKHGEPIPHVEYTAEEIATWKEVYVTLKGLYATHACREHLEGFQLLERYCGYREDSIPQLEDVSRFLKERTGFQLRPVAGLLSARDFLASLAFRVFQCTQYIRHASSPMHSPEPDCCHELLGHVPMLADRTFAQFSQDIGLASLGASDEEIEKLSTVYWFTVEFGLCKQNGELKAYGAGLLSSYGELLHSLSEEPEVRAFDPDTAAVQPYQDQTYQPVYFVSESFNDAKDKLRNYASRIQRPFSVKFDPYTLAIDVLDSPHTIQRSLEGVQDELHTLAHALSAIS"

    #>KR_1020090041874_2
    #seq1 = "MAAVAAAQKNREMFAIKKSYSIENGYPARRRSLVDDARFETLVVKQTKQSVLDEARLRANDSSCEPEIDDVSVQKEQPQTEDDADTGLTEEEVILQNAASENPEDDKVLQKAALILRLREGMSSLSRVLKTVEKHSGSVCHLETRPSKPDSTQLDVLIKVEMSRQNLLQLIKSMRQTSSLANVSLLGEDNISVKTPWFPRHASELDNCNHLMTKYEPDLDMNHPGFADKEYRARRKYIAEVAFSYKYGDAIPFIEYTETENKTWQAVFNTVKELMPKHACAEYNRVFAKLEEAGIFVPDRIPQLEEMCSFLKSQTGFTLRPAAGLLTSRDFLASLAFRIFQSTQYVRHVNSPYHTPEPDCIHELLGHMPLLADPSFAQFSQEIGLASLGASDEEIEKLSTVYWFTVEFGLCKEKGVVKAYGAGLLSAYGELLHALSDKPEHRPFEPASTAVQPYQDQEYQPVYFVAESFEDAKDKFRRWVSTMSRPFEVRFNPHTGRVEMLDSVDQLESLVHQLNTELLHLSNAINKLKKPVL"

    # Calculate similarity
    positives, percent_positives, identity, percent_identities, len_alignment, alignments = seq_similarity(seq1, seq3)

    # Print the similarity score    
    print("score best alignment:"+str(alignments[0].score))
    print("#positives:",positives)
    print("%positives:",percent_positives)
    print("#identities:", identity)
    print("%identities:",percent_identities)

    print("align_length:",len_alignment)
    print("best alignment:\n"+str(alignments[0]))
    
    #print("similarity matrix:")
    #print_similarity_matrix()

if __name__ == "__main__":
    main()
