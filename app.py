from flask import Flask, jsonify, request
from Bio.Align import PairwiseAligner
from Bio.Align import substitution_matrices
from flask import Flask, render_template
import sys
import os

app = Flask(__name__)

def parse_sequence_data(data):
    lines = data.strip().split('\n')
    name_line = next(line for line in lines if line.startswith(">"))
    name = name_line[1:].strip()
    sequence = ''.join(''.join(line.split()[1:]) for line in lines[1:])
    sequence += lines[-1]
    return name, sequence

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/seq_identity', methods=['POST'])
def seq_identity():
    data = request.get_json()
    try:
        seq1_name_data = ">" + data['seq1_name'] + "\n" + data['seq1']
        seq2_name_data = ">" + data['seq2_name'] + "\n" + data['seq2']
        
        seq1_name, seq1 = parse_sequence_data(seq1_name_data)
        seq2_name, seq2 = parse_sequence_data(seq2_name_data)

        blosum62 = substitution_matrices.load("BLOSUM62")
        
        aligner = PairwiseAligner()
        aligner.substitution_matrix = blosum62
        aligner.open_gap_score = -5
        aligner.extend_gap_score = -1
        
        alignments = aligner.align(seq1, seq2)
        best_alignment = alignments[0]
        split_seq = str(best_alignment).split("\n")
        
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
        
        if len(aligned_seq1) == 0:
            return jsonify({'identity percentage': 'N/A'}), 200

        identity = sum(a == b for a, b in zip(aligned_seq1, aligned_seq2)) / len(aligned_seq1) * 100
        identity = round(identity, 2)
        
        response = {
            'identity percentage': '{:.2f}%'.format(identity)  # Add the percentage symbol
        }
        
        return jsonify(response), 200
    except KeyError as e:
        return jsonify({'identity percentage': 'N/A', 'error': 'Invalid request data. Missing key: {}'.format(e)}), 400

@app.route('/seq_similarity', methods=['POST'])
def seq_similarity():
    data = request.get_json()
    try:
        seq1_name_data = ">" + data['seq1_name'] + "\n" + data['seq1']
        seq2_name_data = ">" + data['seq2_name'] + "\n" + data['seq2']
        
        seq1_name, seq1 = parse_sequence_data(seq1_name_data)
        seq2_name, seq2 = parse_sequence_data(seq2_name_data)
        
        blosum62 = substitution_matrices.load("BLOSUM62")
        
        aligner = PairwiseAligner()
        aligner.substitution_matrix = blosum62
        aligner.open_gap_score = -5
        aligner.extend_gap_score = -1
        
        alignments = aligner.align(seq1, seq2)
        best_alignment = alignments[0]
        split_seq = str(best_alignment).split("\n")
        
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
        
        if len(aligned_seq1) == 0:
            return jsonify({'similarity score': 'N/A'}), 200

        similarity = sum(blosum62.get((a, b), -4) for a, b in zip(aligned_seq1, aligned_seq2)) / len(aligned_seq1) * 100
        similarity = round(similarity, 2)
        
        response = {
            'similarity score': similarity
        }
        
        return jsonify(response), 200
    except KeyError as e:
        return jsonify({'similarity score': 'N/A', 'error': 'Invalid request data. Missing key: {}'.format(e)}), 400

@app.route('/seq_modifications', methods=['POST'])
def seq_modifications():
    data = request.get_json()
    try:
        seq1_name_data = ">" + data['seq1_name'] + "\n" + data['seq1']
        seq2_name_data = ">" + data['seq2_name'] + "\n" + data['seq2']
        
        seq1_name, seq1 = parse_sequence_data(seq1_name_data)
        seq2_name, seq2 = parse_sequence_data(seq2_name_data)
        
        substitutions = []
        for i, (a, b) in enumerate(zip(seq1, seq2), start=1):
            if a != b:
                substitutions.append(f"{a}{i}{b}")
        
        response = {
            'modifications': ' '.join(substitutions) if substitutions else 'NONE'  # Check if substitutions list is empty
        }
        
        return jsonify(response), 200
    except KeyError as e:
        return jsonify({'modifications': 'N/A', 'error': 'Invalid request data. Missing key: {}'.format(e)}), 400

@app.route('/seq_alignment', methods=['POST'])
def seq_alignment():
    data = request.get_json()
    try:
        seq1_name_data = ">" + data['seq1_name'] + "\n" + data['seq1']
        seq2_name_data = ">" + data['seq2_name'] + "\n" + data['seq2']
        
        seq1_name, seq1 = parse_sequence_data(seq1_name_data)
        seq2_name, seq2 = parse_sequence_data(seq2_name_data)
        
        blosum62 = substitution_matrices.load("BLOSUM62")
        
        aligner = PairwiseAligner()
        aligner.substitution_matrix = blosum62
        aligner.open_gap_score = -5
        aligner.extend_gap_score = -1
        
        alignments = aligner.align(seq1, seq2)
        best_alignment = alignments[0]
        
        alignment_str = str(best_alignment)
        
        response = {
            'alignment': alignment_str
        }
        
        return jsonify(response), 200
    except KeyError as e:
        return jsonify({'alignment': 'N/A', 'error': 'Invalid request data. Missing key: {}'.format(e)}), 400
    
if __name__ == '__main__':
    app.run(host='0.0.0.0', port=int(os.getenv('PORT', 80)))
