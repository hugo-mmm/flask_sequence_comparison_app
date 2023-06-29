from flask import Flask, jsonify, request, render_template
from Bio.Align import PairwiseAligner
from Bio.Align import substitution_matrices
import sys
import os

app = Flask(__name__)

# Initialize aligner globally
aligner = PairwiseAligner()
blosum62 = substitution_matrices.load("BLOSUM62")
aligner.substitution_matrix = blosum62
aligner.open_gap_score = -5
aligner.extend_gap_score = -1

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
            return jsonify({'identity percentage': 'N/A', 'error': 'No alignment found.'}), 200

        identity = sum(a == b for a, b in zip(aligned_seq1, aligned_seq2)) / len(aligned_seq1) * 100
        identity = round(identity, 2)

        response = {
            'identity percentage': '{:.2f}%'.format(identity)  # Add the percentage symbol
        }

        return jsonify(response), 200
    except KeyError as e:
        return jsonify({'identity percentage': 'N/A', 'error': 'Invalid request data. Missing key: {}'.format(e)}), 400
    except Exception as e:
        return jsonify({'identity percentage': 'N/A', 'error': str(e)}), 500

@app.route('/seq_similarity', methods=['POST'])
def seq_similarity():
    data = request.get_json()
    try:
        seq1_name_data = ">" + data['seq1_name'] + "\n" + data['seq1']
        seq2_name_data = ">" + data['seq2_name'] + "\n" + data['seq2']

        seq1_name, seq1 = parse_sequence_data(seq1_name_data)
        seq2_name, seq2 = parse_sequence_data(seq2_name_data)

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



        if len(aligned_seq1) > 0:
            similarity_seq1 = sum(blosum62.get((a, b), -4) for a, b in zip(aligned_seq1, aligned_seq1))
        else:
            similarity_seq1 = 0

        if len(aligned_seq2) > 0:
            similarity_seq2 = sum(blosum62.get((a, b), -4) for a, b in zip(aligned_seq2, aligned_seq2))
        else:
            similarity_seq2 = 0

        min_similarity = min(similarity_seq1, similarity_seq2)


        if len(aligned_seq1) > 0:
            similarity = (sum(blosum62.get((a, b), -4) for a, b in zip(aligned_seq1, aligned_seq2))) /  min_similarity * 100
            similarity = round(similarity, 2)
        else:
            similarity = 0

        
        response = {
            'similarity score': '{:.2f}%'.format(similarity)  # Add the percentage symbol
        }

        return jsonify(response), 200
    except ValueError as e:
        return jsonify({'error': str(e)}), 400
    except KeyError as e:
        return jsonify({'error': 'Invalid request data. Missing key: {}'.format(e)}), 400
    except Exception as e:
        return jsonify({'error': str(e)}), 500

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
    except Exception as e:
        return jsonify({'modifications': 'N/A', 'error': str(e)}), 500

@app.route('/seq_alignment', methods=['POST'])
def seq_alignment():
    data = request.get_json()
    try:
        seq1_name_data = ">" + data['seq1_name'] + "\n" + data['seq1']
        seq2_name_data = ">" + data['seq2_name'] + "\n" + data['seq2']

        seq1_name, seq1 = parse_sequence_data(seq1_name_data)
        seq2_name, seq2 = parse_sequence_data(seq2_name_data)

        alignments = aligner.align(seq1, seq2)
        best_alignment = alignments[0]

        alignment_str = str(best_alignment)

        response = {
            'alignment': alignment_str
        }

        return jsonify(response), 200
    except ValueError as e:
        return jsonify({'error': str(e)}), 400
    except KeyError as e:
        return jsonify({'error': 'Invalid request data. Missing key: {}'.format(e)}), 400
    except Exception as e:
        return jsonify({'error': str(e)}), 500

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=int(os.getenv('PORT', 80)))
