from flask import Flask, jsonify, request, render_template
from Bio.Align import PairwiseAligner
from Bio.Align import substitution_matrices
from Bio import Align
import sys
import os
import local_align as la

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

@app.route('/maia_align', methods=['POST'])
def maia_align():
    data = request.get_json()
    try:
        #seq1_name_data = ">" + data['seq1_name'] + "\n" + data['seq1']
        #seq2_name_data = ">" + data['seq2_name'] + "\n" + data['seq2']

        seq1_name, seq1 = parse_sequence_data(seq1_name_data)
        seq2_name, seq2 = parse_sequence_data(seq2_name_data)

        positives, percent_positives, identities, percent_identities, substitutions, len_alignment, alignments = la.seq_similarity(seq1, seq2)

        response = {
            'positives': positives,
            'percent_positives': '{:.2f}%'.format(percent_positives),
            'identities': identities,
            'percentage_identities': '{:.2f}%'.format(percent_identities),
            'substitutions': substitutions,
            'length_alignment': len_alignment,
            'alignment': alignments[0]
        }

        return jsonify(response), 200

    except KeyError as e:
        return jsonify({'identity percentage': 'N/A', 'error': 'Invalid request data. Missing key: {}'.format(e)}), 400
    except Exception as e:
        return jsonify({'identity percentage': 'N/A', 'error': str(e)}), 500

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=int(os.getenv('PORT', 80)))
