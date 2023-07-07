from flask import Flask, jsonify, request, render_template
from Bio.Align import PairwiseAligner
from Bio.Align import substitution_matrices
import sys
import os
import local_align as la

app = Flask(__name__)

seq1_name = "Sequence 1"
seq2_name = "Sequence 2"

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

@app.route('/align', methods=['POST'])
def align():
    data = request.get_json()
    #log_filename = 'error.log'
    #logging.basicConfig(filename=log_filename, level=logging.ERROR, format='%(asctime)s %(levelname)s: %(message)s')

    try:
        seq1_name_data = ">" + data['seq1_name'] + "\n" + data['seq1']
        seq2_name_data = ">" + data['seq2_name'] + "\n" + data['seq2']

        seq1_name, seq1 = parse_sequence_data(seq1_name_data)
        seq2_name, seq2 = parse_sequence_data(seq2_name_data)

        positives, percent_positives, identities, percent_identities, substitutions, len_alignment, alignment = la.seq_similarity(seq1, seq2)

        response = {
            'positives': positives,
            'percent_positives': '{:.2f}%'.format(percent_positives),
            'identities': identities,
            'percentage_identities': '{:.2f}%'.format(percent_identities),
            'substitutions': substitutions,
            'len_alignment': len_alignment,
            'alignment': str(alignment)
        }

        return jsonify(response), 200

    except KeyError as e:
        return jsonify({'error': 'Invalid request data. Missing key: {}'.format(e)}), 400
    except Exception as e:
        #logging.exception(e)
        return jsonify({'error': 'An error occurred during sequence comparison.'}), 500



if __name__ == '__main__':
    app.run(host='0.0.0.0', port=int(os.getenv('PORT', 80)), debug=True)

