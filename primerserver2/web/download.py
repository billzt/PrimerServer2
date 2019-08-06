import json

from flask import Blueprint, request

from primerserver2.core import output

bp = Blueprint('download', __name__)
@bp.route('/download_tsv', methods=['POST'])
def download_tsv():
    json_data = json.loads(request.form['json'])
    dbs = request.form['dbs'].split(',')
    print_lines = output.tsv(json_data['primers'], dbs)
    return print_lines
