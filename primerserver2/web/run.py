import os
import json

from flask import Blueprint, request, current_app

from primerserver2.web.config import load
from primerserver2.core import make_sites, make_primers, design_primer, run_blast, sort_primers, output

web_config = load()
db_dir = os.path.join(os.path.dirname(__file__), '../templates/')

bp = Blueprint('run', __name__)
@bp.route('/run', methods=['POST'])
def run():
    ###################  Design primers ###################
    query_string = request.form['query']
    dbs = [db_dir+'/'+x for x in request.form['selected_dbs'].split(',')]
    if request.form['app-type']=='check':
        primers = make_primers.make_primers(query=query_string)
    else:
        if request.form['app-type']=='design':
            sites = make_sites.build(query=query_string, template_file=dbs[0], primer_type=request.form['region_type'], \
                primer_num_return=int(request.form['retain']), size_min=int(request.form['product_size_min']), \
                    size_max=int(request.form['product_size_max']))
        else:
            sites = make_sites.build(query=query_string, template_file=dbs[0], primer_type=request.form['region_type'], \
                primer_num_return=int(request.form['primer_num_return']), size_min=int(request.form['product_size_min']), \
                    size_max=int(request.form['product_size_max']))
        primers = design_primer.multiple(sites, cpu=web_config['cpu'], monitor=False)
    if len(primers)==0:
        return json.dumps({'error': 'Your inputs might be invalid; Check the \
        <a href="javascript:void(0)" data-toggle="modal" data-target="#input-help">manual</a> carefully \
        to ensure that your inputs are in correct formats'}, indent=4)

    ###################  Checking specificity  #############
    if request.form['app-type']!='design':
        primers = run_blast.run_blast_parallel(primers=primers, dbs=dbs, cpu=web_config['cpu'],\
            checking_size_max=int(request.form['checking_size_max']), checking_size_min=int(request.form['checking_size_min']), \
                report_amplicon_seq=bool(int(request.form['report_amplicon_seqs'])), Tm_diff=int(request.form['Tm_diff']), \
                    use_3_end=bool(int(request.form['use_3_end'])), monitor=False)
        primers = sort_primers.sort_rank(primers=primers, dbs=dbs, max_num_return=int(request.form['retain']))

    return json.dumps(primers, indent=4)
    
    # user_input = request.form['user_input']
    # population = request.form.getlist('population')
    