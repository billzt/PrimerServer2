import os
import json
import time

from flask import Blueprint, request, current_app, Response

from primerserver2.web.config import load
from primerserver2.core import make_sites, make_primers, design_primer, run_blast, sort_primers, output, global_var, multiplex

# for running progress
global_var.init()
def progress():
    while True:
        yield "data:" + json.dumps({'complete_count': global_var.complete_count, \
            'all_tasks_num': global_var.all_tasks_num, 'current_task': global_var.current_task}) + '\n\n'
        time.sleep(1)
        # if global_var.complete_count>0 and global_var.complete_count==global_var.all_tasks_num:
        #     break
        if global_var.current_task=='finish':
            break
    yield "data:" + json.dumps({'complete_count': global_var.complete_count, \
        'all_tasks_num': global_var.all_tasks_num, 'current_task': global_var.current_task}) + '\n\n'

bp = Blueprint('run', __name__)
@bp.route('/run', methods=['POST'])
def run():
    
    ###################  init #############################
    web_config = load()
    global_var.init()
    db_dir = web_config['templates_directory']

    ###################  Design primers ###################
    query_string = request.form['query']
    dbs = [db_dir+'/'+x for x in request.form['selected_dbs'].split(',')]
    if request.form['app-type']=='check':
        primers = make_primers.make_primers(query=query_string)
        if 'error' in primers:
            return json.dumps({'error': '<p>'+primers['error']+'</p>'+'<p>Your inputs might be invalid; Check the \
                <a href="javascript:void(0)" data-toggle="modal" data-target="#input-help">manual</a> carefully \
                to ensure that your inputs are in correct formats</p>'}, indent=4)
    else:
        if request.form['app-type']=='design':
            primer_num_return = request.form['retain']
        else:
            primer_num_return = request.form['primer_num_return']
        if make_sites.judge_input_type(query_string)=='pos':
            sites = make_sites.build_by_pos(query=query_string, template_file=dbs[0], primer_type=request.form['region_type'], \
                primer_num_return=int(primer_num_return), size_min=int(request.form['product_size_min']), \
                    size_max=int(request.form['product_size_max']))
        else:
            sites = make_sites.build_by_seq(query=query_string.replace('\r\n', '\n'), primer_type=request.form['region_type'], \
                primer_num_return=int(primer_num_return), size_min=int(request.form['product_size_min']), \
                    size_max=int(request.form['product_size_max']))
        if 'error' in sites:
            return json.dumps({'error': '<p>'+sites['error']+'</p>'+'<p>Your inputs might be invalid; Check the \
                <a href="javascript:void(0)" data-toggle="modal" data-target="#input-help">manual</a> carefully \
                to ensure that your inputs are in correct formats</p>'}, indent=4)
        primers = design_primer.multiple(sites, cpu=web_config['cpu'], monitor=False)

    ###################  Checking specificity  #############
    if request.form['app-type']!='design':
        primers = run_blast.run_blast_parallel(primers=primers, dbs=dbs, cpu=web_config['cpu'],\
            checking_size_max=int(request.form['checking_size_max']), checking_size_min=int(request.form['checking_size_min']), \
                report_amplicon_seq=bool(int(request.form['report_amplicon_seqs'])), Tm_diff=int(request.form['Tm_diff']), \
                    use_3_end=bool(int(request.form['use_3_end'])), monitor=False)
        primers = sort_primers.sort_rank(primers=primers, dbs=dbs, max_num_return=int(request.form['retain']))

    ###################  Checking multiplex  ###############
    dimers = {}
    if int(request.form['multiplex'])==1:
        dimers = multiplex.extract_fake_pair(primers, Tm_diff=int(request.form['Tm_diff']), cpu=web_config['cpu'], monitor=False)

    ###################  Output    #########################
    global_var.current_task = 'finish'
    return json.dumps({'meta':{'mode':request.form['app-type'], 'dbs':dbs, 'region_type': request.form['region_type'], \
         'check_multiplex': bool(int(request.form['multiplex']))}, 'primers':primers, 'dimers':dimers}, indent=4)
    

@bp.route('/monitor')
def monitor():
    time.sleep(1)
    return Response(progress(), mimetype= 'text/event-stream')

@bp.route('/stop')
def stop():
    global_var.stop_run = True
    return ""