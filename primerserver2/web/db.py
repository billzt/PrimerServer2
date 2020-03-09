import json
import os

from flask import Blueprint, request, current_app, make_response

from primerserver2.web.config import load

bp = Blueprint('db', __name__)
@bp.route('/dbselect')
def dbselect():
    web_config = load()
    dbs_in_select = {}
    for (template, meta) in web_config['templates'].items():
        (desc, group, ids) = (meta['description'], meta['group'], meta['IDs'])
        if group not in dbs_in_select:
            dbs_in_select[group] = {}
        if template not in dbs_in_select[group]:
            dbs_in_select[group][template] = {}
        dbs_in_select[group][template]['desc'] = desc
        dbs_in_select[group][template]['IDs'] = ids
    return json.dumps(dbs_in_select, indent=4)

@bp.route('/dbinfo')
def dbinfo():
    web_config = load()
    return json.dumps(web_config['templates'], indent=4)
    
@bp.route('/dbdownload/<template>/')
def dbdownload(template):
    web_config = load()
    db_dir = web_config['templates_directory']
    result = os.popen(f'cut -f 1,2 {db_dir}/{template}.fai').read()
    response = make_response(result)
    response.headers["Content-Disposition"] = f"attachment; filename={template}.fai"
    response.headers["Content-type"] = "application/octet-stream"
    return response