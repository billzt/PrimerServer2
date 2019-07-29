import json

from flask import Blueprint, request, current_app

from primerserver2.web.config import load

bp = Blueprint('db', __name__)
@bp.route('/dbselect')
def dbselect():
    web_config = load()
    dbs_in_select = {}
    for (template, meta) in web_config['templates'].items():
        (desc, group) = (meta['description'], meta['group'])
        if group not in dbs_in_select:
            dbs_in_select[group] = {}
        dbs_in_select[group][template] = desc
    return json.dumps(dbs_in_select, indent=4)
    
