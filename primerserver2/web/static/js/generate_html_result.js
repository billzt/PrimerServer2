function generate_html_result(selected_dbs, db_name_change, data, visualize_mode){
    $('#primers-result').html('');
    var result_data = data; 
    var site_rank = 0;
    for (site_id in result_data) {
        site_rank += 1;

        // site
        var raw_html_site = $('#primers-result-template-site').html();
        $('#primers-result-template-site a.collapsed').attr('href', '#site-'+site_rank)
            .html('Site '+site_rank+' <span class="caret"></span>');
        $('#primers-result-template-site .site-detail').html('Site: '+site_id);
        if (visualize_mode!='check') {
            let i = site_id.replace(/-PARTA$/, '').replace(/-PARTB$/, '').split('-');
            site_seq = i.slice(0,-2).join('-');
            site_pos = i[i.length-2];
            site_len = i[i.length-1];
            $('#primers-result-template-site .site-detail').attr('data-seq', site_seq).attr('data-pos', site_pos)
                .attr('data-length', site_len);
        }
        else {
            site_seq = site_id;
        }
        var primer_num = 0;
        if ('PRIMER_PAIR_NUM_RETURNED_FINAL' in result_data[site_id]) {
            primer_num = result_data[site_id]['PRIMER_PAIR_NUM_RETURNED_FINAL']
        }
        else {
            primer_num = result_data[site_id]['PRIMER_PAIR_NUM_RETURNED']
        }
        $('#primers-result-template-site .site-primer-num').html(primer_num);
        var main_db = selected_dbs.split(',')[0];
        if (main_db in result_data[site_id]) {
            var least_amplicon_num_in_main_db = result_data[site_id][main_db]['PRIMER_PAIR_'
                +result_data[site_id]['PRIMER_PAIR_AMPLICON_NUM_RANK_0']+'_AMPLICONS']
                    .filter(x=> x['plus']['sseqid']==site_seq || x['isoform']!=true).length;
            if (least_amplicon_num_in_main_db==1) {
                $('#primers-result-template-site .row-site-info')
                    .append('<div class="col-md-1"><span class="glyphicon glyphicon-ok"></span></div>');
            }
        }
        $('#primers-result-template-site .panel-site').attr('id', 'site-'+site_rank);
        if (site_rank==1) {
            $('#site-1').addClass('in');
        }
        $('#primers-result').append($('#primers-result-template-site').html());
        $('#primers-result-template-site').html(raw_html_site);

        if (typeof(site_pos) != 'undefined') {
            retrieve_start = site_pos - result_data[site_id]['SEQUENCE_RELATIVE_TARGET_START']
        }
        else {
            retrieve_start = 0
        }
        
        // primer
        for (var i=0; i<primer_num; i++) {
            var raw_html_primer = $('#primers-result-template-primer').html();
            var primer_rank = i+1;
            var raw_rank = 0;
            if ('PRIMER_PAIR_NUM_RETURNED_FINAL' in result_data[site_id]) {
                raw_rank = result_data[site_id]['PRIMER_PAIR_AMPLICON_NUM_RANK_'+i];
            }
            else {
                raw_rank = i;
            }

            $('#primers-result-template-primer h4').attr('id', 'Site'+site_rank+'-Primer'+primer_rank)
                .html('Primer '+primer_rank);
            
            var p_start = result_data[site_id]['PRIMER_LEFT_'+raw_rank][0];
            p_start = p_start==-1 ? '' : p_start+retrieve_start;
            var p_len = result_data[site_id]['PRIMER_LEFT_'+raw_rank][1];
            var p_end = p_start=='' ? '' : p_start+p_len-1;
            $('#primers-result-template-primer .primer-seq-detail')
                .append('<tr><th>Left Primer</th>'
                    +'<td><span class="monospace-style">'+result_data[site_id]['PRIMER_LEFT_'+raw_rank+'_SEQUENCE']+'</span></td>'
                    +'<td>'+p_len+'</td>'
                    +'<td>'+site_seq+':<span class="primer-left-region">'+p_start+'-'+p_end+'</span></td>'
                    +'<td>'+result_data[site_id]['PRIMER_LEFT_'+raw_rank+'_TM'].toFixed(1)+'</td>'
                    +'<td>'+result_data[site_id]['PRIMER_LEFT_'+raw_rank+'_GC_PERCENT'].toFixed(1)+'</td>'
                    +'</tr>')

            p_end = result_data[site_id]['PRIMER_RIGHT_'+raw_rank][0];
            p_end = p_end==-1 ? '' : p_end+retrieve_start;
            p_len = result_data[site_id]['PRIMER_RIGHT_'+raw_rank][1];
            p_start = p_end=='' ? '' : p_end-p_len+1;
            $('#primers-result-template-primer .primer-seq-detail')
                .append('<tr><th>Right Primer</th>'
                    +'<td><span class="monospace-style">'+result_data[site_id]['PRIMER_RIGHT_'+raw_rank+'_SEQUENCE']+'</span></td>'
                    +'<td>'+p_len+'</td>'
                    +'<td>'+site_seq+':<span class="primer-right-region">'+p_start+'-'+p_end+'</span></td>'
                    +'<td>'+result_data[site_id]['PRIMER_RIGHT_'+raw_rank+'_TM'].toFixed(1)+'</td>'
                    +'<td>'+result_data[site_id]['PRIMER_RIGHT_'+raw_rank+'_GC_PERCENT'].toFixed(1)+'</td>'
                    +'</tr>')
            
            if ('PRIMER_INTERNAL_'+raw_rank+'_SEQUENCE' in result_data[site_id] == true) {
                p_start = result_data[site_id]['PRIMER_INTERNAL_'+raw_rank][0]+retrieve_start;
                p_len = result_data[site_id]['PRIMER_INTERNAL_'+raw_rank][1];
                p_end = p_start+p_len-1;
                $('#primers-result-template-primer .primer-seq-detail')
                .append('<tr><th>Internal Oligo</th>'
                    +'<td><span class="monospace-style">'+result_data[site_id]['PRIMER_INTERNAL_'+raw_rank+'_SEQUENCE']+'</span></td>'
                    +'<td>'+p_len+'</td>'
                    +'<td class="primer-internal-region">'+p_start+'-'+p_end+'</td>'
                    +'<td>'+result_data[site_id]['PRIMER_INTERNAL_'+raw_rank+'_TM'].toFixed(1)+'</td>'
                    +'<td>'+result_data[site_id]['PRIMER_INTERNAL_'+raw_rank+'_GC_PERCENT'].toFixed(1)+'</td>'
                    +'</tr>')
            }

            var product_size = result_data[site_id]['PRIMER_PAIR_'+raw_rank+'_PRODUCT_SIZE'];
            product_size = product_size==-1? '-' : product_size;
            $('#primers-result-template-primer .primer-seq-detail')
                .append('<tr><th>Product Size</th><td colspan="5">'+product_size+' bp</td></tr>');
            var penalty = result_data[site_id]['PRIMER_PAIR_'+raw_rank+'_PENALTY'].toFixed(2);
            penalty = penalty==0? '-' : penalty;
            $('#primers-result-template-primer .primer-seq-detail')
                .append('<tr><th>Penalty</th><td colspan="5">'+penalty+'</td></tr>');


            // dbs
            var db_rank = 0;
            for (db of selected_dbs.split(',')) {
                if (db in result_data[site_id] == false) {
                    // $('#primers-result-template-primer .amplicons_table').addClass('hidden');
                    continue;
                }
                db_desc = db in db_name_change ? db_name_change[db] : db
                if (db_rank==0) {
                    $('#primers-result-template-primer .database-list').append('<th>Database: '
                        +db_desc+' <span class="glyphicon glyphicon-star"></span></th>')
                }
                else {
                    $('#primers-result-template-primer .database-list').append('<th>Database: '
                        +db_desc+'</th>')
                }
                // amplicons
                if ('PRIMER_PAIR_'+raw_rank+'_AMPLICONS' in result_data[site_id][db]) {
                    amplicons = result_data[site_id][db]['PRIMER_PAIR_'+raw_rank+'_AMPLICONS'];
                } 
                else {
                    amplicons = [];
                }
                // amplicon num
                if (amplicons.length>0 && 'isoform' in amplicons[0]) {
                    amplicon_num = amplicons.filter(x=> x['plus']['sseqid']==site_seq || x['isoform']!=true).length;
                }
                else {
                    amplicon_num = amplicons.length;
                }
                $('#primers-result-template-primer .amplicons_number').append('<td class="hit-num" data-hit="'
                    +amplicon_num+'">Amplicon Number: '+amplicon_num+'</td>');
                if (amplicon_num==1 && db==main_db) {
                    $('#primers-result-template-primer .list-group-item-primer').addClass('list-group-item-success')
                }
                // print amplicons
                $('#primers-result-template-primer .amplicons_region')
                    .append('<td><ul class="list-group list-group-amplicons-'+db_rank+'"></ul></td>');
                var output_amplicon_num = 0;
                for (var j in amplicons) {
                    $('#primers-result-template-primer .list-group-amplicons-'+db_rank)
                        .append('<li class="list-group-item">'+amplicons[j]['region']+', '+amplicons[j]['product_size']+' bp'
                            +' <a href="javascript:void(0)" data-toggle="modal" data-target="#amplicon-details-modal" data-amplicon='
                            +JSON.stringify(amplicons[j])
                            +'>'
                            +'<span class="glyphicon glyphicon-hand-right"></span></a>');
                    amplicon_li = $('#primers-result-template-primer .list-group-amplicons-'+db_rank+' li:eq('+j+')')
                    if (amplicon_num==1) {
                        amplicon_li.addClass('list-group-item-success')
                    }
                    if (amplicons[j]['isoform']==true) {
                        amplicon_li.append(' <i class="fa fa-random"></i>');
                        if (amplicons[j]['plus']['sseqid']==site_seq) {
                            output_amplicon_num++;
                        }
                    }
                    else {
                        output_amplicon_num++;
                    }
                    if (output_amplicon_num==3) {
                        if (amplicon_num==1) {
                            $('#primers-result-template-primer .list-group-amplicons-'+db_rank)
                                .append('<li class="list-group-item list-group-item-success">...</li>');
                        }
                        else {
                            $('#primers-result-template-primer .list-group-amplicons-'+db_rank)
                                .append('<li class="list-group-item">...</li>');
                        }
                        break;
                    }
                }
                db_rank++;
            }

            $('#site-'+site_rank).find('.list-group-primers').append($('#primers-result-template-primer').html());
            $('#primers-result-template-primer').html(raw_html_primer);
        }

        // no primer
        if (primer_num==0) {
            var el = $('#site-'+site_rank).find('.list-group-primers')
            el.append($('#primers-result-template-error').html());
            el.find('.alert-danger').append('<h5>No primers</h5>');
            el.find('.alert-danger').append('<p>PRIMER_LEFT_EXPLAIN: '+result_data[site_id]['PRIMER_LEFT_EXPLAIN']+'</p>');
            el.find('.alert-danger').append('<p>PRIMER_RIGHT_EXPLAIN: '+result_data[site_id]['PRIMER_RIGHT_EXPLAIN']+'</p>');
            if('PRIMER_INTERNAL_EXPLAIN' in result_data[site_id]) {
                el.find('.alert-danger').append('<p>PRIMER_INTERNAL_EXPLAIN: '+result_data[site_id]['PRIMER_INTERNAL_EXPLAIN']+'</p>');
            }
            el.find('.alert-danger').append('<p>PRIMER_PAIR_EXPLAIN: '+result_data[site_id]['PRIMER_PAIR_EXPLAIN']+'</p>');
        }
    }
}

function generate_dimer_result(data){
    $('#dimers-result tbody').html('');
    $('#dimers-result').next('.alert-danger,.alert-success').remove();
    if ('error' in data) {
        $('#dimers-result').after($('#primers-result-template-error').html());
        $('#dimers-result').next('.alert-danger').append('<p>ERROR: '+data['error']+'</p>');
        return
    }
    if (jQuery.isEmptyObject(data)==true) {
        $('#dimers-result').after($('#primers-result-template-error').html());
        $('#dimers-result').next('.alert-danger').append('<p>No dimers detected. All primers can be used for multiplex PCR &hearts;</p>')
            .removeClass('alert-danger').addClass('alert-success');
        return
    }
    for (site_1 in data) {
        for (primer_1 in data[site_1]) {
            for (site_2 in data[site_1][primer_1]) {
                for (primer_2 in data[site_1][primer_1][site_2]) {
                    dimer = data[site_1][primer_1][site_2][primer_2]
                    $('#dimers-result tbody').append('<tr>'
                        +'<td>'+site_1+'</td>'
                        +'<td>'+primer_1+'</td>'
                        +'<td><span class="monospace-style">'+dimer['seq_1']+'</span></td>'
                        +'<td>'+site_2+'</td>'
                        +'<td>'+primer_2+'</td>'
                        +'<td><span class="monospace-style">'+dimer['seq_2']+'</span></td>'
                        +'<td>'+dimer['Tm']+'</td>'
                        +'<td><pre>'+dimer['ascii_structure']+'</pre></td>'
                        +'</tr>')
                }
            }
        }
    }
}