function tooltip_init() {
    // Tooltip for bootstrap
    $('[data-toggle="tooltip"]').tooltip({html: true});
    $('[data-toggle="popover"]').popover({html: true});
    $('[data-toggle="table"]').on('post-body.bs.table', function () {
        $('[data-toggle="tooltip"]').tooltip({
            container: 'body'
        });
    });
}

function menu_init(data) {
    /* init the bootstrap select menu: $('[name="templates"]')
        return an object db_name_change: template => description
    */
    var db_name_change = new Object;
    var group_data = JSON.parse(data);
    for (group in group_data) {
        $('[name="templates"]').append('<optgroup label="'+group+'">');
        for (template in group_data[group]) {
            $('optgroup[label="'+group+'"]').append('<option data-subtext="seq IDs:'+group_data[group][template]['IDs']
                +'" value="'+template+'">'+group_data[group][template]['desc']+'</option>');
            db_name_change[template] = group_data[group][template]['desc'];
            $('#help-modal-ID-list').append('<a href="'+$SCRIPT_ROOT+'/dbdownload/'+template+'/" class="list-group-item">'
            +group_data[group][template]['desc']+' <span class="glyphicon glyphicon-download"></span></a>');
        }
    }
    return db_name_change
}

function qPCR_filedset(config_data, db, mode, region_type) {
    // junction
    if (mode=='full' || mode=='design') {
        if (db in config_data && config_data[db]['junction']==true && region_type=='SEQUENCE_INCLUDED_REGION') {
            $('#junction_1').parents('.form-group').removeClass('hidden');
            $('#junction_1').prop('checked', true);
        }
        else {
            $('#junction_1').parents('.form-group').addClass('hidden');
            $('#junction_0').prop('checked', true);
        }
    }
    else {
        $('#junction_1').parents('.form-group').addClass('hidden');
        $('#junction_0').prop('checked', true);
    }

    // isoform
    if (mode=='check' || (mode=='full' && region_type=='SEQUENCE_INCLUDED_REGION')) {
        if (db in config_data && config_data[db]['isoform']==true) {
            $('#isoform_1').parents('.form-group').removeClass('hidden');
            $('#isoform_1').prop('checked', true);
        }
        else {
            $('#isoform_1').parents('.form-group').addClass('hidden');
            $('#isoform_0').prop('checked', true);
        }
    }
    else {
        $('#isoform_1').parents('.form-group').addClass('hidden');
        $('#isoform_0').prop('checked', true);
    }
}

function switch_parameter_fieldset(mode) {
    if (mode=='visulize') {
        $('#form-primer,#result').addClass('hidden');
        $('#block-visulization').removeClass('hidden');
    }
    else {
        $('#block-visulization').addClass('hidden');
        $('#form-primer').removeClass('hidden');
        if ($('#primers-result').html()!='') {
            $('#result').removeClass('hidden');
            $('#footer-button').removeClass('hidden');
        }
        switch(mode) {
            case 'full': {
                $('#parameter-design').removeClass('hidden');
                $('#parameter-check').removeClass('hidden');
                $('#textarea-title').html('Input: Design primers and check specificity');
                break;
            }
            case 'design': {
                $('#parameter-design').removeClass('hidden');
                $('#parameter-check').addClass('hidden');
                $('#textarea-title').html('Input: Only design primers');
                break;
            }
            case 'check': {
                $('#parameter-design').addClass('hidden');
                $('#parameter-check').removeClass('hidden');
                $('#textarea-title').html('Input: Only check specificity');
                break;
            }
        }
    }
}

function highlight_changed_field(originalValFor) {
    var inputs = $(':text.save-input');
    for (var i=0; i<inputs.length; i++) {
        var el = inputs[i];
        if (el.value!=originalValFor[el.name]) {
            $(el).css('background-color', '#ffffbf');
        }
        else {
            $(el).css('background-color', 'white');
        }
    }
    inputs = $(':radio.save-input');
    for (var i=0; i<inputs.length; i++) {
        var el = inputs[i];
        if (el.value!=originalValFor[el.name] && el.checked==true) {
            $(el).parent('label').css('color', 'red');
        }
        else {
            $(el).parent('label').css('color', '#333');
        }
    }
}

function show_selected_dbs(str, db_name_change) {
    var display_str = '';
    if (str!='') {
        var dbs = str.split(',');
        for (var i=0; i<dbs.length; i++) {
            if (i==0) {
                display_str += '<li class="list-group-item">' + db_name_change[dbs[i]] + ' <span class="glyphicon glyphicon-star"></span></li>'
            }
            else {
                display_str += '<li class="list-group-item">' + db_name_change[dbs[i]] + '</li>'
            }
        }
    }
    $('#show-selected-templates').html(display_str);
}

function ScrollToResult() {
    $('html,body').animate({
        scrollTop: $('#result').offset().top,
    }, 1000);
};

function move_webpage_when_browsing_result_panels() {
    // In the result panels, improve web page move (for save view)
    $('#primers-result').on('shown.bs.collapse', function () {
        var offset = $(this).find('.collapse.in').prev('.panel-heading');
        if(offset) {
            $('html,body').animate({
                scrollTop: $(offset).offset().top-20
            }, 1000); 
        }
    });
}

function basename(path) {
    return path.replace(/.*\//, '');
}

function visualize(json_data) {
    // data
    var result_data = JSON.parse(json_data);
    $('#result').removeClass('hidden');
    $('#running-modal').modal('hide');
    $('#primers-result').html('');
    ScrollToResult();

    // no results
    if ('error' in result_data) {
        $('#primers-result').append($('#primers-result-template-error').html()).find('.alert-danger')
            .append('<h5><strong>ERROR</strong></h5>'+result_data['error']);
        $('#btn-download-tsv,#btn-download-json').prop('disabled', true);
        return;
    }
    
    // meta and primers
    var visualize_mode = result_data['meta']['mode'];
    var selected_dbs = result_data['meta']['dbs'].map(x=>basename(x)).join(',');
    var primer_data = result_data['primers'];
    let region_type = result_data['meta']['region_type'];

    // valid results
    generate_html_result(selected_dbs, db_name_change, primer_data, visualize_mode);
    if (result_data['meta']['check_multiplex']==false) {
        $('#dimers-result').next('.alert-danger').remove();
        $('#dimers-result').addClass('hidden');
    }
    else {
        generate_dimer_result(result_data['dimers']);
        $('#dimers-result').removeClass('hidden');
    }
    if (visualize_mode=='full') {
        GenerateGraph($('#site-1'), region_type, true);
        $('#primers-result .collapse').on('shown.bs.collapse', function () {
            GenerateGraph($(this), region_type, true);
        });
        $('#primers-result .collapse').on('hidden.bs.collapse', function () {
            $(this).find('.PrimerFigure').html('');
            $(this).find('.PrimerFigureControl').remove();
        });
    }
    else if (visualize_mode=='design') {
        GenerateGraph($('#site-1'), region_type, false);
        $('#primers-result .collapse').on('shown.bs.collapse', function () {
            GenerateGraph($(this), region_type, false);
        });
        $('#primers-result .collapse').on('hidden.bs.collapse', function () {
            $(this).find('.PrimerFigure').html('');
            $(this).find('.PrimerFigureControl').remove();
        });
    }
    move_webpage_when_browsing_result_panels();
    $('#btn-download-tsv,#btn-download-json').prop('disabled', false);
}

var json_data = '';
function AjaxSubmit(selected_dbs, mode) {
    $('#running-modal .modal-body h4').html('<span class="fa fa-spinner fa-spin fa-4x"></span>');
    $('#running-modal .progress-bar').css('width', '0%').html('');
    $('#running-modal').modal('show');
    var currentAjax = $.post($SCRIPT_ROOT + '/run', $('#form-primer').serialize(), function(data){
        json_data = data;
        visualize(json_data);
    });

    // Allow users to stop their running
    $('#stop').click(function(){
        if (currentAjax) {
            currentAjax.abort();
        }
        $.get($SCRIPT_ROOT + '/stop', function(){
            $('#running-modal').modal('hide');
        });
    });
}

// ***********************  Begin main functions  ***********************  
// ********* Page load *********
var db_name_change = new Object;
var config_data = new Object;
var originalValFor = new Object;
$(function () {
    $.get($SCRIPT_ROOT + '/dbinfo', function(data){
        config_data = JSON.parse(data);
    });

    $.get($SCRIPT_ROOT + '/dbselect', function(data){
        // get all the default values
        var inputs = $(':text.save-input');
        for (var i=0; i<inputs.length; i++) {
            var el = inputs[i];
            originalValFor[el.name] = el.defaultValue;
        }
        inputs = $(':radio.save-input');
        for (var i=0; i<inputs.length; i++) {
            var el = inputs[i];
            originalValFor[el.name] = 0;
        }

        // init menu, get db_name_change
        db_name_change = menu_init(data);

        // load last save
        $('.save-input').phoenix();

        // refresh menu (and show databases in paramter column)
        $('select').selectpicker('refresh');

        // refresh running mode
        mode = $("[name='app-type']").val()
        $('a[href="#'+mode+'"]').tab('show');

        highlight_changed_field(originalValFor);

        // test('example.fa', 'primer_full.json', 'full');
    });

    tooltip_init();

});

// ********* Highlight differences when user changes *****************************************
$('input').blur(function(){
    highlight_changed_field(originalValFor);
})

// ********* Switch Running mode (init load change or by user change) and update UI *********
$('a[data-toggle="tab"]').on('shown.bs.tab', function (e) {

    // change the hidden field value : app-type
    var mode = $(e.target).attr('href').replace('#', '');
    $("[name='app-type']").val(mode);

    // change parameters fieldset
    switch_parameter_fieldset(mode);

    // change textarea
    switch (mode) {
        case 'full':
            textarea_color = '#bfffdf';
            break;
        case 'design':
            textarea_color = '#fcd8f7';
            break;
        case 'check':
            textarea_color = '#ffdfbf';
            break;
        default:
            textarea_color = '#bfffdf';
            break;
    }
    $('textarea').css('background-color', textarea_color);

    // qPCR fieldset
    main_db = $('[name="selected_dbs"]').val().split(',')[0];
    region_type = $('[name="region_type"]').val();
    qPCR_filedset(config_data, main_db, mode, region_type);
});

// ********* Switch databases (init load or by user change) and update UI *********
var vals = [];
$('[name="templates"]').on('changed.bs.select refreshed.bs.select', function (event) {
    if (event.type=='refreshed') {  // init load
        selected_dbs = $('[name="selected_dbs"]').val();
    }
    else {  // user change
        for(var i=0; i <$('[name="templates"] option').length; i++) {
            if ($($('[name="templates"] option')[i]).prop('selected') ) {
                if (!vals.includes(i)) {
                    vals.push(i);
                }
            } 
            else if (vals.includes(i)) {
                vals.splice(vals.indexOf(i), 1);
            }
        }
        var selected_dbs = '';
        vals.forEach(function(ele) {
            selected_dbs += $($('[name="templates"] option')[ele]).val() + ',';
        })
        selected_dbs = selected_dbs.replace(/,$/, '');
        $('[name="selected_dbs"]').val(selected_dbs);
    }
    show_selected_dbs(selected_dbs, db_name_change);

    // qPCR fieldset
    main_db = selected_dbs.split(',')[0];
    mode = $("[name='app-type']").val();
    region_type = $('[name="region_type"]').val();
    qPCR_filedset(config_data, main_db, mode, region_type);
});

// ********* Switch region type and update UI ****************************************
$('[name="region_type"]').on('changed.bs.select', function(){
    main_db = $('[name="selected_dbs"]').val().split(',')[0];
    mode = $("[name='app-type']").val();
    region_type = $('[name="region_type"]').val();
    qPCR_filedset(config_data, main_db, mode, region_type);
});

// ********* The reset buttton *******************************************************
$(':reset').click(function(){
    $('[name="templates"]').selectpicker('val', '');
    $('[name="region_type"]').selectpicker('val', 'SEQUENCE_TARGET');
    $('[name="selected_dbs"]').val('');
    $('#show-selected-templates').html('');
    highlight_changed_field(originalValFor);
});

// ********* User submit the form *******************************************************
$('#form-primer').validationEngine('attach', {
    onValidationComplete: function(form, status) {
        if (status) {
            var selected_dbs = $('[name="selected_dbs"]').val();
            var mode = $("[name='app-type']").val();
            AjaxSubmit(selected_dbs, mode);
        }
    }
});

// ********* The running modal *******************************************************
var evtSource;
$('#running-modal').on('shown.bs.modal', function(){
    evtSource = new EventSource($SCRIPT_ROOT + '/monitor');
    evtSource.onmessage = function(e) {
        var progress_data = JSON.parse(e.data);
        if (progress_data.current_task=='design') {
            $('#running-modal .modal-body h4').html('Design Primers: '+ progress_data.all_tasks_num + ' sites');
        }
        if (progress_data.current_task=='blast') {
            $('#running-modal .modal-body h4').html('Waiting for BLAST: '+ progress_data.all_tasks_num + ' tasks');
        }
        if (progress_data.current_task=='multiplex') {
            $('#running-modal .modal-body h4').html('Checking Multiplex Dimers: '+ progress_data.all_tasks_num + ' tasks');
        }
        var progress_per = progress_data.complete_count/progress_data.all_tasks_num*100
        $('#running-modal .progress-bar').css('width', progress_per+'%')
            .html(progress_data.complete_count +' tasks finished');
        if (progress_data.current_task=='finish') {
            $('.progress-bar').removeClass('active').html('Completed. Waiting for generating results...');
            evtSource.close();
        }
    }
});
$('#running-modal').on('hidden.bs.modal', function(){
    evtSource.close();
});

// ********* Download ******************************************************************
$('#btn-download-json').click(function(){
    var blob = new Blob([json_data], {type: "application/json;charset=utf-8"});
    saveAs(blob, "primers.json"); 
});
$('#btn-download-tsv').click(function(){
    $.post($SCRIPT_ROOT + '/download_tsv', {json: json_data, dbs: $('[name="selected_dbs"]').val()}, function(data){
        var blob = new Blob([data], {type: "text/plain;charset=utf-8"});
        saveAs(blob, "primers.txt"); 
    });
});

// ********* Upload and visulization *****************************************************
$('#btn-visulization').click(function(){
    $('#file-visulization').click();
})
$('#file-visulization').on('change', function(){
    files = $('#file-visulization')[0].files;
    var reader = new FileReader();
    reader.readAsText(files[0]);
    reader.onload = function(e){
        json_data = e.target.result;
        visualize(json_data);
    };
})

// ********* amplicon details when click the button ****************************************
function target_seq_format(qseq, sseq){
    qseq = qseq.toUpperCase();
    sseq = sseq.toUpperCase();
    new_sseq = ''
    for (i in sseq) {
        if (sseq[i]==qseq[i]) {
            new_sseq +='.'
        }
        else {
            new_sseq +=sseq[i]
        }
    }
    return new_sseq
}
$('#amplicon-details-modal').on('show.bs.modal', function (event) {
    var button = $(event.relatedTarget);
    var amplicon_details = button.data('amplicon');
    $('#amplicon-details-modal pre').html('');
    $('#amplicon-details-modal pre').append('Template Region: '+amplicon_details['region']+'<br/>')
        .append('Primer Left: '+amplicon_details['plus']['qseqid']+' ('+amplicon_details['plus']['qseq']+')'+'<br/>')
        .append('Primer Right: '+amplicon_details['minus']['qseqid']+' ('+amplicon_details['minus']['qseq']+')'+'<br/>')
        .append('Product Size: '+amplicon_details['product_size']+'<br/>')
        .append('Melting Temperature for Left Primer (°C): '+amplicon_details['plus']['Tm']+'<br/>')
        .append('Melting Temperature for Right Primer (°C): '+amplicon_details['minus']['Tm']+'<br/>')
        .append('Amplicon Seq: "'+amplicon_details['product_seq']+'"<br/><br/>')
        .append('Primer Left      '+' '.repeat(amplicon_details['plus']['sstart'].toString().length-1)+'1 '+amplicon_details['plus']['qseq']
            +' '+amplicon_details['plus']['qseq'].toString().length+'<br/>')
        .append('Binding Template '+amplicon_details['plus']['sstart'].toString()+' '+target_seq_format(amplicon_details['plus']['qseq'], amplicon_details['plus']['sseq'])
            +' '+amplicon_details['plus']['send'].toString()+'<br/><br/>')
        .append('Primer Right     '+' '.repeat(amplicon_details['minus']['sstart'].toString().length-1)+'1 '+amplicon_details['minus']['qseq']
            +' '+amplicon_details['minus']['qseq'].toString().length+'<br/>')
        .append('Binding Template '+amplicon_details['minus']['send'].toString()+' '+target_seq_format(amplicon_details['minus']['qseq'], amplicon_details['minus']['sseq'])
            +' '+amplicon_details['minus']['sstart'].toString()+'<br/>')

});