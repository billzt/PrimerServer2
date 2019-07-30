/***************** Complex functions to display Figure after showing Panel, used after the server return results */
function getMaxOfArray(numArray) {
    return Math.max.apply(null, numArray);
}
function getMinOfArray(numArray) {
    return Math.min.apply(null, numArray);
}
function GenerateGraph(el) {
    // empty the element
    el.find('.PrimerFigure').html('');
    
    // svg
    var svg = d3.select(el.find('.PrimerFigure')[0])
                .append('svg')
                .attr('width', '100%');
    svg.append('style').text('.axis path,.axis line{fill: none;stroke: blue;stroke-width: 2px;shape-rendering: crispEdges;} '
        +'.axis text{font-size:11px;} .primerLabel{font-size:11px;} .primerUniqueLabel{font-size:11px;font-weight:bold;}');
    
    // primers regions
    var primers = el.find('.list-group-item-primer');
    var primers2region = new Object;
    var primers2hit = new Object;
    var allPoses = new Array;
    for (var i=0; i<primers.length; i++) {
        var id = $(primers[i]).find('.list-group-item-heading').attr('id');
        var primer_left = $(primers[i]).find('.primer-left-region');
        var primer_right = $(primers[i]).find('.primer-right-region');
        var region_1 = $(primer_left[0]).html().split('-');
        var region_2 = $(primer_right[0]).html().split('-');
        primers2region[id] = [region_1, region_2];
        Array.prototype.push.apply(allPoses, region_1);
        Array.prototype.push.apply(allPoses, region_2);
        var hitNum = $(primers[i]).find('.hit-num').data('hit');
        primers2hit[id] = hitNum;
    }
    
    // axis
    var axisStart;
    var axisEnd;
    var targetPos = el.prev().find('.site-detail').data('pos');
    var targetLen = el.prev().find('.site-detail').data('length');
    if ($('[name="region_type"]:checked').val()=='SEQUENCE_INCLUDED_REGION') {
        axisStart = targetPos-Math.round(targetLen/5)>0 ? targetPos-Math.round(targetLen/5) : 1;
        axisEnd = targetPos+targetLen+Math.round(targetLen/5);
    }
    else if ($('[name="region_type"]:checked').val()=='SEQUENCE_TARGET' || 
    $('[name="region_type"]:checked').val()=='FORCE_END') {
        axisStart = getMinOfArray(allPoses);
        axisEnd = getMaxOfArray(allPoses);
    }
    var axisScale = d3.scaleLinear().domain([axisStart, axisEnd]).range([0, 1000]);
    var axis = d3.axisTop().scale(axisScale).ticks(10);
    svg.append('g').attr('class', 'axis').call(axis); // axis: translate(x,y) is no longer needed as PanZoom can do it
    
    // Text
    var template = el.prev().find('.site-detail').data('seq');
    svg.append('text').attr('x','0').attr('y','-30').text('Template '+template).attr('font-size', '120%');
    
    // target region
    var rectHight = 60;
    svg.append('rect').attr('x', axisScale(targetPos))  
       .attr('y', -rectHight/2)  // axis:y-rect:height/2
       .attr('width', axisScale(targetPos+targetLen)-axisScale(targetPos)).attr('height', rectHight)
       .attr("fill", "none").attr('stroke', 'red').attr('stroke-width', '3');
    
    // Primer Group
    var colorScale = d3.scaleLinear().domain([1, 100]).range([0, 32]);
    function AddPrimer(LprimerStart, LprimerEnd, RprimerStart, RprimerEnd, i, h, id) {
        var primerGroup = svg.append('a').attr('xlink:href','#'+id).attr('class', 'primerGroup')
                        .attr('title', 'Primer '+i).append('g');
        var baseY = rectHight+30*(i-1);
        var lineFunction = d3.line()
            .x(function(d) { return Math.round(axisScale(d.x)); })
            .y(function(d) { return d.y; });
        var color = 'rgb('+Math.round(colorScale(h)*8)+','+Math.round(colorScale(h)*8)+','+Math.round(colorScale(h)*8)+')';
        
        // Left Primer
        var Llength = LprimerEnd-LprimerStart+1;
        
        var LlineData = [ { "x": LprimerStart, "y": baseY-5},  { "x": LprimerStart+Math.round(Llength/3*2), "y": baseY-5},
                         { "x": LprimerStart+Math.round(Llength/3*2), "y": baseY-10}, {"x": LprimerEnd, "y": baseY},
                         { "x": LprimerStart+Math.round(Llength/3*2), "y": baseY+10}, {"x": LprimerStart+Math.round(Llength/3*2), "y": baseY+5},
                         { "x": LprimerStart, "y": baseY+5},  { "x": LprimerStart, "y": baseY-5}];
        
        primerGroup.append('path').attr('d', lineFunction(LlineData)).attr("fill", color).attr('stroke', color);
        
        
        // Right Primer
        var Rlength = RprimerEnd-RprimerStart+1;
        
        var RlineData = [ {"x": RprimerStart, "y": baseY}, {"x": RprimerStart+Math.round(Rlength/3*1), "y": baseY-10},
                          {"x": RprimerStart+Math.round(Rlength/3*1), "y": baseY-5}, {"x": RprimerEnd, "y": baseY-5},
                          {"x": RprimerEnd, "y": baseY+5}, {"x": RprimerStart+Math.round(Rlength/3*1), "y": baseY+5},
                          {"x": RprimerStart+Math.round(Rlength/3*1), "y": baseY+10}, {"x": RprimerStart, "y": baseY}];
        primerGroup.append('path').attr('d', lineFunction(RlineData)).attr("fill", color).attr('stroke', color);
        
        // Center Line
        var LineData = [{"x": LprimerEnd, "y": baseY}, {"x": RprimerStart, "y": baseY}];
        primerGroup.append('path').attr('d', lineFunction(LineData)).attr("fill", color).attr('stroke', color);
        
        // Text
        var primerLabel = h==1 ? '[Primer'+i+'] '+h+' Amplicon' : '[Primer'+i+'] '+h+' Amplicons';
        if (h==1) {
            primerGroup.append('text').attr("x", axisScale(LprimerStart)-120)
            .attr("y", baseY+5).attr('class', 'primerUniqueLabel').attr('fill', 'red').text(primerLabel);
        }
        else {
            primerGroup.append('text').attr("x", axisScale(LprimerStart)-120)
            .attr("y", baseY+5).attr('class', 'primerLabel').text(primerLabel);
        }
    }

    var primerRank = 1;
    for ( id in primers2region) {
        var LprimerStart=primers2region[id][0][0];
        var LprimerEnd=primers2region[id][0][1];
        var RprimerStart=primers2region[id][1][0];
        var RprimerEnd=primers2region[id][1][1];
        var h = primers2hit[id];
        AddPrimer(LprimerStart*1, LprimerEnd*1, RprimerStart*1, RprimerEnd*1, primerRank*1, h*1, id);
        primerRank++;
        //break;
    }
    
    // extend svg height if there are too many primers
    if (primerRank>3) {
        svg.attr('height', (primerRank-3)*30+150);
    }
    
    // Pan and Zoom
    if ($('.PrimerFigure svg').length>0) {
        var zoomObj = svgPanZoom('.PrimerFigure svg');
        if (el.find('.PrimerFigureControl').length==0) {
            el.find('.PrimerFigure').before('<div class="PrimerFigureControl btn-group">'
                +'<button type="button" class="btn btn-default zoom-in" title="Zoom in" data-toggle="tooltip"><i class="fa fa-search-plus"></i></button>'
                +'<button type="button" class="btn btn-default zoom-out" title="Zoom out" data-toggle="tooltip"><i class="fa fa-search-minus"></i></button>'
                +'<button type="button" class="btn btn-default zoom-reset" title="Reset" data-toggle="tooltip"><i class="fa fa-refresh"></i></button>'
                +'<button type="button" class="btn btn-default svg-download" title="Download SVG" data-toggle="tooltip"><i class="fa fa-download"></i></button>'
                +'</div>');                
        }
        $('[data-toggle="tooltip"]').tooltip({html: true});
        $(window).resize(function(){
            zoomObj.resize();
            zoomObj.fit();
            zoomObj.center();
        });
        $('.zoom-in').click(function(){
            zoomObj.zoomIn();
        });
        $('.zoom-out').click(function(){
            zoomObj.zoomOut();
        });
        $('.zoom-reset').click(function(){
            zoomObj.resetZoom().resetPan();
        });
        $('.svg-download').click(function(){
            var download_text = '<svg width="100%" xmlns="http://www.w3.org/2000/svg" '
            +'xmlns:xlink="http://www.w3.org/1999/xlink" xmlns:ev="http://www.w3.org/2001/xml-events"> '
            + $('.PrimerFigure svg').html() + '</svg>';
            var blob = new Blob([download_text], {type: "text/svg;charset=utf-8"});
            saveAs(blob, "primer.svg"); 
        });
    }
}

/***************** Complex functions finished  ***********************************************************/

// function for save results
function toURL(value, row) {
    return '<a target="_blank" href="save/'+row.url+'">'+value+'</a>';
}

// function to show selected databases
function dbShow(str, group_data) {
    if (!str) {
        return '';
    }
    var dbs = str.split(',');
    var display_str = '';
    // db file name to descriptions
    var db_name_change = new Object;
    for (group in group_data){
        for (template in group_data[group]) {
            db_name_change[template] = group_data[group][template]['desc']
        }
    }
    for (var i=0; i<dbs.length; i++) {
        if (i==0) {
            display_str += '<li class="list-group-item">' + db_name_change[dbs[i]] + ' <span class="glyphicon glyphicon-star"></span></li>'
        }
        else {
            display_str += '<li class="list-group-item">' + db_name_change[dbs[i]] + '</li>'
        }
    }
    return display_str;
}

$(function () {
    // Tooltip for bootstrap
    $('[data-toggle="tooltip"]').tooltip({html: true});
    $('[data-toggle="popover"]').popover({html: true});
    $('[data-toggle="table"]').on('post-body.bs.table', function () {
        $('[data-toggle="tooltip"]').tooltip({
            container: 'body'
        });
    });
    
    // select template: options
    var originalValFor = new Object;
    var group_data = new Object;
    $.get('/dbselect', function(data){
        group_data = JSON.parse(data);
        for (group in group_data) {
            $('[name="templates[]"]').append('<optgroup label="'+group+'">');
            for (template in group_data[group]) {
                $('optgroup[label="'+group+'"]').append('<option data-subtext="seq IDs:'+group_data[group][template]['IDs']
                    +'" value="'+template+'">'+group_data[group][template]['desc']+'</option>');
            }
        }

        // get all the default values
        var inputs = $(':text.save-input');
        for (var i=0; i<inputs.length; i++) {
            var el = inputs[i];
            originalValFor[el.name] = el.defaultValue;
        }

        // Load user's last saved inputs
        $('.save-input').phoenix({
            saveInterval: 1000,
        });
        $('[name="templates[]"]').selectpicker('refresh');
        
        // Highlight Changed Field
        for (var i=0; i<inputs.length; i++) {
            var el = inputs[i];
            if (el.value!=originalValFor[el.name]) {
                $(el).css('background-color', '#ffffbf');
            }
        }
        
        // Highlight Changed Field after users' input
        $(inputs).blur(function(){
            for (var i=0; i<inputs.length; i++) {
                var el = inputs[i];
                if (el.value!=originalValFor[el.name]) {
                    $(el).css('background-color', '#ffffbf');
                }
                else {
                    $(el).css('background-color', 'white');
                }
            }
        })
        
        // In help modal, list IDs (After the options completely loaded)
        var options = $('[name="templates[]"] option');
        for (var i=0; i<options.length; i++) {
            var option_val = $(options[i]).val();
            if (option_val!='' && option_val!='custom') {
                var option_text = $(options[i]).html();
                $('#help-modal-ID-list').append('<a href="/dbdownload/'+option_val+'" class="list-group-item">'
                +option_text+' <span class="glyphicon glyphicon-download"></span></a>');
            }
        }
        
        // App type: design OR check
        // when the document is ready, show the last App type select by users
        $('a[href="#'+$("[name='app-type']").val()+'"]').tab('show');
        var app_type_title_save = ' ';
        switch($("[name='app-type']").val()) {
            case 'design': app_type_title_save = 'Design Primers';break;
            case 'check': app_type_title_save = 'Check Primers';break;
        }
        $('#save-webpage input').val(app_type_title_save);
        
        // Paramter panel collapse
        if ($("[name='parameter-design']").val()=='hide') {
            $('#parameter-design').collapse('hide');
        }
        else {
            $('#parameter-design').collapse('show');
        }
        if ($("[name='parameter-check']").val()=='hide') {
            $('#parameter-check').collapse('hide');
        }
        else {
            $('#parameter-check').collapse('show');
        }
    });
    
    
    // modify reset button to satisfy selector
    $(':reset').click(function(){
        $('[name="templates[]"]').selectpicker('val', '');
        localStorage.setItem('primer-templates:'+location.href, '');
        $('#show-selected-templates').html('');
    });
    
    // change App type when user change it
    $('a[data-toggle="tab"]').on('shown.bs.tab', function (e) {
        var type = $(e.target).attr('href').replace('#', '');
        $("[name='app-type']").val(type);
        if (type=='saved') {
            $('#parameter-panel').addClass('hidden');
            $(':submit').addClass('hidden');
            $(':reset').addClass('hidden');
            $('#result').addClass('hidden');
            $('#footer-button').addClass('hidden');
        }
        else {
            var app_type_title_save = ' ';
            switch(type) {
                case 'full': {
                    $('#parameter-design').removeClass('hidden');
                    $('#parameter-check').removeClass('hidden');
                    break;
                }
                case 'design': {
                    $('#parameter-design').removeClass('hidden');
                    $('#parameter-check').addClass('hidden');
                    app_type_title_save = 'Design Primers';
                    break;
                }
                case 'check': {
                    $('#parameter-design').addClass('hidden');
                    $('#parameter-check').removeClass('hidden');
                    app_type_title_save = 'Check Primers';
                    break;
                } 
            }
            $('#save-webpage input').val(app_type_title_save);
            $('#parameter-panel').removeClass('hidden');
            $(':submit').removeClass('hidden');
            $(':reset').removeClass('hidden');
            $('#result').removeClass('hidden');
            if ($('#result').html()!='') {
                $('#footer-button').removeClass('hidden');
            }
        }
    });
    
    // change panel collapse state when user change it
    $('#parameter-design').on('shown.bs.collapse', function(){
        $("[name='parameter-design']").val('show');
        $("[data-target='#parameter-design'] i").removeClass('fa-caret-up').addClass('fa-caret-down');
    });
    $('#parameter-design').on('hidden.bs.collapse', function(){
        $("[name='parameter-design']").val('hide');
        $("[data-target='#parameter-design'] i").removeClass('fa-caret-down').addClass('fa-caret-up');
    });
    $('#parameter-check').on('shown.bs.collapse', function(){
        $("[name='parameter-check']").val('show');
        $("[data-target='#parameter-check'] i").removeClass('fa-caret-up').addClass('fa-caret-down');
    });
    $('#parameter-check').on('hidden.bs.collapse', function(){
        $("[name='parameter-check']").val('hide');
        $("[data-target='#parameter-check'] i").removeClass('fa-caret-down').addClass('fa-caret-up');
    });
    
    // Remove flanking blanks after text input; If it is blank, fill original value for it
    $(':text').blur(function(){
        var val = $.trim($(this).val());
        if (val!='') {
            $(this).val(val);
        }
        else {
            var el = $(this);
            $(this).val(originalValFor[el[0].name]);
        }
    });

    // If users select (Or inintially load) custom template, then showing custom template FASTA sequence input textarea
    var selectedDBNum;
    $('[name="templates[]"]').on('changed.bs.select refreshed.bs.select', function (event) {
        var selectedOptions = event.target.selectedOptions;
        var hideCustomDB = 1;
        selectedDBNum = selectedOptions.length;
        for (var i=0; i<selectedDBNum; i++) {
            if (selectedOptions[i].value=='custom') {
                hideCustomDB = 0;
            }
        }
        if (hideCustomDB==0) {
            $('[name="custom-db-sequences"]').parent().removeClass('hidden');
        }
        else {
            $('[name="custom-db-sequences"]').parent().addClass('hidden');
        }
    });

    // Show initial selected databases
    var selected_dbs = localStorage.getItem('primer-templates:'+location.href);
    $('#show-selected-templates').html(dbShow(selected_dbs, group_data));
    $('[name="templates"]').val(selected_dbs);
    if (!selected_dbs) {    // inintialize
        $('[name="templates[]"]').on('refreshed.bs.select', function (event) {
            var selectedOptions = event.target.selectedOptions;
            selected_dbs = '';
            for (var i=0; i<selectedDBNum; i++) {
                selected_dbs += selectedOptions[i].value+',';
            }
            selected_dbs = selected_dbs.replace(/,$/, '');
            localStorage.setItem('primer-templates:'+location.href, selected_dbs);
            $('#show-selected-templates').html(dbShow(selected_dbs, group_data));
            $('[name="templates"]').val(selected_dbs);
        });
    }
    
    // Change selected databases
    var vals = [];
    $('[name="templates[]"]').change(function (event) {
        for(var i=0; i <$('[name="templates[]"] option').length; i++) {
            if ($($('[name="templates[]"] option')[i]).prop('selected') ) {
                if (!vals.includes(i)) {
                    vals.push(i);
                }
            } 
            else if (vals.includes(i)) {
                vals.splice(vals.indexOf(i), 1);
            }
        }
        selected_dbs = '';
        vals.forEach(function(ele) {
          selected_dbs += $($('[name="templates[]"] option')[ele]).val() + ',';
        })
        selected_dbs = selected_dbs.replace(/,$/, '');
        $('#show-selected-templates').html(dbShow(selected_dbs, group_data));
        $('[name="templates"]').val(selected_dbs);
        localStorage.setItem('primer-templates:'+location.href, selected_dbs);
    });
    

    // form validation & submit
    function ScrollToResult() {
        $('html,body').animate({
            scrollTop: $('#result').offset().top,
        }, 1000);
    };
    function AjaxSubmit() {
        $('#running-modal .modal-body h4').html('<span class="fa fa-spinner fa-spin fa-4x"></span>');
        $('#running-modal .progress-bar').css('width', '0%').html('');
        $('#running-modal').modal('show');
        var currentAjax = $.post('script/primer.php', $('#form-primer').serialize(), function(data){
            $('#result').html(data);
            $('#running-modal').modal('hide');
            ScrollToResult();
            // call Complex functions if we are in design & check mode
            if ($('[name="app-type"]').val()=='design') {
                GenerateGraph($('#site-1'));
                $('#primers-result').find('.collapse').on('shown.bs.collapse', function (e) {
                    GenerateGraph($(this));
                });
                $('#primers-result').find('.collapse').on('hidden.bs.collapse', function (e) {
                    $(this).find('.PrimerFigure').html('');
                    $(this).find('.PrimerFigureControl').remove();
                });                    
            }

            // show download area
            if ($('#primers-result').length>0) {
                $('#footer-button').removeClass('hidden');
            }
            else {
                $('#footer-button').addClass('hidden');
            }
            
            // remove primer amplicon links if no hits were found
            $('[data-hit=0]').children('[data-toggle="modal"]').addClass('hidden');
            
            // In the result panels, improve web page move (for save view)
            $('#primers-result').on('shown.bs.collapse', function () {
                var offset = $(this).find('.collapse.in').prev('.panel-heading');
                if(offset) {
                    $('html,body').animate({
                        scrollTop: $(offset).offset().top-20
                    }, 1000); 
                }
            });
        });
        
        // Allow users to stop their running
        $('#stop').click(function(){
            if (currentAjax) {
                currentAjax.abort();
            }
            $.get('script/modal_stop.php', function(){
                $('#running-modal').modal('hide');
            });
        });
    };
    $('#form-primer').validationEngine('attach', {
        autoHidePrompt: true,
        autoHideDelay: 5000,
        onFieldFailure: function(field) {
            if (field) {    // A bug: it will return an extra undef field
                field.parents('.panel').find('[data-target]').removeAttr('data-toggle');
                field.parents('.panel').find('i').css('cursor', 'not-allowed');
            }
        },
        onFieldSuccess: function(field) {
            if (field) {
                field.parents('.panel').find('[data-target]').attr('data-toggle', 'collapse');
                field.parents('.panel').find('i').css('cursor', 'pointer');
            }
        },
        onValidationComplete: function(form, status) {
            if (status) {
                AjaxSubmit();
            }
        }
    });

    // Showing Specificity result in modal dynamically by Ajax 
    $('#specificity-check-modal').on('show.bs.modal', function (event) {
        var button = $(event.relatedTarget); // Button that triggered the modal
        var fileName = button.data('whatever');
        var modal = $(this);
        var target_band_size = button.data('targetsize');
        $.get('script/modal_specificity_result.php', {file: fileName}, function(data) {
            var result_data = JSON.parse(data);
            modal.find('.modal-body .fa-spinner').addClass('hidden');
            modal.find('.modal-body pre').html(result_data.file);
            
            /*************** Virtual electrophoresis ***************************/
            var sizes = result_data.sizes.map(function(x){return x*1});
            var gel = electrophoresis();
            $("#ve").html("");
            var svg = d3.select("#ve").append("svg")
              .attr("width", 240)
              .attr("height", 600)
              .call(gel.makeGel); // Make black background first.
            var words = [[100, 250, 500, 750, 1000, 2000],sizes];
            var names = ["Marker", "Amplicons"];
            // You can set every parameter by dictionary or method-chain. Example: {DNA: text} or .DNA(text).
            gel = electrophoresis().lane_number(2).tooltip_name_offsetX(-10).duration(1000).enzymes(words)
            .names(names).tooltip_band_offsetX(40).tooltip_band_offsetY(10);
            svg.call(gel);
          
            var bands = $('.gel-band text');
            for (var i=0; i<bands.length; i++) {
                var this_size = $(bands[i]).html();
                $(bands[i]).html(this_size+ ' bp');
                if (this_size == target_band_size) {
                    $(bands[i]).attr('fill','red').attr('font-weight','bold');
                }
            }
            /*************** Virtual electrophoresis finished ***************************/
        });
    })
    
    // When running, showing a progress bar
    var timer = $.timer(function(){
        $.get('script/modal_progress.php', function(data) {
            var progress_data = JSON.parse(data); // total, finished, percent
            if (progress_data.total>0) {
                $('#running-modal .modal-body h4').html(progress_data.total+' primer sequences generated. <br/>Waiting for BLAST: ' + progress_data.total
                + ' primer sequences &times; ' + selectedDBNum + ' databases' +' <span class="fa fa-spinner fa-spin"></span>');
                $('#running-modal .progress-bar').css('width', progress_data.percent+'%').html(progress_data.finished+' primer sequences finished');
                if (progress_data.percent>=100) {
                    timer.stop();
                    $('.progress-bar').removeClass('active').html('Completed. Waiting for generating results...');
                }                    
            }
        });
    });
    $('#running-modal').on('shown.bs.modal', function(){
        timer.set({ time : 1000, autostart : true });
    });
    $('#running-modal').on('hidden.bs.modal', function(){
        timer.stop();
    });
    
    // Download primers from Web UI
    $('#download-primer button.btn-primary').click(function(){
        // Get User Input
        var download_hit = $(':radio[name="download_hit"]:checked').val();   // 1: only download unique primers or 0
        var download_site = $(':radio[name="download_site"]:checked').val(); // 1: only download best primer per site or 0

        var download_text = "#Site_ID\tPrimer_Rank\tPenalty_Score\tTarget_Product_Size\tPossible_Hit_Num\tPrimer_Seq\r\n";
        var sites = $('#primers-result .panel');
        for (var i=0; i<sites.length; i++) {
            var site_id = $(sites[i]).find('.panel-heading').find('small');
            var primers = $(sites[i]).find('.panel-body').find('.list-group-item-primer');
            if (primers.length==0) {    // No primer for this site, only appear in design & check
                download_text += site_id.data('seq')+'-'+site_id.data('pos')+'-'+site_id.data('length')+"\tNo_Primer\r\n";
                continue;
            }
            if ($(sites[i]).find('.list-group-item-success').length==0 && download_hit==1) {    // No unique primers for this site
                if (site_id.data('seq')) {  // design & check
                    download_text += site_id.data('seq')+'-'+site_id.data('pos')+'-'+site_id.data('length')+"\t";
                }
                else {  // check only
                    download_text += site_id.html() + "\t";
                }
                download_text += "No_Unique_Primers\r\n";
                continue;
            }
            for (var j=0; j<primers.length; j++) {
                // deside whether to print this primer or not
                var hit_num = $(primers[j]).find('.hit-num').data('hit');
                if (hit_num==0 && download_hit==1) {
                    continue;
                }
                var primer_id = $(primers[j]).find('.list-group-item-heading').html();
                if (primer_id!='Primer 1' && download_site==1) {
                    continue;
                }
                
                // print site ID
                if (site_id.data('seq')) {  // design & check
                    download_text += site_id.data('seq')+'-'+site_id.data('pos')+'-'+site_id.data('length')+"\t";
                }
                else {  // check only
                    download_text += site_id.html() + "\t";
                }

                // print primers

                var primer_seqs = $(primers[j]).find('.list-group-item-text').find('.monospace-style');
                var penalty = $(primers[j]).find('.penalty');
                var product_size = $(primers[j]).find('a[href="javascript:void(0)"]').data('targetsize');
                
                download_text += primer_id + "\t";
                if (penalty.length>0) {
                    download_text += penalty.html() + "\t";
                }
                else {
                    download_text += "\t";
                }
                if (product_size) {
                    download_text += product_size + "\t";
                }
                else {
                    download_text += "\t";
                }
                download_text += hit_num + "\t";
                for (var k=0; k<primer_seqs.length; k++) {
                    download_text += $(primer_seqs[k]).html() + "\t";
                }
                download_text += "\r\n";
            }
        }
        var blob = new Blob([download_text], {type: "text/plain;charset=utf-8"});
        saveAs(blob, "primer.list.txt"); 
    });
    
    // Save webpage
    // system info
    $.get('script/current_time.php', function(data){
        $('#saved .alert-info').html(data);
    });
    
    // render table from localStorage
    var rows = [];
    var storage = window.localStorage;
    var storage_key_pattern = new RegExp('^PrimerServer.*'+location.href+'$');
    for (var i=0; i<storage.length; i++) {
        var storage_key = storage.key(i);
        if (storage_key_pattern.test(storage_key)) {
            var storage_key_array = storage_key.split('.');
            var thisURL = storage_key.replace('.'+location.href, '');
            rows.push({
                time: storage_key_array[1].replace(/-/g,'/')+' '+storage_key_array[2].replace(/-/g,':'),
                title: localStorage.getItem(storage_key),
                file_status: '<span class="file-status glyphicon" data-url="save/'+thisURL+'"></span>',
                url: thisURL,
            });
        }
    }
    
    // show file status in table
    $('#saved_table').on('post-body.bs.table', function (){
        var file_status = $('.file-status');
        for (var i=0; i<file_status.length; i++) {
            var fileURL = $(file_status[i]).data('url');
            $.ajax({
                url : fileURL,
                type : 'HEAD',
                statusCode: {
                    200: function() {
                            var thisURL = $(this)[0].url;
                            $('[data-url="'+thisURL+'"]').removeClass('glyphicon-remove').addClass('glyphicon-ok');
                         },
                    404: function() {
                            var thisURL = $(this)[0].url;
                            $('[data-url="'+thisURL+'"]').removeClass('glyphicon-ok').addClass('glyphicon-remove');
                         }
                }
            });
        }
    });
    
    $('#saved_table').bootstrapTable({data:rows});

    
    // control selection
    var selections = [];
    $('#saved_toolbar button').prop('disabled', true);
    $('#saved_table').on('check.bs.table check-all.bs.table uncheck.bs.table uncheck-all.bs.table', function (e, rows) {
        selections = $('#saved_table').bootstrapTable('getAllSelections');
        if (selections.length>0) {
            $('#saved_toolbar button').prop('disabled', false);
        }
        else {
            $('#saved_toolbar button').prop('disabled', true);
        }
    });
    
    // remove items
    $('#saved_toolbar button').click(function(){
        var remove_rows = $.map($('#saved_table').bootstrapTable('getSelections'), function (row) {
            return row.url;
        });
        $('#saved_table').bootstrapTable('remove', {
            field: 'url',
            values: remove_rows
        });
        for (var i=0; i<remove_rows.length; i++) {
            localStorage.removeItem(remove_rows[i]+'.'+location.href);
        }
        $.get('script/remove_webpage.php', {urls:remove_rows.join(' ')}, function(data){
            1;
        });
    });
    
    // store results
    $('#save-webpage button').click(function(){
        var app_type = $("[name='app-type']").val();
        var region_type = $('[name="region_type"]:checked').val();
        var save_title = $('#save-webpage input').val();
        $.get('script/save_webpage.php', {type1:app_type, type2:region_type}, function(data){
            localStorage.setItem(data+'.'+location.href, save_title); // data: PrimerServer.$date.$time.$session_id.html
            $('#save-result-modal').modal('show');
            var storage_key_array = data.split('.');
            $('#saved_table').bootstrapTable('prepend', [{
                time: storage_key_array[1].replace(/-/g,'/')+' '+storage_key_array[2].replace(/-/g,':'),
                title: save_title,
                url: data,
            }]);
        });
    });
    $('#save-result-modal-button').click(function(){
        $('#save-result-modal').modal('hide');
        $('a[href="#saved"]').tab('show');
    });
    
    // In the result panels, improve web page move (for save view)
    $('#primers-result').on('shown.bs.collapse', function () {
        var offset = $(this).find('.collapse.in').prev('.panel-heading');
        if(offset) {
            $('html,body').animate({
                scrollTop: $(offset).offset().top-20
            }, 1000); 
        }
    });
    
    // remove temporary files when use close or refresh the page
    $(window).bind("beforeunload", function() { 
        $.post('script/remove_tmp_files.php');
    });
});

