/***************** Complex functions to display Figure after showing Panel, used after the server return results */
function getMaxOfArray(numArray) {
    return Math.max.apply(null, numArray);
}
function getMinOfArray(numArray) {
    return Math.min.apply(null, numArray);
}
function GenerateGraph(el, region_type, specificity) {
    // empty the element
    el.find('.PrimerFigure').html('');

    // whether has primers?
    if (el.find('.alert-danger').length>0) {
        return;
    }
    
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
        var primer_internal = $(primers[i]).find('.primer-internal-region');
        var region_1 = $(primer_left[0]).html().split('-');
        var region_2 = $(primer_right[0]).html().split('-');
        primers2region[id] = [region_1, region_2];
        if (primer_internal.length>0) {
            var region_3 = $(primer_internal[0]).html().split('-');
            primers2region[id] = [region_1, region_2, region_3];
        }
        Array.prototype.push.apply(allPoses, region_1);
        Array.prototype.push.apply(allPoses, region_2);
        if (specificity==true) {
            var hitNum = $(primers[i]).find('.hit-num').data('hit');
            primers2hit[id] = hitNum;
        }
    }

    // axis
    var axisStart;
    var axisEnd;
    var targetPos = el.prev().find('.site-detail').data('pos');
    var targetLen = el.prev().find('.site-detail').data('length');
    if (region_type=='SEQUENCE_INCLUDED_REGION') {
        axisStart = targetPos-Math.round(targetLen/5)>0 ? targetPos-Math.round(targetLen/5) : 1;
        axisEnd = targetPos+targetLen+Math.round(targetLen/5);
    }
    else if (region_type=='SEQUENCE_TARGET' || region_type=='FORCE_END') {
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
    function AddPrimer(LprimerStart, LprimerEnd, RprimerStart, RprimerEnd, IprimerStart, IprimerEnd, i, h, id) {
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

        // optional internal probe
        if (IprimerStart != -1) {
            var Ilength = IprimerEnd-IprimerStart+1;
            var IlineData = [ {"x": IprimerStart, "y": baseY-2}, {"x": IprimerEnd, "y": baseY-2},
                              {"x": IprimerEnd, "y": baseY+2}, {"x": IprimerStart, "y": baseY+2},];
            primerGroup.append('path').attr('d', lineFunction(IlineData)).attr("fill", color).attr('stroke', color);
        }
        
        // Text
        if (h==1) {
            primerLabel = '[Primer'+i+'] '+h+' Amplicon';
            primerGroup.append('text').attr("x", axisScale(LprimerStart)-120)
            .attr("y", baseY+5).attr('class', 'primerUniqueLabel').attr('fill', 'red').text(primerLabel);
        }
        else if (h>1) {
            primerLabel = '[Primer'+i+'] '+h+' Amplicons';
            primerGroup.append('text').attr("x", axisScale(LprimerStart)-120)
            .attr("y", baseY+5).attr('class', 'primerLabel').text(primerLabel);
        }
        else {
            primerLabel = '[Primer'+i+'] ';
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
        var h = specificity==true ? primers2hit[id] : -1;
        var IprimerStart = -1;
        var IprimerEnd = -1;
        if (primers2region[id].length>2) {
            IprimerStart = primers2region[id][2][0];
            IprimerEnd = primers2region[id][2][1];
        }

        AddPrimer(LprimerStart*1, LprimerEnd*1, RprimerStart*1, RprimerEnd*1, IprimerStart*1, IprimerEnd*1, primerRank*1, h*1, id);
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
            + $('.PrimerFigure svg').html().replace(/NS\d+:href/gi, 'xlink:href') + '</svg>';
            var blob = new Blob([download_text], {type: "text/svg;charset=utf-8"});
            saveAs(blob, "primer.svg"); 
        });
    }
}

/***************** Complex functions finished  ***********************************************************/