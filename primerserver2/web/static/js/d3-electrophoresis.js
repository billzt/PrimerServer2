
function electrophoresis(config){
  // ---------------Dict arguments ----------------
  var config = (config === undefined) ? {} : config;
  config.DNA = (config.DNA === undefined) ? "AGCGCCCAATACGCAAACCGCCTCTCCCCGCGCGTTGGCCGATTCACGTTTACGAGTTGGAAACGA" : config.DNA;
  config.enzymes = (config.enzymes === undefined) ? [[1,5,10,15,20,30,40,50,70,100],"CAAA","GTTGG","GGATCC", "GA",["GT","CCT"],"ACC", "GAAA"] : config.enzymes;
  config.names = (config.names === undefined) ? ["marker", "1", "2", "3", "4", "5", "6", "7"]: config.names;
  config.scale = (config.scale === undefined) ? d3.scaleSqrt(): config.scale;
  config.lane_number = (config.lane_number === undefined) ? 8 : config.lane_number;
  config.gel_margin = (config.gel_margin === undefined) ? {top: 40, right: 20, bottom: 30, left: 20} : config.gel_margin;
  config.duration = (config.duration === undefined) ? 6000: config.duration;
  config.band_width = (config.band_width === undefined) ? 38 : config.band_width;
  config.band_blur = (config.band_blur === undefined) ? 2 : config.band_blur;
  config.band_thick_min = (config.band_thick_min === undefined) ? 1 : config.band_thick_min;
  config.band_thick_rate = (config.band_thick_rate === undefined) ? 0.04 : config.band_thick_rate;
  config.tooltip_name_on = (config.tooltip_name_on === undefined) ? "on" : config.tooltip_name_on;
  config.tooltip_name_size = (config.tooltip_name_size === undefined) ? 17 : config.tooltip_name_size;
  config.tooltip_name_offsetX = (config.tooltip_name_offsetX === undefined) ? 0 : config.tooltip_name_offsetX;
  config.tooltip_name_offsetY = (config.tooltip_name_offsetY === undefined) ? 20 : config.tooltip_name_offsetY;
  config.tooltip_band_on = (config.tooltip_band_on === undefined) ? "on" : config.tooltip_band_on;
  config.tooltip_band_size = (config.tooltip_band_size === undefined) ? 13 : config.tooltip_band_size;
  config.tooltip_band_offsetX = (config.tooltip_band_offsetX === undefined) ? 0 : config.tooltip_band_offsetX;
  config.tooltip_band_offsetY = (config.tooltip_band_offsetY === undefined) ? 0 : config.tooltip_band_offsetY;


  function elect(selection){
    // ---------------- Define functiions---------------
    function makeLenList(fragments, enzyme){
      var len_list = [];
      var enzyme_half1 = Math.round(enzyme.length/2);
      var enzyme_half2 = enzyme.length - enzyme_half1;

      fragments.forEach(function(fragment){
        var no_blank = fragment.replace(/\s*/g, "");
        len_list.push(fragment.length);
        });

      // Distribute Enzyme strings at each cut end to fill the deleted parts by split().
      if (len_list.length > 1){
        len_list.forEach(function(d,i){
          if (i != 0){
            len_list[i] = len_list[i] + enzyme_half1;
          }
          if (i != len_list.length-1){
            len_list[i] += enzyme_half2;
          }
        })
      }
      return len_list;
    }

    function cutDNA(DNA_seq, enzymes){
      var len_list = [];
      var fragments = [];
      if(typeof DNA_seq === "string"){
        DNA_seq = DNA_seq.toUpperCase();
      }

      // If enzyme is a list.
      if (enzymes.constructor === Array){
        var len_list_temp = [];
        fragments = [DNA_seq];
        enzymes.forEach(function(enzyme){
          if (typeof enzyme == "number"){
            len_list = len_list.concat(enzyme);
          }
          // Cut string by multiple strings.
          if (typeof enzyme == "string"){
            var fragments_temp = [];
            fragments = fragments.map(function(frag){
              enzyme = enzyme.toUpperCase();
              return frag.split(enzyme);
            });
            // Flatten lists.
            fragments = [].concat.apply([], fragments);
            len_list_temp = len_list_temp.concat(makeLenList(fragments, enzyme));
          }
        });

        len_list = len_list.concat(len_list_temp);
        return len_list;
      } // If enzyme is not a list.
      else {
        if (typeof enzymes === "number"){
          len_list.push(enzymes);
          return len_list;
        }
        if (typeof enzymes === "string"){
          enzymes = enzymes.toUpperCase();
          fragments = DNA_seq.split(enzymes);
          len_list = makeLenList(fragments, enzymes);
          return len_list;
        }
      }
    }

    // --------------- Make SVG ------------------
    function electrify(len_list, name_list, lane_number, band_width, band_duration,
      yScale, gel_margin, blur, thick_min, thick_rate,
      tooltip_name_on, tooltip_name_size, tooltip_name_offsetX, tooltip_name_offsetY,
      tooltip_band_on, tooltip_band_size, tooltip_band_offsetX,tooltip_band_offsetY) {

      var gel_width = selection.attr("width") - gel_margin.left - gel_margin.right;
      var gel_height = selection.attr("height") - gel_margin.top - gel_margin.bottom;

      var max_value = d3.max(d3.merge(len_list));
      yScale.domain([0, max_value]).range([gel_height,0]);

      // Filter for blur.
      var lane_filter = selection.append("defs").attr("class", "gel")
        .append("filter").attr("id", "gel-filter").attr("height","500%").attr("widht","100%")
        .append("feGaussianBlur").attr("in", "SourceGraphic").attr("stdDeviation", blur).attr("result", "blur")
        .append("feOffset").attr("result","offsetBlur").attr("dx",0).attr("dy",-5)
        .append("feBrend").attr("in","blur").attr("in2","offsetBlur").attr("mode","normal");

      gel= selection.append("g")
        .attr("class", "gel")
        .attr("transform", "translate(" + gel_margin.left + "," + gel_margin.top + ")");

      var lanes = gel.selectAll("g .gel-lane")
        .data(len_list).enter()
        .append("g")
        .attr("class", "gel-lane")
        .attr("transform", function(d,i){return "translate(" + i*(gel_width/lane_number) + "," + 0 + ")";});

      var rects = lanes.selectAll("g .gel-band")
        .data(function(d){return d;})
        .enter()
        .append("g")
        .attr("class","gel-band");

      var rectAttr = rects.append("rect")
        .attr("x", function(d){return 2;})
        .attr("y", function(d){return 0;})
        .attr("width", band_width)
        .attr("height", 12)
        .attr("fill", "white")
        .attr("filter", "url(#gel-filter)") // Apply blur filter.
        .transition()
        .attr("height", function(d){ return thick_min + ((gel_height-yScale(d)) * thick_rate); })
        .attr("y", function(d){ return yScale(d);})
        .duration(band_duration);

      if(tooltip_name_on == "on"){
        var tooltips_name = lanes.append("text")
          .attr("class","gel-tooltip")
          .attr("fill","white")
          .attr("font-size",tooltip_name_size)
          .attr("text-anchor", "middle")
          .attr("x", (gel_width/lane_number/3) + tooltip_name_offsetX)
          .attr("y",-(gel_margin.top - tooltip_name_offsetY))
          .text(function(d,i){ return name_list[i]; })
          .attr("opacity", 0);

          tooltips_name.transition()
          .duration(1000)
          .delay(band_duration)
          .attr("opacity", 1);
      }
      if(tooltip_band_on == "on"){
        var tooltips_band = rects.append("text")
          .attr("class","gel-tooltip")
          .attr("fill","white")
          .attr("font-size",tooltip_band_size)
          .attr("x",tooltip_band_offsetX)
          .attr("y",function(d){ return yScale(d) + tooltip_band_offsetY;})
          .text(function(d){return d;})
          .attr("opacity", 0);

          tooltips_band.transition()
          .duration(1500)
          .delay(band_duration+800)
          .attr("opacity", 1);
      }
    }

    // ----------------- Execute ------------------
    var len_list = [];
    config.enzymes.forEach(function(enzyme){
      len_list.push(cutDNA(config.DNA, enzyme));
    });
    electrify(len_list, config.names,
      lane_number=config.lane_number,
      band_width=config.band_width,band_duration=config.duration,
      yScale=config.scale, gel_margin=config.gel_margin, blur=config.band_blur,
      thick_min=config.band_thick_min, thick_rate=config.band_thick_rate,
      tooltip_name_on=config.tooltip_name_on, tooltip_name_size=config.tooltip_name_size,
      tooltip_name_offsetX=config.tooltip_name_offsetX, tooltip_name_offsetY=config.tooltip_name_offsetY,
      tooltip_band_on=config.tooltip_band_on, tooltip_band_size=config.tooltip_band_size,
      tooltip_band_offsetX=config.tooltip_band_offsetX, tooltip_band_offsetY=config.tooltip_band_offsetY);
  }

  // --------------- Method Chain --------------
  // Make black background by .call() before make bands.
  elect.makeGel = function(selection) {
    selection.append("rect")
        .attr("width", "100%")
        .attr("height", "100%")
        .attr("fill", "#222222");
    return elect;
  };
  elect.DNA = function(value) {
    if (!arguments.length) return config.DNA;
    config.DNA = value;
    return elect;
  };
  elect.enzymes = function(value) {
    if (!arguments.length) return config.enzymes;
    config.enzymes = value;
    return elect;
  };
  elect.names = function(value) {
    if (!arguments.length) return config.names;
    config.names = value;
    return elect;
  };
  elect.scale = function(value) {
    if (!arguments.length) return config.scale;
    config.scale = value;
    return elect;
  };
  elect.lane_number = function(value) {
    if (!arguments.length) return config.lane_number;
    config.lane_number = value;
    return elect;
  };
  elect.gel_margin = function(value) {
    if (!arguments.length) return config.gel_margin;
    config.gel_margin = value;
    return elect;
  };
  elect.duration = function(value) {
    if (!arguments.length) return config.duration;
    config.duration = value;
    return elect;
  };
  elect.band_width = function(value) {
    if (!arguments.length) return config.band_width;
    config.band_width = value;
    return elect;
  };
  elect.band_blur = function(value) {
    if (!arguments.length) return config.band_blur;
    config.band_blur = value;
    return elect;
  };
  elect.band_thick_min = function(value) {
    if (!arguments.length) return config.band_thick_min;
    config.band_thick_min = value;
    return elect;
  };
  elect.band_thick_rate = function(value) {
    if (!arguments.length) return config.band_thick_rate;
    config.band_thick_rate = value;
    return elect;
  };
  elect.tooltip_name_on = function(value) {
    if (!arguments.length) return config.tooltip_name_on;
    config.tooltip_name_on = value;
    return elect;
  };
  elect.tooltip_name_size = function(value) {
    if (!arguments.length) return config.tooltip_name_size;
    config.tooltip_name_size = value;
    return elect;
  };
  elect.tooltip_name_offsetX = function(value) {
    if (!arguments.length) return config.tooltip_name_offsetX;
    config.tooltip_name_offsetX = value;
    return elect;
  };
  elect.tooltip_name_offsetY = function(value) {
    if (!arguments.length) return config.tooltip_name_offsetY;
    config.tooltip_name_offsetY = value;
    return elect;
  };
  elect.tooltip_band_on = function(value) {
    if (!arguments.length) return config.tooltip_band_on;
    config.tooltip_band_on = value;
    return elect;
  };
  elect.tooltip_band_size = function(value) {
    if (!arguments.length) return config.tooltip_band_size;
    config.tooltip_band_size = value;
    return elect;
  };
  elect.tooltip_band_offsetX = function(value) {
    if (!arguments.length) return config.tooltip_band_offsetX;
    config.tooltip_band_offsetX = value;
    return elect;
  };
  elect.tooltip_band_offsetY = function(value) {
    if (!arguments.length) return config.tooltip_band_offsetY;
    config.tooltip_band_offsetY = value;
    return elect;
  };

  return elect;
}
