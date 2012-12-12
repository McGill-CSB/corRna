function main() {
	 //applet();
	$(document).ready(function() {
		//ie edition
		$('#body').html('<div id="log"></div><div id="def"><p></p></div><div id="table" class="table"></div><div id="panel"></div>');
		if(navigator.appName == "Microsoft Internet Explorer") {
			var x = output[0].type;
			if(x != null) {
				TOKEN = NODE[0];
				var tmp = output;
				output = [];
				/// FILTER
				for(var i = 0; i<tmp.length;i++) {
					var bad = false;
					for(var j=0;j<tmp[i].seq.length;j++) {
						var y = tmp[i].seq.charAt(j);
						if(y.match(/\./)) { bad = true; break;}
					}
					if(bad == false)
						output.push(tmp[i]);
				}
				wildtype = tinyWild(output); 
				panel();
				print(output);
				highlight();
			//	bind();
			} else if (output[0].method != undefined) {
				if(output[0].method == "mutation") {
					TOKEN = NODE[2];
					mutation_main(output[0]);
				} else {
					TOKEN = NODE[1];
					structural_main(output[0]);
				}
			}
		} else {
			if(output[0].type != undefined) {
				TOKEN = NODE[0];
				var tmp = output;
				output = [];
				/// FILTER
				for(var i = 0; i<tmp.length;i++) {
					var bad = false;
					for(var j=0;j<tmp[i].seq.length;j++) {
						if(tmp[i].seq[j].match(/\./)) { bad = true; break;}
					}
					if(bad == false)
						output.push(tmp[i]);
				}
				wildtype = tinyWild(output); 
			//	$('#body').html('<div id="table" class="table"></div><div id="panel"></div>');
				panel();
				print(output);
				highlight();
			//	bind();
			} else if (output[0].method != undefined) {
				if(output[0].method == "mutation") {
					TOKEN = NODE[2];
					mutation_main(output[0]);
				} else {
					TOKEN = NODE[1];
					structural_main(output[0]);
				}
			}
		}
	});
}
var wildtype = "";
var SIG = [];
var DEBUG = false;
var truth = "";
var TOKEN = "";
var NODE = ["None", "Structural", "Mutational"];

//Structural
function structural_main(output) {
	$('#def p').html('Wildtype: '+output.wild+'<br>Structure: '+output.struct+'<br>');
	var out = output.seq;
	SIG = out;
	wildtype = output.wild;
	panel();
	printStructuralHeader();
	printStructuralBody(out);
}

function printStructuralHeader() {
	var HTML = "<table id='results'><thead><td>Sequence</td><td>Structure</td><td><label id='struct_corr' class='hyper'>Correlation<input type='hidden' value='corr'/></label></td><td><label id='struct_mfe' class='hyper'>MFE<input type='hidden' value='mfe'/></label></td><td><label id='struct_breakN' class='hyper'>BreakNumber Set<input type='hidden' value='breakNumber'/></label></td><td>Structure Set</td><td><label id='struct_boot'class='hyper'>Significant<input type='hidden' value='boot'/></label></td></thead><tbody></tbody></table>";
	$('#table').html(HTML);
	//struct_switch();
	highlight();
}

function printStructuralBody(out) {
	var HTML = "";	
	for( var i=0;i<out.length;i++) {
		var STR = out[i];
		if(i!=0 && i%2 != 0) {
			HTML += "<tr class='alt'>";
		} else {
			HTML += "<tr>";
		}	
		if(navigator.appName != "Microsoft Internet Explorer") {
			HTML+='<td><a href="#'+i+'" onclick='+"'"+'return _applet("'+STR.seq+'","'+STR.struct+'","'+tinySeq(STR.seq)+'")'+"'"+'><div id="'+i+'">'+tinySeq(STR.seq)+'</div></a></td><td>'+STR.struct+'</td><td>'+round(STR.corr)+'</td><td>'+STR.mfe+'</td><td>'+STR.breakNumber+'</td><td>'+STR.breakStruct+'</td><td>'+round(STR.boot)+'</td></tr>';
		} else {
			var str = STR.seq;
			HTML+='<td><a href="#'+i+'" onclick='+"'"+'return _applet("'+STR.seq+'","'+STR.struct+'","'+tinySeq(str)+'")'+"'"+'><div id="'+i+'">'+tinySeq(str)+'</div></a></td><td>'+STR.struct+'</td><td>'+round(STR.corr)+'</td><td>'+STR.mfe+'</td><td>'+STR.breakNumber+'</td><td>'+STR.breakStruct+'</td><td>'+round(STR.boot)+'</td></tr>';
		}
	}
	$('#results tbody').html(HTML);
}

function round(x){
	return Math.round(x*1000000)/1000000
}
function highlight() {
	$('label.hyper').click(function() {
		$('label.hyper').removeClass('activeHyper');
		$(this).addClass('activeHyper');
		val = $(this).parent().find("input").val();
			if(TOKEN == NODE[0]) {
				if(val == truth) { 
					printBody(SIG.sort(sortBy(val,true, parseFloat)));
					truth = "";
				} else { 
					printBody(SIG.sort(sortBy(val,false, parseFloat)));	
					truth = val;
				}

			} else if (TOKEN == NODE[1]) {
				if(val == truth) { 
					printStructuralBody(SIG.sort(sortBy(val,true, parseFloat)));
					truth = "";
				} else { 
					printStructuralBody(SIG.sort(sortBy(val,false, parseFloat)));	
					truth = val;
				}

			} else if (TOKEN == NODE[2]) {
				if(val == truth) { 
					printMutationalBody(SIG.sort(sortBy(val,true, parseFloat)));
					truth = "";
				} else { 
					printMutationalBody(SIG.sort(sortBy(val,false, parseFloat)));	
					truth = val;
				}

			}
	});
}
var sortBy = function(field, reverse, primer) {
	reverse = (reverse)? -1: 1;
	return function (a,b) {
		a = a[field];
		b = b[field];
		if(typeof(primer) != 'undefined') {
			a = primer(a);
			b = primer(b);
		}
		if(a<b) return reverse * -1;
		if(a>b) return reverse * 1;
		return 0;
	}
}
///

/// MUTATION
function mutation_main(output) {
	var none = false;
	$('#def p').html('<div id="mutation" class="mutation"></div><div id="mutation_sig"><div id="mutation_sig_header"></div>');		
	var HTML="<p><span class='mutation_settings'>Detection Settings:<br>Wildtype: "+output.wildtype+"&nbsp;&nbsp;<pre style='display:none'>Average: "+output.setAvg+"&nbsp;&nbsp;</pre> Threshold:";
	wildtype = output.wildtype;
	if(output.setThres == "-1.0") {
		HTML+= " None (Unfiltered Results)<br></span>";	
		none = true;
	} else {
		HTML+= " "+output.setThres+"<br></span>";	
	}
	HTML += "<a href='#' id='viewRawData'>Show Raw Data</a><br><div class='viewRaw' style='display:none'>";
	HTML += "&lt;Under Avg>&gt; :: Under "+output.setAvg+" across 13 given confirmed HCV from our DB <br>"; 
	out = output.data;
	//var SIG = [];
	for(var i=0;i<out.length;i++) {		
		seq = out[i].seq;	
		for(var j=0;j<seq.length;j++) {
			var STR = seq[j];
			if(none == true) SIG.push(seq[j]);
			HTML += "<p class='sig'>"+STR.seq+"&nbsp;&nbsp;"+STR.struct+"&nbsp;&nbsp;"+STR.corr+"&nbsp;&nbsp;"+STR.mfe; 
			if(STR.avg == 't') {
				HTML+= " &lt;Under Avg&gt;";
			}
			if(STR.navg == 't') {
				SIG.push(seq[j]);				/////CHECK SYNTAX	
				HTML+= " &lt;Under Threshold&gt;";
			}
			HTML+="</p>";
		}
		HTML+="<p>Mutation Space: "+out[i].node+"&nbsp; Avg Correlation: "+out[i].avg+"&nbsp; Time: "+out[i].time+"</p>";
	}	
	HTML += "</div></p>";
	$('#mutation').html(HTML);
	SIG_head = "<p>Total Seq: "+output.tSeq+"&nbsp; Total Time: "+output.tTime+"&nbsp; Avg Time: "+output.tAvgTime+"&nbsp; Avg Correlation: "+output.tAvgCorr+"<br>Total Significants Found:"+output.tThres+" &nbsp Under Threshold:"+output.setThres+"</p>";
	$('#mutation_sig_head').html(SIG_head);
	viewRawDataHandler();
	panel();
	print_mutation_sig(SIG);
}

function viewRawDataHandler() {
	$('#viewRawData').click(function() {
		if($('#viewRawData').html() == "Show Raw Data") {
			$('.viewRaw').show();
			$('#viewRawData').html('Hide Raw Data');
		} else {
			$('.viewRaw').hide();
			$('#viewRawData').html('Show Raw Data');
		}
	});
}

function print_mutation_sig(SIG) {
	var HTML  = "<table id='results'><thead><tr><td>Sequence</td><td>Structure</td><td><label class='hyper' id='mutation_corr'>Correlation<input type='hidden' value='corr'/></label></td><td><label class='hyper' id='mutation_mfe'>MFE<input type='hidden' value='mfe'/></label></tr></thead><tbody></tbody></table>";
	$("#table").html(HTML);
	printMutationalBody(SIG);
	highlight();
}

function printMutationalBody(SIG) {
	var HTML = "";
	for(var i=0;i<SIG.length;i++) {
		if(i!=0 && i%2 != 0) {
			HTML += "<tr class='alt'>";
		} else {
			HTML += "<tr>";
		}	
		if(navigator.appName != "Microsoft Internet Explorer") {
			HTML+="<td><a href='#"+i+"' onclick="+'"'+"return _applet('"+SIG[i].seq.toUpperCase()+"','"+SIG[i].struct+"','"+tinySeq(SIG[i].seq)+"');"+'"'+"><div id='"+i+"'>"+tinySeq(SIG[i].seq)+"</div></a></td><td>"+SIG[i].struct+"</td><td>"+SIG[i].corr+"</td><td>"+SIG[i].mfe+"</td></tr>";	
		} else {
			var str = SIG[i].seq;
			HTML+="<td><a href='#"+i+"' onclick="+'"'+"return _applet('"+SIG[i].seq.toUpperCase()+"','"+SIG[i].struct+"','"+tinySeq(str)+"');"+'"'+"><div id='"+i+"'>"+tinySeq(str)+"</div></a></td><td>"+SIG[i].struct+"</td><td>"+SIG[i].corr+"</td><td>"+SIG[i].mfe+"</td></tr>";	
		}
	}
	$("#results tbody").html(HTML);
}
//END of mutation

///   NO STRUCT
function tinySeq(seq) {
	if(navigator.appName == "Microsoft Internet Explorer") {
		//$("#log").append(seq+" "+wildtype+"<br>");
		var str = "";
		seq = seq.toUpperCase();
		for(var i=0;i<seq.length;i++) {
			//$("#log").append(seq.charAt(i)+"="+wildtype.charAt(i)+"<>");
			if(wildtype.charAt(i) != seq.charAt(i)) {
				if(str!="") str+="_";
				str+=wildtype.charAt(i)+(i+1)+seq.charAt(i).toUpperCase();	
			}
		}
		//$("#log").append(str+"<br>");
		if (str=="") str = "wildtype";
		return str;
	
	} else {
		var str = "";
		seq = seq.toUpperCase();
		for(var i=0;i<seq.length;i++) {
			if(wildtype[i] != seq[i]) {
				if(str!="") str+="_";
				str+=wildtype[i]+(i+1)+seq[i].toUpperCase();	
			}
		}
		if (str=="") str = "wildtype";
		return str;
	}
}

function tinyWild(out) {
	if(out[0].type != 'R' && out[0].mutation != "0") {
		out.sort(sort_corr_asc);
	}
	return out[0].seq;
}

function print(out) {
	var HTML = "Method Definition F: fixedCGSampling R: RNAmutant<br>Wildtype: "+wildtype;
	
	$("#def p").html(HTML);

	HTML = "<table id='results'><thead><tr><td class='type'>Method</td><td>Sequence</td><td>Secondary Structure</td><td><label class='hyper' id='mutation'>No of Mutations<input type='hidden' value='mutation'/></label></td><td><label class='hyper' id='MFE'>MFE<input type='hidden' value='MFE'/></label></td><td><label class='hyper' id='corr_l'>Correlation<input type='hidden' value='corr'/></label></td><td><label id='boot' class='hyper'>Significant<input type='hidden' value='boot'/></label></td></tr></thead><tbody></tbody></table>";
	$('#table').html(HTML);
	SIG = out;
	printBody(out);
}
	
function printBody(SIG) {
	var HTML = "";
	for(var i=0;i<SIG.length;i++) {
		var j = SIG[i];
		if(i!=0 && i%2 != 0) {
			HTML += "<tr class='alt'>";
		} else {
			HTML += "<tr>";
		}	
		if(navigator.appName != "Microsoft Internet Explorer") {
			HTML+="<td class='type'>"+j.type+"</td><td><a href='#"+i+"' onclick='return _applet("+'"'+j.seq+'","'+j.struct+'","'+tinySeq(j.seq)+'"'+");'><div id='"+i+"'>"+tinySeq(j.seq)+"</div></a></td><td>"+j.struct+"</td><td>"+j.mutation+"</td><td>"+j.MFE+"</td><td>"+j.corr+"</td><td>"+j.boot+"</td></tr>";
		} else {
			var str = j.seq;
			HTML+="<td class='type'>"+j.type+"</td><td><a href='#"+i+"' onclick='return _applet("+'"'+j.seq+'","'+j.struct+'","'+tinySeq(str)+'"'+");'><div id='"+i+"'>"+tinySeq(str)+"</div></a></td><td>"+j.struct+"</td><td>"+j.mutation+"</td><td>"+j.MFE+"</td><td>"+j.corr+"</td><td>"+j.boot+"</td></tr>";
		}
	}
	$("#results tbody").html(HTML);	
	/*
	var top = $('table thead').offset().top - parseFloat($('table thead').css('marginTop').replace(/auto/,0));
	$(window).scroll(function (event) {
		var y = $(this).scrollTop();
			
			//YUI().use('console', function(Y) {
			//	var Y = new YUI({debug : true});
			//	Y.log(y);
			if(y >= top) {
		//		Y.log("added");
				$('table thead').addClass('fixed');
			} else {
				$('table thead').removeClass('fixed');
		//		Y.log("removed");
			}
		//	});
	});
	*/
/*
	$('.jquery_columnSizing').click(function() {
		 var y = $(this).css("td");
			YUI().use('console', function(Y) {
				var Y = new YUI({debug : true});
				Y.log(y);
			});
	}); */
}
var counter = 0;
function _applet(seq, struct, tiny) {
	$('#app').append('<span class="applet"><label id="remove" class="remove">X</label><applet code="VARNA.class" codebase="../../" archive="VARNA.jar" width="250" height="250"> <param name="sequenceDBN" value="'+seq+'" />' +
	' <param name="structureDBN" value="'+struct+'" /> <param name="title" value="'+tiny+'" /></applet></span>');
	counter+=1;
	$('.remove').click(function() {
		$(this).parent().html('');	
		counter-=1;
		if(counter == 0) {
			$('#panel').hide();
		}
	});
	if($('#app').html() != '') {
		$('#panel').show();
	}
}
/* <---function depreciate--->
function applet() {
	$('#app').html('<applet code="VARNA.class" codebase="../../" archive="VARNA.jar" width="250" height="250"> <param name="sequenceDBN" value="GGGGCCAAUAUGGCCAUCC" />' +
' <param name="structureDBN" value="((((((.....))))..))" /> <param name="title" value="" /></applet>');
}
*/
function activePanel() {
	if($('#hide').html() == 'Hide') {
		$('#app').hide();
		$('#controlPanel').addClass('nonActive');
		$('#hide').html('Show');
	} else {
		$('#app').show();
		$('#controlPanel').removeClass('nonActive');
		$('#hide').html('Hide');
	}
}
function panel() {
	$("#panel").html("<div id='controlPanel'><span id='removeAll'>Clear All</span>&nbsp;&nbsp;<span id='hide' onClick='return activePanel();'>Hide</div></div><div id='app' class='app'></div>")
	$('#panel').hide();
	$('#removeAll').click(function() {
		$('#app').html('');
		$('#hide').html('Hide');
		$('#panel').hide();
		
	});
}

function interactAPI() {
	var head = document.getElementsByTagName("head")[0];
	
	var script = document.createElement('script');
	script.id = 'yui';
	script.type = 'text/javascript';
	script.src = 'http://yui.yahooapis.com/3.2.0/build/yui/yui-min.js';
	head.appendChild(script);
	/*
	var script = document.createElement('script');
	script.id = 'dimension';
	script.type = 'text/javascript';
	script.src = "../../API/jquery.dimensions.pack.js";
	head.appendChild(script);	
	var script = document.createElement('script');
	script.type = 'text/javascript';
	script.id = 'cookies';
	script.src = "../../API/jquery.cookies.pack.js";
	head.appendChild(script);
	var script = document.createElement('script');
	script.type = 'text/javascript';
	script.id = 'iutil';
	script.src = "../../API/jquery.iutil.pack.js";
	head.appendChild(script);
	var script = document.createElement('script');
	script.type = 'text/javascript';
	script.id = 'idrag';
	script.src = "../../API/jquery.idrag.js";
	head.appendChild(script);
	var script = document.createElement('script');
	script.type = 'text/javascript';
	script.id = 'columnSizing';
	script.src = "../../API/jquery.grid.columnSizing.js";
	head.appendChild(script);
	*/
	var css = document.createElement('link');
	css.type= 'text/css';
	css.id = 'css-non-ie';
	css.href = '../../appStyle-non.css';
	css.rel ='stylesheet';
	head.appendChild(css);
	
	
}

function resizeColumn() {
	//reference http://www.ita.es/jquery/jquery.grid.columnSizing.htm
	if(navigator.appName != "Microsoft Internet Explorer") {
		eval("var script = document.getElementById('columnSizing')");
		if(script==null) {
			//alert('test');
		} else 
			$("table").columnSizing({
				selectCells : "tr:first>*:not(:first)"
			});
	}
	
}
///// END OF STRUCT / NO STRUCT

interactAPI();
main();
