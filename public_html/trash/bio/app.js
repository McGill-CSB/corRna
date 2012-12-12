function main() {
	 //applet();
	$(document).ready(function() {
		$('#body').html('<div id="table" class="table"></div><div id="panel"></div>');
		panel();
		print(output);
		bind()
	});
 }

function bind() {
	$('#corr_l').click(function() {
		sort(1);
	});
	$('#mutation').click(function(){
		sort(2);
	});
	$('#MFE').click(function(){
		sort(3);
	});
	$('#boot').click(function(){
		sort(4);
	});
	resizeColumn();	
}
function sort(n) {
	switch(n) {
		case 1: //correlation
			if($('#corr_h').val() == '1') {
				print(output.sort(sort_corr_desc));
				$('#corr_h').val('2');
			} else if($('#corr_h').val() == '2') {
				print(output.sort(sort_corr_asc));
				$('#corr_h').val('1');
			} else {
				print(output.sort(sort_corr_desc));
				$('#corr_h').val('2');
			}
			break;
		case 2: //mutation
			if($('#mutation_h').val() == '1') {
				print(output.sort(sort_mutation_desc));
				$('#mutation_h').val('2');
			} else if($('#mutation_h').val() == '2') {
				print(output.sort(sort_mutation_asc));
				$('#mutation_h').val('1');
			} else {
				print(output.sort(sort_mutation_desc));
				$('#mutation_h').val('2');
			}
			break;
		case 3: //MFE
			if($('#MFE_h').val() == '1') {
				print(output.sort(sort_MFE_desc));
				$('#MFE_h').val('2');
			} else if($('#MFE_h').val() == '2') {
				print(output.sort(sort_MFE_asc));
				$('#MFE_h').val('1');
			} else {
				print(output.sort(sort_MFE_desc));
				$('#MFE_h').val('2');
			}
			break;
		case 4: //Sig
			if($('#boot_h').val() == '1') {
				print(output.sort(sort_Significant_desc));
				$('#boot_h').val('2');
			} else if($('#boot_h').val() == '2') {
				print(output.sort(sort_Significant_asc));
				$('#boot_h').val('1');
			} else {
				print(output.sort(sort_Significant_desc));
				$('#boot_h').val('2');
			}
			break;
		default:
			break;
	}
	bind();
}
function sort_Significant_asc(a,b){
	return parseFloat(a.boot) - parseFloat(b.boot);
}
function sort_Significant_desc(a,b) {
	return parseFloat(b.boot) - parseFloat(a.boot);
}
function sort_MFE_asc(a,b){
	return parseFloat(a.MFE) - parseFloat(b.MFE);
}
function sort_MFE_desc(a,b) {
	return parseFloat(b.MFE) - parseFloat(a.MFE);
}
function sort_mutation_asc(a,b){
	return parseInt(a.mutation) - parseInt(b.mutation);
}
function sort_mutation_desc(a,b) {
	return parseInt(b.mutation) - parseInt(a.mutation);
}
function sort_corr_asc(a,b){
	return parseFloat(a.corr) - parseFloat(b.corr);
}
function sort_corr_desc(a,b) {
	return parseFloat(b.corr) - parseFloat(a.corr);
}
function print(out) {
	HTML = "<p class='def'>Method Definition F: fixedCGSampling R: RNAmutant</p><div class='table_wrapper'<table><thead><tr><td class='type'>Method</td><td>Sequence</td><td>Secondary Structure</td><td><label class='hyper' id='mutation'>No of Mutations</label><input type='hidden' id='mutation_h' value='0'/></td><td><label class='hyper' id='MFE'>MFE</label><input type='hidden' value='0'/></td><td><label class='hyper' id='corr_l'>Correlation</label><input type='hidden' id='corr_h' value='0'/></td><td><label id='boot' class='hyper'>Significant</label><input type='hidden' id='boot_h' value='0'/></td></tr></thead><tbody>";
	for(var i=0;i<out.length;i++) {
		var j = out[i];
		if(i!=0 && i%2 != 0) {
			HTML += "<tr class='alt'>";
		} else {
			HTML += "<tr>";
		}	
		HTML+="<td class='type'>"+j.type+"</td><td><a href='#"+i+"' onclick='return _applet("+'"'+j.seq+'","'+j.struct+'"'+");'><div id='"+i+"'>"+j.seq+"</div></a></td><td>"+j.struct+"</td><td>"+j.mutation+"</td><td>"+j.MFE+"</td><td>"+j.corr+"</td><td>"+j.boot+"</td></tr>";
	}
	HTML+="</tbody></table></div>";
$('#table').html(HTML);
	var top = $('table thead').offset().top - parseFloat($('table thead').css('marginTop').replace(/auto/,0));
	$(window).scroll(function (event) {
		var y = $(this).scrollTop();
			/*
			YUI().use('console', function(Y) {
				var Y = new YUI({debug : true});
				Y.log(y);
			*/
			if(y >= top) {
		//		Y.log("added");
				$('table thead').addClass('fixed');
			} else {
				$('table thead').removeClass('fixed');
		//		Y.log("removed");
			}
		//	});
	});
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
function _applet(seq, struct) {
	$('#app').append('<span class="applet"><label id="remove" class="remove">X</label><applet code="VARNA.class" codebase="../../" archive="VARNA.jar" width="250" height="250"> <param name="sequenceDBN" value="'+seq+'" />' +
	' <param name="structureDBN" value="'+struct+'" /> <param name="title" value="'+seq+'" /></applet></span>');
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
	
}
function resizeColumn() {
	//reference http://www.ita.es/jquery/jquery.grid.columnSizing.htm
	eval("var script = document.getElementById('columnSizing')");
	if(script==null)
		alert('test');
	else
		$("table").columnSizing({
			selectCells : "tr:first>*:not(:first)"
		});
	
}

interactAPI();
main();
