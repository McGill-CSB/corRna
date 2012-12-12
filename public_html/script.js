function isValid() {
	if(!validate_seq()) {
		return false;
	}
	var email = $('#email').val();
	email = jQuery.trim(email);
	var seq = $('#sequence').val();
	seq = jQuery.trim(seq);
	var bootstrap = "";
	var data = "e="+email+"&s="+seq
	if($('#Mutation_Heuristic').is(':checked')) {
		var cycle = $.trim($('#cycle').val());
		var threshold = $.trim($('#threshold').val());
		var mutation = $.trim($("#mDepth").val());
		var mtop = $.trim($('#top').val());
		if(cycle == "") {
			return false;
		} else {
			data+="&cycle="+cycle;
		}
		if(threshold !="") {
			data+="&threshold="+threshold;	
		}
		data+="&meth=RNAmutant";
		if(mutation == "") {
			return false;
		} else {
			data+="&mutation="+mutation;
		}
		if(mtop == "") {
			return false;
		} else {
			data+="&t="+mtop;	
		}
		data+="&constraints=mutation";
		
	} else if($("#Structure_Heuristic").is(':checked')) {
		var mutation = $.trim($("#mDepth").val());
		var breakN = $.trim($('#breakN').val()); 
		if(mutation == "") {
			return false;
		} else {
			data+="&mutation="+mutation;
		}
		if(breakN == "") {
			return false;
		} else {
			data+="&breakN"+breakN
		}
		data+="&constraints=structural";
		
	}else {
		if($('#bootstrap').is(':checked')) {
			bootstrap = 't';
		} else {
			bootstrap = 'f';
		}
		var meth = "";
		if($('#RNAmutant').is(':checked')) {
			meth = "RNAmutant";
			
		} else if($('#fixedCGSampling').is(':checked')) {
			meth = "fixedCGSampling";
		} else {
			meth = "RNAfixed";
		}
		var mutation = $("#mDepth").val();
		var dangling = "";
		if($('#d0').is(':checked')) {
			dangling = "d0";
		} else if($('#d1').is(':checked')) {
			dangling = "d1";
		} else if($('#d2').is(':checked')) {
			dangling = "d2";
		} else if($('#d3').is('checked')) {
			dangling = "d3";
		}
		var eValue = $('#e').val();
		var mtop = $('#top').val();
		data+="&constraints=none&b="+bootstrap+"&t="+mtop+"&meth="+meth+"&eValue="+eValue+"&dangling="+dangling+"&mutation="+mutation;
	}
	$('#wemail').html('');
	$('#wseq').html('');
	$('#note').html('Loading');
	$('#results').html('');
	//alert(data);
	$.ajax({
		type: "POST",
		url: "Submit.php",
		data: data,
		success: function(data) {
		//$('#note').html("Please check ur email, for data will be sent to");	
		$('#note').html(data);
	}});
	return false;
}

function procedure(op) {
	if(op=='RNAmutant') {
		$('.extra_settings').hide();
	} else {
		$('.extra_settings').show();
	}
}

function validate_e() {
	value = $('#id').val();
	value = $.trim(value);
	if(!isNumeric(value)) {
		$('#we').html("Numbers only!");			
		return false;
	} else {
		$('#we').html("");
		return true;
	}
}

function validate_top() {
	value = $('#top').val();
	value = $.trim(value);
	
	//if(value != 'all') {
		if(!isNumeric(inputVal)) {
			$('#wtops').html("Numbers only!");
			return false;
		} else {
			$('#wtop').html("");
			return true;
		}
	/*} else {
		$('#wtop').html("");
		return true;
	}*/
}

function validate_mDepth() {
	value = $('#mDepth').val();
	value = $.trim(value);
		if(!isNumeric(inputVal)) {
			$('#wmDepth').html('Numbers only!');
			return false;
		} else {
			$('wmDepth').html("");
			return true;
		}
}

function validate_seq() {
	raw_value = $('#sequence').val();
	str = $.trim(raw_value);
	sub = str.split("\n");
	$('#wseq').html('');
//	if($('#Structure_Heuristic').is(':checked') == false) {
		if(sub.length == 1) {
			var str = sub[0];
			for(var i=0;i<str.length;i++) {
				if(!str[i].match(/A|U|C|G|a|u|c|G/)) {
					$('#wseq').html('Please enter an actual sequence');
					return false;
				}
			}
			$('#wseq').html();
			return true;
		} else {
			$('#wseq').html('Incorrect Entry');
			return false;
		}
	/*
	} else {
		if(sub.length == 3) {
			bSeq = false;
			bMarker1 = false;
			bMarker2 = false;
			if(sub.length != 3) {
				$('#wseq').html('Incorrect Entry');
				return false;
			} else if(sub[0].length != sub[1].length || sub[0].length != sub[2].length || sub[1].length != sub[2].length) {
				if(sub[0].length != sub[1].length && sub[0].length != sub[2].length && sub[1].length != sub[2].length) {
					$('#wseq').html('Neither of the lines\' length match');
					return false;
				} else if(sub[0].length == sub[1].length && sub[2].length != sub[0].length) {
					$('#wseq').html('Line 3\'s length is the odd one out');
					return false;
				} else if(sub[0].length == sub[2].length && sub[1].length != sub[0].length) {
					$('#wseq').html('Line 2\'s length is the odd one out');
					return false;
				} else if(sub[1].length == sub[2].length && sub[2].length != sub[0].length) {
					$('#wseq').html('Line 1\'s length is the odd one out');
					return false;
				}
			}
			for(var i=0;i<sub.length;i++) {
				var str = sub[i];
				if(str.match(/A|U|C|G|a|u|c|G/) && bSeq == false) {
					for(var j=0;j<str.length;j++) {
						if(!str[j].match(/A|U|C|G|a|u|c|G/)) {
							$('#wseq').append('There is a typo on the seqeunce!');
							return false;
						}					
					}
					bSeq = true;
				} else if(str.match(/o|X/) && bMarker1 == false) {
					for(var j=0;j<str.length;j++) {
						if(!str[j].match(/o|X/)) {
							$('#wseq').append('There is a typo on the o / X marker');
							return false;
						}
					}
					bMarker1 = true;
				} else if(bMarker2 == false && str.match(/\.|\?/)) {
					for(var j=0;j<str.length;j++) {
						if(!str[j].match(/\.|\?/)) {
							$('#wseq').append('There is a typo on the . / ? marker');
							return false;
						}
					}
					bMarker2 = true;
				} else {
					$('#wseq').html('Incorrect Entry');
					return false;
				}
			}
			return true;
		}
	}*/
	return false;
}

function remove_str() {
	raw_str = $('#sequence').val();
	str = $.trim(raw_str); 
	if(str == 'Enter sequence and/or constrains here') {
		$('#sequence').val('');
	}
}

function wseq(str){
	$('#wseq').html(str);
}

function isNumeric(inputVal) {
     if (isNaN(parseFloat(inputVal))) {
          return false;
     }
     return true
}

function main() {	

	helpContent = "Where the manual goes~";
	HTML = '<ul class="tabs">'+
		'<li><a href="#compute">Compute</a></li>'+
		'<li><a href="#intro">Intro</a></li>'+
		'<li><a href="#help">Manual</a><li>'+
		'<li><a href="#sample">Example</a><li>'+
		'<li><a href="#contacts">Contacts</a></li>'+
		'</ul>'+
		'<div class="tab_container">'+
			'<div id="compute" class="tab_content">'+
			'</div>'+
			'<div id="intro" class="tab_content">' +
			'</div>' +
			'<div id="help" class="tab_content">' +
			'</div>' +
			'<div id="sample" class="tab_content">'+
			'</div>' +
			'<div id="contacts" class="tab_content">'+
			'</div>' +
		'</div>';
	$(document).ready(function() {
		$("#header").html("<h1>corRna: Predict and analyze deleterious mutations</h1>")
		$('#body').html(HTML);

		if(navigator.appName == "Microsoft Internet Explorer") {
			//document.write("*Note: Sorry for the inconvinence, Currently Site Maintanence");
			YUI().use('node', function(Y) {
				Y.all(".tab_content").setStyle('display', 'none');
				Y.one("ul.tabs li").addClass("active").setStyle('display', null);
				Y.one(".tab_content").setStyle('display', null);
				Y.all("ul.tabs li").each(function(node) {
						node.on("click",function(x) {
						Y.all("ul.tabs li").removeClass("active");
						this.addClass("active");
						Y.all(".tab_content").setStyle('display', 'none');
						var activeTab = this.one("a").get("href");
						activeTab = activeTab.replace(/^.*#/,"");
						Y.one("#"+activeTab).setStyle('display',null);
						return false;
					});
				});
			});	
			
		} else {
		//When page loads...
			$(".tab_content").hide(); //Hide all content
			$("ul.tabs li:first").addClass("active").show(); //Activate first tab
			$(".tab_content:first").show(); //Show first tab content

			//On Click Event
			$("ul.tabs li").click(function() {

				$("ul.tabs li").removeClass("active"); //Remove any "active" class
				$(this).addClass("active"); //Add "active" class to selected tab
				$(".tab_content").hide(); //Hide all tab content

				var activeTab = $(this).find("a").attr("href"); //Find the href attribute value to identify the active tab + content
				$(activeTab).fadeIn(); //Fade in the active ID content
				return false;
			});
		} 
	
		setContent("Content/sampleContent.html","sample");
		setContent("Content/compute.html","compute");
		setContent("Content/manual.html","help");
		setContent("Content/intro.html", "intro");
		setContent("Content/contacts.html", "contacts");

	});
}
function setContent(loc, append) {
	$.ajax({
		type: 'post',
		url: loc,
		success: function (data) {
			$("#"+append).html(data);
		}
	});
}
function methology(op) {
	if(op == 'structure' || op == 'mutation') {
		$('.methology_active').hide();
		$('.extra_settings').hide();
		$('#RNAmutant').attr("checked", "checked");	
	} else {
		$('.methology_active').show();
	}

	if(op=='mutation') {
		$('.Mutational').show();
		$('.Structural').hide();
		$('.Struct').hide();
		$('.nnn').show();
	} else if(op=='structure') {
		$('.Mutational').hide();
		$('.Structural').hide();
		$('.Struct').show();
		$('.nnn').hide();
	} else {
		$('.Mutational').hide();
		$('.Structural').show();
		$('.Struct').hide();
		$('.nnn').show();
	}
}

main();
