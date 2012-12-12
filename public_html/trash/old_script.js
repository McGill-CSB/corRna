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
	if($('#Mutation_Heurisitc').is(':checked')) {
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
		
	} else {
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
		data+="&b="+bootstrap+"&t="+mtop+"&meth="+meth+"&eValue="+eValue+"&dangling="+dangling+"&mutation="+mutation;
	}
	$('#wemail').html('');
	$('#wseq').html('');
	$('#note').html('Loading');
	$('#results').html('');
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
	} else if(sub.length == 3) {
		return false;
		for(var i=0;i<sub.length;i++) {
			var str = sub[i];
			for(var j=0;j<str.length;j++) {
			
			}
		}
	} else {
		return false;
	}
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
		'<li><a href="#sample">Sample</a><li>'+
		'<li><a href="#bugs">Bugs</a></li>'+
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
			'<div id="bugs" class="tab_content">'+
			'</div>' +
		'</div>';
	$(document).ready(function() {
		$("#header").html("<h1>Cor-RNA: Predict and analyze deleterious mutations</h1>")
		$('#body').html(HTML);

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
	
		setContent("Content/sampleContent.html","sample");
		setContent("Content/compute.html","compute");
		setContent("Content/Manual.html","help");
		setContent("Content/intro.html", "intro");
		setContent("Content/bugs.html", "bugs");

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
	} else {
		$('.Mutational').hide();
		$('.Structural').show();
	}
}

main();
