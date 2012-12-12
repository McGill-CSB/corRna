<?php
	$DEBUG = false;
	if($DEBUG) {
		ini_set('display_errors', '1');
	}
	function genRandomString() {
	    $length = 10;
	    $characters = "0123456789abcdefghijklmnopqrstuvwxyz";
	    $string = "";    

	    for ($p = 0; $p < $length; $p++) {
		$string .= $characters[mt_rand(0, strlen($characters))];
	    }

	    return $string;
	}
	
	function getDirectory() {
		$loc = genRandomString();
		if(is_dir($loc)) {
			return getDirectory();
		}
		if(!mkdir("results/$loc",0777)) {
			die('Failed to create folder...');
		}
		$file = fopen("./results/$loc/index.php","w")
			or die("<br>Cant open file<br>");
		$pass = fwrite($file, "Waiting in Queue");
		if(!$pass) {
			echo "Failed to write <br>";
		}   
		fclose($file);
		chmod("./results/$loc/index.php", 0777);
		return $loc;
	}

	if(!empty($_POST['s'])) {
		if($_POST['s']!= "") {
			$email = $_POST['e'];
			$seq = $_POST['s'];	//have to sort
			$param = "";
			if($_POST['constraints'] == "none") {
				//we accepted a sequence
				///	if detect structural mutation....

				///	
				$eValue = $_POST['eValue'];	///numeric
				$method = $_POST['meth'];
				$bootstrap = $_POST['b']; // f | t
				$nSize = $_POST['t'];	// all | numeric
				$mutation = $_POST['mutation'];
				$dangling = $_POST['dangling'];


				$param+="./script.py -j -s ".$seq;
				if($method == "fixedCGSampling") {
					$param.= ' -g';
					if($eValue!= '')
						$param.= ' -e '.$eValue;
				} else if($method == "RNAmutant") {
					$param.= ' -r';
				} else {
				$param.= ' -r -g';
					if($eValue!= '')
						$param.= ' -e '.$eValue;
				}
				if($bootstrap == 't')
					$param.= ' -b';
				if($nSize != '') 
					$param.= ' -n '.$nSize;
				if($mutation != '') 
					$param.=' -m '.$mutation;
				if($dangling != '')
					$param.=' -d '.$dangling;
				if($email != "")
					$param.=" -a $email";
				if($DEBUG) {
					echo "<br> param:: $param";
				} else {
					$loc = getDirectory();
					$param.= " -u http://corrna.cs.mcgill.ca/results/$loc";
					echo "<b> <a href='results/".$loc."'>To View Status of result, Please click here</a><b>";
					$res = exec("./client.py $loc $param");
				}
			} else if($_POST['constraints'] == "mutation") {
				$cycle = $_POST['cycle'];
				$threshold = $_POST['threshold'];		
				$mutation = $_POST['mutation'];
				$top = $_POST['t'];
				$param+="./wrapper.py -s ".$seq;
				if($mutation == "" || $mutation < 1) {
					$param.=" -m 2";	
				} else {
					$param.=" -m ".$mutation;	
				}
				if($cycle == "" || $cycle < 1) {
					$param.=" -l 10";
				} else {
					$param.=" -l ".$cycle;
				}
				if($top == "" || $top < 1) {
					$param.=" -n 10";
				} else {
					$param.=" -n ".$top;
				}
				if($threshold == "" || $threshold < 0) {
					$param.=" -t -1";
				} else {
					$param.=" -t ".$threshold;
				}
				if($email != "")
					$param.=" -a $email";
				if($DEBUG) {
					echo "<br> param :: $param";
				} else {
					$loc = getDirectory();
					$param.= " -u http://corrna.cs.mcgill.ca/results/$loc";
					echo "<b> <a href='results/".$loc."'>To View Status of result, Please click here</a><b>";
					$res = exec("./client.py $loc $param");
				}
			}
		} 	
	}
?>
