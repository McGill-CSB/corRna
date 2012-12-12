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
			//we accepted a sequence
			$seq = $_POST['s'];	//have to sort
			$email = $_POST['e'];
			$eValue = $_POST['eValue'];	///numeric
			$method = $_POST['meth'];
			$bootstrap = $_POST['b']; // f | t
			$nSize = $_POST['t'];	// all | numeric
			$mutation = $_POST['mutation'];
			$dangling = $_POST['dangling'];

			$loc = getDirectory();

			$param="./script.py -j -s ".$seq;
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
			$res = exec("./client.py $loc $param");
			echo "<b> <a href='results/".$loc."'>To View Status of result, Please click here</a><b>";
			if($DEBUG) {
				echo "<br> param:: $param";
			}
			
			
		} 	
	}
?>
