<?php
	if($_POST['s']!= "") {
		//we accepted a sequence
		$seq = $_POST['s'];	//have to sort
		$email = $_POST['e'];
		$eValue = $_POST['eValue'];	///numeric
		$method = $_POST['meth'];
		$bootstrap = $_POST['b']; // f | t
		$nSize = $_POST['t'];	// all | numeric
	
		if($method == "RNAmutant") { 

		} else if($method == "fixedCGSampling") {
			//fixedCG property eValue

		} else {  //run both RNAmutant and fixedCG

		}

		echo "Processing\n";
	}	

class timer {
	public $time = 0, $start, $stop;
	public 	function mtime() {
		$mtime = microtime();
		$mtime = explode(' ', $mtime);
		$mtime = $mtime[1]+$mtime[0];
		return  $mtime;
	}

	public function start() {
		$this->start = $this->mtime();
	}
	
	public function stop() {
		$this->stop = $this->mtime();
		$this->time = $this->stop - $this->start; 
	}
}

?>
