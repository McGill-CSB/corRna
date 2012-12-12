<?php
ini_set('display_errors', '1');

$queue = new jobQueue();
$queue->active();

class Process {
	private $seq, $email, $error, $method, $boot, $nSize;
	function __construct($SEQ, $EMAIL, $ERROR, $METHOD, $BOOT, $NSIZE) {
		$this->seq = $SEQ;
		$this->email = $EMAIL;
		$this->error = $ERROR;
		$this->method = $METHOD;
		$this->boot = $BOOT;
		$this->nSize = $NSIZE;
	}

	function start() {
		//set parameters and start
		if($method == "RNAmutant") { 

		} else if($method == "fixedCGSampling") {
			//fixedCG property eValue

		} else {  //run both RNAmutant and fixedCG

		}
		return $out;
	}
}
class jobQueue {
	private $jobs = array();
	private $page = array();
	public function __construct() {
		//initialize	
		$this->monitor();
	}	
	public function active() {

	}
	public function qSize() {
		return count($this->jobs);
	}
	public function add(Process $info) {
		//make page
		array_unshift($this->jobs, $info);
	} 
	private function monitor() {
		while(1){
			if($this->ServerLoad() < 300) {
				if(count($this->jobs) > 0) {
					$current = array_pop($this->jobs);	
					$page = array_pop($this->pages);
					$out = $current->start();
					$this->pHTML($out, $page);
				}
			}
		}
	}
	private function pHTML($out, $page) {
		
	}
	private function ServerLoad() {
		$stats = exec("ps aux|awk 'NR > 0 { s +=$3 }; END {print s}'");
		return $stats;
	}
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
